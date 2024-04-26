import argparse
import logging
import pathlib
import random
import warnings
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from typeguard import typechecked

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s"
)
sc.settings.verbosity = 2
sc.logging.print_header()


# https://stackoverflow.com/questions/55063560/replace-empty-strings-on-argparse-to-none
def nullable_string(val):
    if not val:
        return None
    return val


def parse_bool(val):
    return val.lower() == "true"


@typechecked
def split_data_by_category(data: AnnData, category: str):
    category_counts = data.obs[category].value_counts()
    category_data_dict = {
        c: data[data.obs[category] == c] for c in category_counts.index
    }
    return category_counts, category_data_dict


@typechecked
def preprocess_batch(
    batch: AnnData,
    n_genes_by_counts_upper: int,
    pct_counts_mt_upper: float,
    filter_hvg: bool,
    filter_hvg_n: int,
) -> AnnData:
    sc.pp.calculate_qc_metrics(
        batch, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    batch = batch[
        (batch.obs.n_genes_by_counts < n_genes_by_counts_upper)
        & (batch.obs.pct_counts_mt < pct_counts_mt_upper)
    ]
    if filter_hvg:
        sc.pp.highly_variable_genes(
            adata=batch,
            n_top_genes=filter_hvg_n,
        )
        batch = batch[:, batch.var.highly_variable]
    return batch


@typechecked
def main(
    read_path: str,
    save_path: str,
    organism: str,
    gene_id_column: Optional[str],
    batch_column: Optional[str],
    min_genes: int,
    min_cells: int,
    n_genes_by_counts_upper: int,
    pct_counts_mt_upper: float,
    cell_cycle_genes_reference: str,
    filter_hvg: bool,
    filter_hvg_n: int,
    target_column: Optional[str],
    use_gene_id: bool,
    na_string: Optional[str] = "NA",
    doublet_column: Optional[str] = None,
    doublet_value: Optional[str] = None,
    use_index_barcode: bool = True,
    barcode_column: Optional[str] = None,
    attrname: Optional[str] = None,
) -> int:
    data = ad.read_h5ad(read_path)
    if doublet_column is not None and doublet_value is not None:
        data = data[data.obs[doublet_column].values != doublet_value, :]
    if use_index_barcode:
        data.obs["cell_id"] = data.obs.index
        data.obs.reset_index(drop=True, inplace=True)
    elif barcode_column != "cell_id":
        data.obs["cell_id"] = data.obs[barcode_column].values
        data.obs.drop(labels=barcode_column, axis=1, inplace=True)
    if "counts" not in data.layers.keys():
        data.layers["counts"] = data.X
    logging.info(f"data.shape: {data.shape}")
    sc.pp.filter_cells(data, min_genes=min_genes)
    if target_column is not None and target_column != "y":
        data.obs["y"] = data.obs[target_column].values
        data.obs.drop(labels=target_column, axis=1, inplace=True)
    if target_column is not None and na_string is not None:
        data = data[data.obs["y"].values != na_string, :]

    if batch_column is None:
        warnings.warn(
            "`batch_column` was passed as `None`. Assuming that all samples belong to the same batch.",
            UserWarning,
        )
        random.seed(read_path + save_path + organism)
        batch_id = random.randrange(int(1e10))
        data.obs["batch_id"] = np.array(
            [batch_id for i in range(data.obs.shape[0])]
        ).astype(str)
    elif batch_column not in data.obs.columns:
        warnings.warn(
            "`batch_column` was passed as `None`. Assuming that all samples belong to the same batch.",
            UserWarning,
        )
        random.seed(read_path + save_path + organism)
        batch_id = random.randrange(int(1e32))
        data.obs["batch_id"] = np.array([batch_id for i in range(data.obs.shape[0])])
    elif batch_column != "batch_id":

        data.obs["batch_id"] = data.obs[batch_column].values
        data.obs.drop(batch_column, axis=1, inplace=True)
    # Keep genes that appear in more than min_cells cells in every batch
    batch_dict = {i: data[data.obs.batch_id == i] for i in data.obs.batch_id.unique()}
    filtered_genes = []
    for _, v in batch_dict.items():
        l = sc.pp.filter_genes(v, min_cells=min_cells, inplace=False)[0]
        filtered_genes.append(l)
    filtered_genes_all = np.all(filtered_genes, axis=0)
    logging.info(
        f"Filtered out {len(filtered_genes_all)} genes that are not expressed in at least {min_cells} cells in every batch."
    )
    data = data[:, filtered_genes_all]
    logging.info(f"data.shape: {data.shape}")

    mt_gene_id = sc.queries.mitochondrial_genes(
        organism,
        chromosome="mitochondrion_genome" if organism == "dmelanogaster" else "MT",
        attrname="ensembl_gene_id" if use_gene_id else attrname,
    )

    if gene_id_column is None:
        warnings.warn(
            "`gene_id_column` was passed as None. Assuming that gene names are in var.index."
        )
        data.var["gene_id"] = np.array(data.var.index.copy())
    if isinstance(gene_id_column, str):
        if gene_id_column != "gene_id":
            data.var["gene_id"] = data.var[gene_id_column].values
            data.var.drop(gene_id_column, axis=1, inplace=True)
    data.var["mt"] = data.var["gene_id"].isin(
        mt_gene_id["ensembl_gene_id" if use_gene_id else attrname]
    )
    logging.info(f"Identified {sum(data.var.mt)} mitochondrial genes.")

    # Score cell cycle
    cell_cycle_genes_ref = pd.read_csv(cell_cycle_genes_reference)
    s_genes_ref = cell_cycle_genes_ref.loc[cell_cycle_genes_ref.phase == "S"][
        "geneID" if use_gene_id else "geneName"
    ]
    g2m_genes_ref = cell_cycle_genes_ref.loc[cell_cycle_genes_ref.phase == "G2/M"][
        "geneID" if use_gene_id else "geneName"
    ]
    cc_col = "gene_id"
    s_genes = data.var.loc[data.var.gene_id.isin(s_genes_ref), cc_col]
    g2m_genes = data.var.loc[data.var.gene_id.isin(g2m_genes_ref), cc_col]
    data.var_names = data.var[cc_col].values
    sc.tl.score_genes_cell_cycle(
        data, s_genes=s_genes, g2m_genes=g2m_genes, use_raw=False
    )
    data.obs["cell_cycle_diff"] = data.obs["S_score"] - data.obs["G2M_score"]
    logging.info(f"Scored cell cycle and calculated 'cell_cycle_diff'.")
    batch_counts, batch_dict = split_data_by_category(data, category="batch_id")
    logging.info(f"batch_counts: {batch_counts}")

    batch_pp_list = []
    for batch_, batch_data in batch_dict.items():
        batch_pp = preprocess_batch(
            batch_data,
            n_genes_by_counts_upper,
            pct_counts_mt_upper,
            filter_hvg,
            filter_hvg_n,
        )
        logging.info(f"batch_{batch_}_pp.shape = {batch_pp.shape}")
        batch_pp_list.append(batch_pp)
    tissue_pp = ad.concat(batch_pp_list, merge="same")
    logging.info(f"tissue_pp.shape = {tissue_pp.shape}")
    new_anndata = ad.AnnData(X=tissue_pp.X)
    new_anndata.obs = tissue_pp.obs[
        ["cell_id", "batch_id", "pct_counts_mt", "cell_cycle_diff"]
    ]
    new_anndata.var = tissue_pp.var[["gene_id"]]
    new_anndata.layers["counts"] = tissue_pp.layers["counts"]
    if target_column is not None:
        new_anndata.obs["y"] = tissue_pp.obs[["y"]]
    new_anndata.write(save_path)
    logging.info(f"Preprocessed data saved to {save_path}")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Set the parameters for preprocessing.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--read_path",
        type=str,
        help="Path to read the data file",
    )
    parser.add_argument(
        "--save_path",
        type=str,
        help="Path to save the preprocessed data",
    )
    parser.add_argument(
        "--organism",
        type=str,
        default="dmelanogaster",
        help="Organism (string) to query the mitochondrial gene symbols in ensembl biomart",
    )
    parser.add_argument(
        "--gene_id_column",
        type=nullable_string,
        nargs="?",
        default=None,
        help="Column of gene id in var. If not provided, will assume var.index contains gene id or name.",
    )
    parser.add_argument(
        "--batch_column",
        type=nullable_string,
        nargs="?",
        default=None,
        help="Column name of batches. \nIf not provided (either missing arg or empty string) will assume all cells are in the same batch.",
    )
    parser.add_argument(
        "--min_genes",
        type=int,
        default=200,
        help="Keep cells with >= (min_genes) genes",
    )
    parser.add_argument(
        "--min_cells",
        type=int,
        default=3,
        help="Keep genes that are found in >= (min_cells) cells",
    )
    parser.add_argument(
        "--n_genes_by_counts_upper",
        type=int,
        default=10000,
        help="Set the upper limit of the number of genes expressed in the counts matrix (n_genes_by_counts)",
    )
    parser.add_argument(
        "--pct_counts_mt_upper",
        type=float,
        default=5.0,
        help="Set the upper limit of the percentage of counts in mitochondrial genes",
    )
    parser.add_argument(
        "--cell_cycle_genes_reference",
        type=str,
        help="The path to read the reference file of cell cycle genes.",
    )
    parser.add_argument(
        "--filter_hvg",
        type=parse_bool,
        default=False,
        help="Filter highly variable genes?",
    )
    parser.add_argument(
        "--filter_hvg_n",
        type=int,
        default=5000,
        help="Number of highly variable genes to keep.",
    )
    parser.add_argument(
        "--target_column",
        type=nullable_string,
        nargs="?",
        help="Which column contains the labels. Don't pass or pass an empty string if no labels are available for this dataset (i.e., dataset is for inference.)",
        default=None,
    )

    parser.add_argument(
        "--use_gene_id",
        type=parse_bool,
        help="Whether the `gene_id_column` contains ENSEMBL gene ids. If false, external gene names are expected for drosophila and HGNC symbols for human.",
        default=False,
    )
    parser.add_argument(
        "--na_string",
        type=nullable_string,
        nargs="?",
        help="Which string indicates NAs in the targets.",
        default=None,
    )

    parser.add_argument(
        "--doublet_column",
        type=nullable_string,
        nargs="?",
        help="Which column contains doublet indicators. Pass an empty string to not perform doublet removal.",
        default=None,
    )

    parser.add_argument(
        "--doublet_value",
        type=nullable_string,
        nargs="?",
        help="Which value in the doublet column indicates doublets.",
        default=None,
    )

    parser.add_argument(
        "--use_index_barcode",
        type=parse_bool,
        help="Whether the obs index contains the cell barcode.",
        default=False,
    )

    parser.add_argument(
        "--barcode_column",
        type=nullable_string,
        nargs="?",
        help="Which column contains the cell barcode. If set, `use_index_barcode` should be False.",
        default=None,
    )

    parser.add_argument(
        "--attrname",
        type=nullable_string,
        nargs="?",
        help="Biomart attribute field to return if gene ids are not used (otherwise this does not have to be set). See https://scanpy.readthedocs.io/en/stable/generated/scanpy.queries.mitochondrial_genes.html for details. Typically `external_gene_name` for drosophila and `hgnc_symbol` for other species.",
        default=None,
    )

    args = parser.parse_args()
    main(
        read_path=args.read_path,
        save_path=args.save_path,
        organism=args.organism,
        gene_id_column=args.gene_id_column,
        batch_column=args.batch_column,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        n_genes_by_counts_upper=args.n_genes_by_counts_upper,
        pct_counts_mt_upper=args.pct_counts_mt_upper,
        cell_cycle_genes_reference=args.cell_cycle_genes_reference,
        filter_hvg=args.filter_hvg,
        filter_hvg_n=args.filter_hvg_n,
        target_column=args.target_column,
        use_gene_id=args.use_gene_id,
        na_string=args.na_string,
        doublet_column=args.doublet_column,
        doublet_value=args.doublet_value,
        use_index_barcode=args.use_index_barcode,
        barcode_column=args.barcode_column,
        attrname=args.attrname,
    )
