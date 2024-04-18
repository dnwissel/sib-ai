import argparse

import numpy as np
import pandas as pd
import scanpy
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import make_pipeline
from typeguard import typechecked


@typechecked
def main(
    read_path: str,
    master_path: str,
    master_tissue_column: str,
    master_gene_column: str,
    deploy_decide_gene_column: str,
    file_path_pre: str,
    file_path_post: str,
    pca_n_components: int,
) -> 0:
    pipe = make_pipeline(
        # PCA(n_components=pca_n_components, svd_solver="arpack"),
        RandomForestClassifier(class_weight="balanced"),
    )
    master = scanpy.read_h5ad(master_path)
    data = scanpy.read_h5ad(read_path)
    if master_gene_column == "index":
        master.var["gene_name"] = master.var.index
        master.var.reset_index(drop=True, inplace=True)
    elif master_gene_column != "gene_name":
        master.var["gene_name"] = master.var[master_gene_column].values
        master.drop(labels=[master_gene_column], axis="columns", inplace=True)

    if deploy_decide_gene_column == "index":
        data.var["gene_name"] = data.var.index
        data.var.reset_index(drop=True, inplace=True)
    elif deploy_decide_gene_column != "gene_name":
        data.var["gene_name"] = data.var[deploy_decide_gene_column].values
        data.drop(labels=[deploy_decide_gene_column], axis="columns", inplace=True)

    joint_genes = np.intersect1d(data.var.gene_name.values, master.var.gene_name.values)
    if joint_genes.shape[0] < 1000:
        raise ValueError(
            f"Expected significant overlap (1000+) between master and data gene names. Only {joint_genes.shape[0]} overlaps found."
        )
    data_gene_ix = np.array(
        [np.where(gene == data.var.gene_name)[0] for gene in joint_genes]
    ).squeeze()
    master_gene_ix = np.array(
        [np.where(gene == master.var.gene_name)[0] for gene in joint_genes]
    ).squeeze()
    print(data_gene_ix.shape)
    print(master_gene_ix.shape)
    data = data[:, data_gene_ix]
    master = master[:, master_gene_ix]

    pipe.fit(master.X, master.obs[master_tissue_column].values)
    selected_tissue = pipe.predict(data.X)
    selected_tissue = np.unique(selected_tissue, return_counts=True)[0][
        np.argmax(np.unique(selected_tissue, return_counts=True)[1])
    ]
    file_path = file_path_pre + str(selected_tissue[0]) + file_path_post + ".h5ad"
    data.write(file_path)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Set the parameters for inference batch correction.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--read_path",
        type=str,
        help="Path to read the deployment file for which the tissue is to be determined from. Must be h5ad.",
    )

    parser.add_argument(
        "--master_path",
        type=str,
        help="Path to read the master file which is used to train the tissue model from. Must be h5ad.",
    )

    parser.add_argument(
        "--master_tissue_column",
        type=str,
        help="In which column of .obs of the master file the tissue information is stored.",
    )

    parser.add_argument(
        "--master_gene_column",
        type=str,
        help="In which column of .var of the master file the gene names are stored.",
    )

    parser.add_argument(
        "--deploy_decide_gene_column",
        type=str,
        help="In which column of .var of the query file the gene names are stored.",
    )

    parser.add_argument(
        "--file_path_pre",
        type=str,
        help="Prefix of the output path.",
    )

    parser.add_argument(
        "--file_path_post",
        type=str,
        help="Postfix of the output path.",
    )

    parser.add_argument(
        "--pca_n_components",
        type=int,
        help="Number of components to run PCA with.",
    )

    args = parser.parse_args()
    main(
        read_path=args.read_path,
        master_path=args.master_path,
        master_tissue_column=args.master_tissue_column,
        master_gene_column=args.master_gene_column,
        deploy_decide_gene_column=args.deploy_decide_gene_column,
        file_path_pre=args.file_path_pre,
        file_path_post=args.file_path_post,
        pca_n_components=args.pca_n_components,
    )
