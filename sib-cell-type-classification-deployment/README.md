# SIB hierarchical cell-type classification (deployment)

## General idea

This snakemake workflow takes an input scRNA-seq training dataset without cell-type annotations and a corresponding trained cell-type classification model. After preprocessing and batch correction on the input dataset, probabilities for the CL are predicted on the deployment datasets, along with one chosen path through the CL.

## How to run

- Ensure that the `config.yaml` file reflects your desired configuration (see below for details). Recall that any config parameter can also be overwritten on the command line as is usual for snakemake.
- Start the workflow using `snakemake --use-conda --use-singularity --cores 1 --singularity-args "-B $HOME/projects/cr_sib:$HOME/projects/cr_sib"`
- The following parameters do not have defaults in the config file and must thus always be bound on the command line:
  - `output_path_json`
  - `output_path_tsv`
  - `input_h5ad`
  - `output_preprocessed_h5ad`
  - `id`
  - `organism`
  - `tissue`
  - `attrname`
  - `input_batch_model_dir`
  - `output_embedding_path`
  - `input_model_path_leaf`
  - `input_model_path_non_leaf`
  - `input_model_path_binarizer`
  - `cell_cycle_genes_reference`
  - `hierarchy_path`

Note: Since this workflow is designed PER dataset, it does not natively (i.e., using Snakemake) parallelize and you have to parallelize manually across multiple datasets by spawning multiple Snakemake runs.

## Requirements

- `singularity-ce`
- `snakmake>=7.30.1`

## General assumptions:

- The `.X` of the h5ad input file contains counts and not any kind of transformed data.

## Config file

- `gene_id_column`: In which column of var of the h5ad the gene names or gene IDs are located. If an empty string is passed, assume that gene IDs/gene names are located in var.index.
- `batch_column`: In which column of abs the batch ID is. If there are no batches either pass all the same batch ID or pass an empty string.
- `min_genes`: min number of counts per gene per batch to keep a gene (set to zero since I assume you may have filtered before).
- `min_cells`: min number of cells where a gene is expressed per batch to keep a gene (set to zero since I assume you may filter before).
- `n_genes_by_counts_upper`: Upper bound on the number of genes with at least 1 count in a cell. Used for cell filtering. Should be set to a very large value to disable.
- `pct_counts_mt_upper`: Upper bound on percentage mitochondrial counts per cell. Used for cell filtering. Should be set to a very large value to disable.
- `cell_cycle_genes_reference`: Path for cell cycle reference genes. Need to add additional ones with the same naming scheme in the config for other organisms. Please note that both gene IDs and gene names are required in the file (see example on GDrive) to enable working with either gene IDs or gene names.
- `filter_hvg`: Flag whether HVG selection should be performed before embedding.
- `filter_hvg_n`: Number of HVGs to select.
- `mambaforge_container_version`: Mambaforge version to use within singularity containers. You probably never need to touch this but can bump to the latest version if you don't care. I didn't use the latest for strict reproducibility.
- `covariate_keys`: Covariates to correct for during embedding, in addition to batch id.
- `seed`: Random seed for reproducibility.
- `hierarchy_path`: Path to the hierarchy as a TSV.
- `check_val_every_n_epoch_deploy_batch`: How often (every X epochs) the validation loss is checked for early stopping during the integration of any new batch ids into the embedding model.
- `max_epochs_deploy_batch`: Maximum total number of epochs to train during the integration of any new batch IDs into the embedding model.
- `thresh`: Threshold to shorten the path JSON. Only nodes with predicted probability larger than `thresh` are retrained in the output JSON per cell.
- `use_gene_id`: whether gene IDs are used to identify genes. If false, assume that gene names are used.
- `use_index_barcode_deploy`: Same as use_index_barcode but for deployment datasets (since it may differ).
- `barcode_column_deploy`: Same as barcode_column but for deployment datasets (since it may differ).
- `doublet_column_deploy`: Same as doublet_column but for deployment datasets (since it may differ).
- `doublet_value_deploy`: Same as doublet_value but for deployment datasets (since it may differ).
- `gene_id_column_deploy`: Same as gene_id but for deployment datasets (since it may differ).
- `batch_column_deploy`: Same as batch_column but for deployment datasets (since it may differ).
- `tissue`: Which tissue the model is being deployed on.
- `organism`: Which organism the model is being deployed on.
- `id`: For which ID the model is being deployed.
- `output_path_json`: Path to which to write the predicted output path as a JSON.
- `output_path_tsv`: Path to which to write the probabilities for the full ontology as a TSV.
- `input_h5ad`: Path from which to read the input h5ad for which cell types are to be determined.
- `output_preprocessed_h5ad`: Path to which to write the preprocessed h5ad.
- `attrname`: Biomart attribute field to return if gene ids are not used (otherwise this does not have to be set). See https://scanpy.readthedocs.io/en/stable/generated/scanpy.queries.mitochondrial_genes.html for details. Typically `external_gene_name` for drosophila and `hgnc_symbol` for other species.
- `input_batch_model_dir`: Input directory where the batch correction model is located.
- `output_embedding_path`: Output path to which the embedded input dataset should be written.
- `input_model_path_leaf`: Input path where the `.pkl` file for the leaf model is located.
- `input_model_path_non_leaf`: Input path where the `.pkl` file for the non-leaf model is located.
- `input_model_path_binarizer`: Input path where the `.pkl` file for the binarizer is located.

## Per rule results

There are three rules:

- `preprocess_deploy`: This preprocesses deploys datasets and writes a preprocessed (and stripped of anything unnecessary for training) h5ad for deployment. Differs from - preprocess_train in only very few ways, mostly related to the handling of annotation columns (since deploy typically has none).
- `deploy_batch_correct`: Adds any new batches to the scANVI model (via architecture surgery) and then embeds the deployment data and outputs a CSV of the embedded deployment data.
- `deploy_model`: Deploys the multi-label model together with the multi-class model to choose a path through the hierarchy. Violations (cases where P(ancestor) < P(child)) are corrected using isotonic regression. Outputs a TSV containing probabilities for each node in the ontology and a JSON that encodes a path through the ontology for each cell, thresholded (see `thresh` in the config file).
