# SIB tissue classifier

## General idea

This snakemake workflow takes a scRNA-seq dataset without cell-type annotations and determines which tissue type the cells likely originate from. This is done via supervised ML based on the tissue types contained in a reference atlas such as the Fly Cell Atlas.

## How to run

- Ensure that the `config.yaml` file reflects your desired configuration (see below for details). Recall that any config parameter can also be overwritten on the command line as is usual for snakemake.
- Start the workflow using `snakemake --use-conda --use-singularity --cores 1 --singularity-args "-B $HOME/projects/cr_sib:$HOME/projects/cr_sib"`
- The following parameters do not have defaults in the config file and must thus always be bound on the command line:
  - `output_h5ad`
  - `input_h5ad`
  - `input_master`
  - `seed`
  - `organism`
  - `id`

Note: Since this workflow is designed PER dataset, it does not natively (i.e., using Snakemake) parallelize and you have to parallelize manually across multiple datasets by spawning multiple Snakemake runs.

## Requirements

- `singularity-ce`
- `snakmake>=7.30.1`

## General assumptions

- The `.X` of the h5ad input files (both the master and the file to be predicted) contain an appropriate quantity (e.g., counts, log-counts, or a transformation thereof) for determining the tissue type.
- The `.X` between the master and the file to be predicted are the same (e.g., both counts).
- The master_gene_column and deploy_decide_gene_column contain features of the same type, i.e., both contain gene IDs or both contain gene names.

## Limitations

- The tissue classifier was added quite late in the process of this pipeline and has not been extensively tested. Thus, results should be treated with extra care before using them.
- The current setup retrains the master model for each dataset on which it is deployed to ensure maximum overlap of genes. This increases computation time somewhat.

# Config file

- `mambaforge_container_version`: Mambaforge version to use within singularity containers. You probably never need to touch this but can bump to the latest version if you don't care. I didn't use the latest for strict reproducibility.
- `master_tissue_column`: The column in the `.obs` of the master h5ad that contains the labels, that is the tissue labels.
- `master_gene_column`: The column in the `.var` of the master h5ad that contains the feature names, that is gene names or gene IDs.
- `deploy_decide_gene_column`: The column in the `.var` of the query h5ad that contains the feature names, that is gene names or gene IDs.
- `pca_n_components`: How many components are to be kept after a truncated SVD on the original feature matrix of the master and query file.
- `seed`: Random seed for the pseudo-random number generator.
- `output_h5ad`: Output file name for the output h5ad with added tissue information.
- `input_h5ad`: Input file name for the dataset for which the tissue is to be determined.
- `input_master`: Input file name for the dataset which is used to train the model to perform the tissue typing.
- `organism`: Which organism we are working on. Only used to make the log path informative, can be arbitrary.
- `id`: ID for the sample. Only used to make the log path informative, can be arbitrary.

## Per rule results

There is one rule:

- `decide_unknown_tissue`: This determines the likely tissue type of an input h5ad file and outputs the input h5ad file unmodified but with a new `tissue` column in the `.obs` data frame that contains the predicted tissue type. Note that we enforce a uniform tissue type across all cells of a dataset.
