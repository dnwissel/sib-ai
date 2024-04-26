# SIB hierarchical cell-type classification (training)

## General idea

This snakemake workflow takes an input scRNA-seq training dataset with cell-type annotations. After preprocessing and batch correction on the input dataset, models to predict probabilities for the CL are trained such that they can easily be deployed on other datasets later.

## How to run

- Ensure that the `config.yaml` file reflects your desired configuration (see below for details). Recall that any config parameter can also be overwritten on the command line as is usual for snakemake.
- Start the workflow using `snakemake --use-conda --use-singularity --cores 1 --singularity-args "-B $HOME/projects/cr_sib:$HOME/projects/cr_sib"`
- The following parameters do not have defaults in the config file and must thus always be bound on the command line:
  - `tissue`
  - `organism`
  - `id`
  - `output_model_path_leaf`
  - `output_model_path_non_leaf`
  - `output_model_path_binarizer`
  - `input_h5ad`
  - `output_preprocessed_h5ad`
  - `output_batch_model_dir`
  - `output_embedding_path`
  - `cell_cycle_genes_reference`
  - `hierarchy_path`
  - `attrname`

Note: Since this workflow is designed PER dataset, it does not natively (i.e., using Snakemake) parallelize and you have to parallelize manually across multiple datasets by spawning multiple Snakemake runs.

## Requirements

- `singularity-ce`
- `snakemake>=7.30.1`

## General assumptions

- The `.X` of the h5ad input file contains counts and not any kind of transformed data.

## Config file

- `gene_id_column`: In which column of `.var` of the h5ad the gene names or gene IDs are located. If an empty string is passed, assume that gene IDs/gene names are located in `.var`.index.
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
- `n_latent`: Number of latent dimensions of the embedding.
- `n_layers`: Number of hidden layers for the embedding model.
- `n_epochs_scvi`: Number of epochs to train the scVI embedding model.
- `n_epochs_scanvi`: Number of epochs to train the scANVI embedding model.
- `seed`: Random seed for reproducibility.
- `hierarchy_path`: Path to the hierarchy as a TSV.
- `target_column`: Column of `.obs` that contains the cell-type labels in the training dataset.
- `use_gene_id`: whether gene IDs are used to identify genes. If false, assume that gene names are used.
- `na_string`: Which string indicates missing annotations in target_column. Can be set to any value if no missing annotations are present.
- `doublet_column`: Which column contains doublet information.
- `doublet_value`: Which value in doublet_column indicates a doublet.
- `use_index_barcode`: Whether the index of `.obs` contains cell barcodes.
- `barcode_column`: Which column of `.obs` contains cell IDs. Only applicable if use_index_barcode is false.
- `tissue`: Which tissue the model is being deployed on.
- `organism`: Which organism the model is being deployed on.
- `id`: For which ID the model is being deployed.
- `output_model_path_leaf`: Ouput path where the .pkl file for the leaf model should be written to.
- `output_model_path_non_leaf`: Output path where the .pkl file for the non-leaf model should be written to.
- `output_model_path_binarizer`: Output path where the .pkl file for the binarizer is should be written to.
- `input_h5ad`: Input path to the training dataset as h5ad.
- `output_preprocessed_h5ad`: Output path to write the intermediate preprocessed h5ad of the training dataset.
- `output_batch_model_dir`: Output directory in which to place the trained torch model for creating embeddings.
- `output_embedding_path`: Output path to which to write the embedding of the input dataset used for training.
- `attrname`: Biomart attribute field to return if gene ids are not used (otherwise this does not have to be set). See https://scanpy.readthedocs.io/en/stable/generated/scanpy.queries.mitochondrial_genes.html for details. Typically `external_gene_name` for drosophila and `hgnc_symbol` for other species.
- `n_jobs`: How many threads to parallelize over when fitting the per-node classifiers.
- `C`: Inverse regularization parameter for `sklearn.linear_model.LogisticRegression`.
- `tol`: Convergence tolerance parameter for `sklearn.linear_model.LogisticRegression`.
- `class_weight`: Class weight parameter for `sklearn.linear_model.LogisticRegression`.
- `max_iter`: Max iteration parameter for `sklearn.linear_model.LogisticRegression`.

## Per rule results

There are three rules:

- `preprocess_train`: Preprocesses the input data and writes a preprocessed (and stripped of anything unnecessary for training) h5ad for training.
- `train_batch_correct`: Trains an embedding model for dimensionality reduction and covariate removal (including batch). Outputs an h5ad that contains the embedded data for training and a torch model that can be used to update and later embed the deployment data.
- `train_model`: Trains a multi-label model for the whole ontology and a multi-class model for the leaf nodes to choose the path through the ontology. Outputs both models and a binarizer to recover the cell type labels.
