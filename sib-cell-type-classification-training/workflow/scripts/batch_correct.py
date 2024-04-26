import argparse
from typing import List

import scanpy as sc
import scvi
from anndata import AnnData
from pandas.api.types import is_numeric_dtype
from typeguard import typechecked


@typechecked
def train_scanvi(
    train: AnnData,
    batch_key: str,
    covariate_keys: List[str],
    model_path: str,
    n_latent: int,
    n_layers: int,
    n_epochs_scvi: int,
    n_epochs_scanvi: int,
    label_column: str,
    embedding_path: str,
) -> None:
    categorical_keys = []
    continuous_keys = []
    for key in covariate_keys:
        if not is_numeric_dtype(train.obs[key]):
            categorical_keys.append(key)
        else:
            continuous_keys.append(key)
    scvi.model.SCVI.setup_anndata(
        train,
        layer="counts",
        batch_key=batch_key,
        categorical_covariate_keys=categorical_keys,
        continuous_covariate_keys=continuous_keys,
    )
    vae_ref = scvi.model.SCVI(train, n_latent=n_latent, n_layers=n_layers)
    vae_ref.train(n_epochs_scvi)
    train.obs["labels_scanvi"] = train.obs[label_column].values
    vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
        vae_ref,
        adata=train,
        labels_key="labels_scanvi",
        unlabeled_category="Unknown",
    )
    vae_ref_scan.view_anndata_setup(train)
    vae_ref_scan.train(max_epochs=n_epochs_scanvi)
    vae_ref_scan.save(model_path, overwrite=True)
    train.obsm["X_scANVI"] = vae_ref_scan.get_latent_representation(train)
    train.write(embedding_path)
    return None


@typechecked
def main(
    data_path: str,
    batch_key: str,
    covariate_keys,
    model_path: str,
    n_latent: int,
    n_layers: int,
    n_epochs_scvi: int,
    n_epochs_scanvi: int,
    label_column: str,
    seed: int,
    embedding_path: str,
) -> int:
    scvi.settings.seed = seed
    data = sc.read_h5ad(data_path)
    train_scanvi(
        train=data,
        batch_key=batch_key,
        covariate_keys=covariate_keys,
        model_path=model_path,
        n_latent=n_latent,
        n_layers=n_layers,
        n_epochs_scvi=n_epochs_scvi,
        n_epochs_scanvi=n_epochs_scanvi,
        label_column=label_column,
        embedding_path=embedding_path,
    )
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Set the parameters for batch correction.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        type=str,
        help="Path to read the training data for scANVI from.",
    )

    parser.add_argument(
        "--model_path",
        type=str,
        help="Path to dump the trained scANVI model to.",
    )

    parser.add_argument(
        "--batch_key",
        type=str,
        help="Column name of the batch id.",
    )
    parser.add_argument(
        "--label_key",
        type=str,
        help="Column name of the label (i.e., cell-type).",
    )

    parser.add_argument(
        "--embedding_path",
        type=str,
        help="Path to write the embedding of the training data to.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        help="Initial state of the pseudo-random number generator.",
    )
    parser.add_argument(
        "--covariate_keys",
        type=str,
        nargs="*",
        help="Covariate keys to be corrected for using scANVI.",
    )
    parser.add_argument(
        "--n_latent",
        type=int,
        help="Number of latent dimensions of the scANVI model.",
        default=30,
    )

    parser.add_argument(
        "--n_layers",
        type=int,
        help="Number of hidden layers of the scANVI model.",
        default=2,
    )

    parser.add_argument(
        "--n_epochs_scvi",
        type=int,
        help="Number of epochs to train the scVI model (used to warm-start scANVI).",
        default=100,
    )

    parser.add_argument(
        "--n_epochs_scanvi",
        type=int,
        default=25,
        help="Number of epochs to train the scANVI model.",
    )

    args = parser.parse_args()
    main(
        data_path=args.data_path,
        batch_key=args.batch_key,
        covariate_keys=args.covariate_keys,
        model_path=args.model_path,
        n_latent=args.n_latent,
        n_layers=args.n_layers,
        n_epochs_scvi=args.n_epochs_scvi,
        n_epochs_scanvi=args.n_epochs_scanvi,
        label_column=args.label_key,
        seed=args.seed,
        embedding_path=args.embedding_path,
    )
