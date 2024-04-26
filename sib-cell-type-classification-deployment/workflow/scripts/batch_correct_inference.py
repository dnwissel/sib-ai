import argparse

import pandas as pd
import scanpy
import scvi
from typeguard import typechecked


@typechecked
def main(
    data_path: str,
    embedding_path: str,
    model_path: str,
    seed: int,
    accelerator: str,
    check_val_every_n_epoch: int,
    max_epochs: int,
) -> int:
    scvi.settings.seed = seed
    data = scanpy.read_h5ad(data_path)
    scvi.model.SCANVI.prepare_query_anndata(data, model_path)
    vae_query = scvi.model.SCANVI.load_query_data(
        data,
        model_path,
    )
    vae_query.train(
        max_epochs=max_epochs,
        plan_kwargs={"weight_decay": 0.0},
        check_val_every_n_epoch=check_val_every_n_epoch,
        accelerator=accelerator,
    )
    X_test = pd.DataFrame(vae_query.get_latent_representation(data))
    X_test["cell_id"] = data.obs["cell_id"].values
    X_test.to_csv(embedding_path, index=False)

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Set the parameters for inference batch correction.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        type=str,
        help="Path to read the test data for scANVI from.",
    )

    parser.add_argument(
        "--embedding_path",
        type=str,
        help="Path to dump the test embedding to.",
    )

    parser.add_argument(
        "--model_path",
        type=str,
        help="Path where the trained scANVI model is saved.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for reproducibility.",
    )

    parser.add_argument(
        "--accelerator",
        type=str,
        help="Which accelerator to use. Must be cpu or gpu.",
        default="cpu",
    )

    parser.add_argument(
        "--check_val_every_n_epoch",
        type=int,
        help="How often to check validation loss.",
        default=10,
    )

    parser.add_argument(
        "--max_epochs",
        type=int,
        help="How many epochs to train at a maximum.",
        default=100,
    )

    args = parser.parse_args()
    main(
        data_path=args.data_path,
        embedding_path=args.embedding_path,
        model_path=args.model_path,
        seed=args.seed,
        accelerator=args.accelerator,
        check_val_every_n_epoch=args.check_val_every_n_epoch,
        max_epochs=args.max_epochs,
    )
