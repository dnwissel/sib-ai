import argparse

import joblib
import networkx as nx
import numpy as np
import pandas as pd
import scanpy
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MultiLabelBinarizer, StandardScaler
from typeguard import typechecked


@typechecked
def load_full_hier(path: str):
    hier = pd.read_csv(path, sep="\t")
    edges = hier.to_records(index=False).tolist()
    G = nx.DiGraph(edges).reverse()
    roots_label = [v for v, d in G.in_degree() if d == 0]
    if len(roots_label) > 1:
        for node in roots_label:
            G.add_edge("root", node)
        roots_label += ["root"]
    G.add_edge("root", "FBbt:00005073")
    return G, roots_label


class MultiLabelLogReg(LogisticRegression):
    def fit(self, X, y, sample_weight=None):
        if np.all(y == 1):
            self.coef_ = np.zeros((2, X.shape[1]))
            self.intercept_ = np.array([-100000, 100000])
            self.classes_ = np.unique(y)
        elif np.all(y == 0):
            self.coef_ = np.zeros((2, X.shape[1]))
            self.intercept_ = np.array([100000, -100000])
            self.classes_ = np.unique(y)
        else:
            super().fit(X, y, sample_weight)
        return None


@typechecked
def main(
    data_path: str,
    model_path_leaf_nodes: str,
    model_path_non_leaf_nodes: str,
    model_path_binarizer: str,
    hierarchy_path: str,
    seed: int,
    label_key: str,
    n_jobs: int = 1,
    C: float = 100.0,
    tol: float = 0.01,
    class_weight: str = "balanced",
    max_iter: int = 5000,
) -> int:
    mlb = MultiLabelBinarizer()
    data = scanpy.read_h5ad(data_path)
    features = np.array(data.obsm["X_scANVI"])
    y = data.obs[label_key].values
    hierarchy, _ = load_full_hier(hierarchy_path)
    y_multi_label = [
        tuple(set(nx.ancestors(hierarchy, node)).union(set([node]))) for node in y
    ]

    y_multi_label = mlb.fit_transform(y_multi_label)

    leaf_node_pipeline = make_pipeline(
        VarianceThreshold(),
        StandardScaler(),
        LogisticRegression(
            class_weight=class_weight,
            random_state=seed,
            # Irrelevant when `multi_class`==`multinomial`
            n_jobs=1,
            max_iter=max_iter,
            tol=tol,
            C=C,
        ),
    )
    non_leaf_node_pipeline = make_pipeline(
        VarianceThreshold(),
        StandardScaler(),
        MultiOutputClassifier(
            estimator=MultiLabelLogReg(
                class_weight=class_weight,
                random_state=seed,
                max_iter=max_iter,
                tol=tol,
                C=C,
                # Nothing to parallelize on the inside
                n_jobs=1,
            ),
            n_jobs=n_jobs,
        ),
    )

    leaf_node_pipeline.fit(features, y)
    non_leaf_node_pipeline.fit(features, y_multi_label)
    _ = joblib.dump(leaf_node_pipeline, model_path_leaf_nodes)
    _ = joblib.dump(non_leaf_node_pipeline, model_path_non_leaf_nodes)
    _ = joblib.dump(mlb, model_path_binarizer)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Set the parameters for training the classification model.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        type=str,
        help="Path to read the training data for the classifier model from.",
    )

    parser.add_argument(
        "--model_path_leaf_nodes",
        type=str,
        help="Path to dump the trained classifier model to.",
    )

    parser.add_argument(
        "--model_path_non_leaf_nodes",
        type=str,
        help="Path to dump the trained classifier model to.",
    )

    parser.add_argument(
        "--model_path_binarizer",
        type=str,
        help="Path to dump the trained binarizer to.",
    )

    parser.add_argument(
        "--hierarchy_path",
        type=str,
        help="Path to read the hierarchy TSV from.",
    )
    parser.add_argument(
        "--label_key",
        type=str,
        help="Column name of the label (i.e., cell-type).",
    )

    parser.add_argument(
        "--seed",
        type=int,
        help="Initial state of the pseudo-random number generator.",
    )

    parser.add_argument(
        "--n_jobs",
        type=int,
        help="How many threads to parallelize over when fitting the per-node classifiers.",
    )

    parser.add_argument(
        "--C",
        type=float,
        help="Inverse regularization parameter for `sklearn.linear_model.LogisticRegression`.",
    )

    parser.add_argument(
        "--tol",
        type=float,
        help="Convergence tolerance parameter for `sklearn.linear_model.LogisticRegression`.",
    )

    parser.add_argument(
        "--class_weight",
        type=str,
        help="Class weight parameter for `sklearn.linear_model.LogisticRegression`.",
    )

    parser.add_argument(
        "--max_iter",
        type=float,
        help="Max iteration parameter for `sklearn.linear_model.LogisticRegression`.",
    )

    args = parser.parse_args()
    main(
        data_path=args.data_path,
        model_path_leaf_nodes=args.model_path_leaf_nodes,
        model_path_non_leaf_nodes=args.model_path_non_leaf_nodes,
        model_path_binarizer=args.model_path_binarizer,
        hierarchy_path=args.hierarchy_path,
        label_key=args.label_key,
        seed=args.seed,
        n_jobs=args.n_jobs,
        C=args.C,
        tol=args.tol,
        class_weight=args.class_weight,
        max_iter=args.max_iter,
    )
