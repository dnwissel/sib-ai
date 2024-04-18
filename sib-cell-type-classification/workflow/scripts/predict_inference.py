import argparse
import json
from operator import itemgetter

import joblib
import networkx as nx
import numpy as np
import obonet
import pandas as pd
from quadprog import solve_qp
from sklearn.linear_model import LogisticRegression
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


# @typechecked
# def load_full_hier(path: str):
#     G = obonet.read_obo(path)
#     return G


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
    output_path_json: str,
    output_path_tsv: str,
    thresh: float,
) -> int:

    leaf_node_pipeline = joblib.load(model_path_leaf_nodes)
    non_leaf_node_pipeline = joblib.load(model_path_non_leaf_nodes)
    mlb = joblib.load(model_path_binarizer)

    data = pd.read_csv(data_path)
    cell_id = data["cell_id"].values
    data.drop(labels=["cell_id"], axis="columns", inplace=True)
    y = np.unique(mlb.classes_)
    features = data.to_numpy()

    predicted_total = pd.DataFrame(
        np.stack(
            [i[:, 1] for i in non_leaf_node_pipeline.predict_proba(features)], axis=1
        ),
        columns=mlb.classes_,
    )
    predicted_leaf = pd.DataFrame(
        leaf_node_pipeline.predict_proba(features),
        columns=leaf_node_pipeline[2].classes_,
    )
    hierarchy, _ = load_full_hier(hierarchy_path)
    # Partially adapted from CellO: https://github.com/deweylab/CellO/blob/master/cello/models/isotonic_regression.py
    constraint_size = np.sum(
        [len(tuple(nx.ancestors(hierarchy, node))) for node in np.unique(y)]
    )
    constraint_matrix = np.zeros((y.shape[0], constraint_size))
    b = np.zeros(constraint_size)

    Q = np.eye(predicted_total.shape[1])
    total_constraints = 0
    for node in predicted_total.columns:
        node_ix = np.where(node == predicted_total.columns)[0]
        for ancestor in tuple(nx.ancestors(hierarchy, node)):
            ancestor_ix = np.where(ancestor == predicted_total.columns)[0]

            constraint_matrix[ancestor_ix, total_constraints] = 1.0
            constraint_matrix[node_ix, total_constraints] = -1.0
            total_constraints += 1

    for row in range(predicted_total.shape[0]):
        sol = solve_qp(Q, predicted_total.iloc[row, :].values, constraint_matrix, b)
        predicted_total.iloc[row, :] = sol[0]
    predicted_leaf.index = cell_id
    predicted_total.to_csv(output_path_tsv, index=True, sep="\t")
    json_list = [None] * predicted_total.shape[0]
    leaf_ix = np.array(
        [
            np.where(leaf_name == predicted_total.columns)[0]
            for leaf_name in leaf_node_pipeline.classes_
        ]
    )
    ordered_columns = np.squeeze(predicted_total.columns.to_numpy()[leaf_ix])
    for row in range(predicted_total.shape[0]):
        chosen_leaf = ordered_columns[predicted_leaf.iloc[row,].argmax()]
        chosen_path = list(
            set(nx.ancestors(hierarchy, chosen_leaf)).union({chosen_leaf})
        )
        # chosen_path = [*nx.topological_sort(hierarchy.subgraph(chosen_path))]
        chosen_row = predicted_total[chosen_path].iloc[row, :]
        chosen_row = chosen_row[chosen_row > thresh]
        json_list[row] = sorted(
            [
                [
                    chosen_row.index[i],
                    chosen_row[i],
                    len(
                        nx.shortest_path(
                            hierarchy, source="root", target=chosen_row.index[i]
                        )
                    )
                    - 1,
                ]
                for i in range(chosen_row.shape[0])
            ],
            key=itemgetter(2),
        )

    with open(output_path_json, "w") as file:
        json.dump(
            obj={cell_id[i]: json_list[i] for i in range(len(cell_id))},
            fp=file,
            indent=4,
        )

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Set the parameters for inference.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        type=str,
        help="Path to read the test data for the classifier model from.",
    )

    parser.add_argument(
        "--model_path_leaf_nodes",
        type=str,
        help="Path to read the trained classifier model from.",
    )

    parser.add_argument(
        "--model_path_non_leaf_nodes",
        type=str,
        help="Path to read the trained classifier model from.",
    )

    parser.add_argument(
        "--model_path_binarizer",
        type=str,
        help="Path to read the binarizer from.",
    )

    parser.add_argument(
        "--hierarchy_path",
        type=str,
        help="Path to read the hierarchy TSV from.",
    )
    parser.add_argument(
        "--output_path_json",
        type=str,
        help="Path to write the JSON path file to.",
    )

    parser.add_argument(
        "--output_path_tsv",
        type=str,
        help="Path to write the TSV file to.",
    )

    parser.add_argument(
        "--thresh",
        type=float,
        default=0.5,
        help="Minimum probability to be contained in the output path.",
    )

    args = parser.parse_args()
    main(
        data_path=args.data_path,
        model_path_leaf_nodes=args.model_path_leaf_nodes,
        model_path_non_leaf_nodes=args.model_path_non_leaf_nodes,
        model_path_binarizer=args.model_path_binarizer,
        hierarchy_path=args.hierarchy_path,
        output_path_json=args.output_path_json,
        output_path_tsv=args.output_path_tsv,
        thresh=args.thresh,
    )
