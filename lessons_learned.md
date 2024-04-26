# Lessons learned

## Models need to be trusted

Our first major takeaway concerns model trust and how to achieve it for cell-type annotation via supervised methods. In particular, the successful deployment of supervised ML methods that are supposed to work in conjunction with domain experts requires careful calibration of the method and the expectations of domain experts.

In our case, the domain experts are primarily curators and the supervised ML task is the automated cell-type annotation via reference scRNA-seq datasets. Throughout the project, we have found, and we now argue that it is crucial for automated cell-type annotation methods to produce predictions that are not wildly unrealistic in the curator's eyes.

For curators to work successfully together with cell-type annotation methods, curators should give some degree of trust to the classifier outputs, while keeping in mind that they are far from always correct. We briefly outline several aspects that can undermine curator trust, why they are important in our context, and how they were remedied in our project.

Inconsistent probabilities in the ontology: Inconsistent probabilities in a DAG-structured hierarchical classification problem are well-known in the literature [1] and have also been briefly explored for cell-type annotation before [2]. Briefly, inconsistent probabilities occur when for a node n, we have an arbitrary ancestor n’ for which we have P(n) > P(n’). From a human point of view, this is of course an impossible scenario since a node being annotated as n implies n’, thus P(n) must be less than or equal to P(n’) in all cases. Hence, cells in which these inconsistencies occur must be corrected using appropriate methods [1, 3], to prevent a loss of trust in the model predictions.

We explored several solutions to the inconsistency problem, including post-hoc correction methods [1], global [3], and local methods [4]. After a benchmark, we found that similar to previous work, relatively simple post-hoc methods such as isotonic regression worked best for correcting inconsistencies, keeping in mind both computation time and prediction performance [1].

Inconsistent paths through the ontology: The second problem relates to how a path is chosen through the ontology for a particular cell. The multi-path multi-label formulation is more common in the ML literature: For each sample, we choose one or more paths through the ontology, where each node is annotated with a probability.

Our problem, on the other hand, is that of single-path multi-label hierarchical classification: For each cell, we want to choose a single path through the ontology, meaning at a particular height multiple nodes can be chosen, but they must all lead to a common descendant.

While this problem has been relatively sparsely explored in the literature, it offers an interesting methodological challenge. A typical method to choose a path for this problem thresholds all probabilities using a user-defined cutoff, followed by taking the argmax of all leaf nodes in this newly induced (through the cutoff) sub-graph [2]. This approach suffers from a central problem called non-monotonicity: When the threshold is increased, there is no guarantee that the newly chosen set of nodes (induced by the chosen path) is a subset of the previously chosen set of nodes. This however, is a very desirable property in practice: While very severe violations of this sort are rare, it is theoretically possible for a model to predict a neuron subtype for a threshold of 0.5 and then predict a glial cell for a threshold of 0.8. Once again, this behavior is very illogical from a human perspective and could impede model trust.

To remedy this issue, we proposed several simple solutions that enforce monotonicity either post-hoc or during training:
Training a second leaf-node multi-class classification model. We then take the predicted path as the path induced by the predicted leaf node for each cell and simply threshold the nodes induced by the predicted leaf node for each cutoff level.
Choosing the path as the argmax over all leaf nodes. Similar to the first strategy, we take the predicted path as the path induced by the argmax leaf node for each cell and again threshold the nodes induced by the predicted leaf node for each cutoff level.
Employing a straight-through gumbel softmax (ST GS) during training. Using an ST GS, we can enforce all probabilities except those belonging to one selected path to be zero. In particular, we inputted the leaf node logits into the ST GS during training, thus jointly training the probabilities throughout the ontology and which path is most likely to be correct.
When benchmarking these three proposals, we found the first and arguably simplest method to work the best. Future work might further explore different strategies to enforce monotonicity.

Model calibration: Lastly, model calibration is paramount for user trust in a cell-type annotation model. Informally, model calibration refers to whether a model “knows what it does not know”. Said another way, for a binary classification problem, if we have 10 samples all of which have P(Y=1)=0.6, are 4 of them actually Y=0?

We benchmarked model calibration at each ontology node for our trained models using the top-label Expected Calibration Error and found that our models were mostly very well-calibrated.

That said, there were some outliers for which calibration was very bad, but we believe this to be mostly dataset-related, especially since we use relatively simple linear models, which are not known to typically suffer from calibration issues.

## Modern deep-learning and other non-linear methods do not outperform

Our second big takeaway was that non-linear methods such as deep learning and gradient boosting did not outperform linear methods such as logistic regression. This mirrors previous benchmark papers, which, roughly speaking, showed that as long as your reference scRNA-seq dataset is good, almost all methods perform very well, independent of model class [5].

Due to this, we ended up deploying a relatively simple model overall, which is both fast to train and inference and can be easily maintained since its dependencies are minimal (we use scikit-learn, but it could also easily be written by hand using a standard L-BFGS optimizer from e.g., Scipy).

In addition, interestingly, modern deep-learning-based techniques for correcting probability inconsistencies (see above) also did not outperform simple and much older techniques such as isotonic regression [1, 3]. This suggests that such simple post-hoc techniques should also be tried for the multi-path multi-label hierarchical classification problem for which deep-learning-based methods are now standard in the ML literature [3].

## Interface and integration details are crucial

Lastly, via discussions had throughout the project, we came to the conclusion that the interface and integration details are crucial for the successful deployment of an automated cell-type annotation method.

Regarding the interface, we are speaking primarily about how model predictions are presented to curators or other model end-users. We quickly found that displaying full model predictions throughout the ontology was prohibitive since this would require displaying a matrix that has n_cells rows and n_terms columns, where typically n_terms > 100. Thus, we settled on primarily interfacing the predicted path with the user, displaying only the predicted terms (typically < 10) for each cell to the user, along with their probabilities. This also allows the users to change the probability thresholds, lengthening or shortening the path (although we do not have an implemented interface for this yet).

On the integration side, our main takeaway and conclusion was to keep things as simple and modular as possible. To that effect, we implemented three self-contained Snakemake workflows, that serve the following purposes:

Determine the most likely tissue for scRNA-seq datasets for which we do not know the tissue
Train and save a model for an annotated scRNA-seq input dataset for a particular ID, organism, and tissue
Produce predicted cell types for an unannotated scRNA-seq input dataset for a particular ID, organism, and tissue

All three of our workflows have minimal dependencies and can be run fully containerized via singularity within Snakemake. Thus, we had hoped to ease the integration within pipelines and other use cases within ASAP and BGEE.

# Bibliography

[1] Obozinski, Guillaume, et al. "Consistent probabilistic outputs for protein function prediction." Genome Biology 9 (2008): 1-19.

[2] Bernstein, Matthew N., et al. "CellO: Comprehensive and hierarchical cell type classification of human cells with the Cell Ontology." Iscience 24.1 (2021).

[3] Giunchiglia, Eleonora, and Thomas Lukasiewicz. "Coherent hierarchical multi-label classification networks." Advances in neural information processing systems 33 (2020): 9662-9673.

[4] Silla, Carlos N., and Alex A. Freitas. "A survey of hierarchical classification across different application domains." Data mining and knowledge discovery 22 (2011): 31-72.

[5] Abdelaal, Tamim, et al. "A comparison of automatic cell identification methods for single-cell RNA sequencing data." Genome biology 20 (2019): 1-19.
