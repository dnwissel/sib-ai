# Cell-type classification tools for the SIB AI project

## Example use case
Please find a tarball containing full example data on [Zenodo](TODO-Sunday). The example use case trains a head tissue model based on FCA data and deploy it one dataset known to contain cells coming from head tissue and one dataset predicted to come from head tissue.

The example use case only has four dependencies:

```
scanpy
singularity
typeguard
snakemake
```

You may reproduce the example use case by executing the following command in the root directory of the tarball.

```
bash deploy.sh
```

## Documentation
For documentation, please see the README files within each of the lorem:

- [Training](https://github.com/dnwissel/sib-ai/sib-cell-type-classification-training/README.md)
- [Deployment](https://github.com/dnwissel/sib-ai/sib-cell-type-classification-deployment/README.md)
- [Tissue classification](https://github.com/dnwissel/sib-ai/sib-tissue-classifier/README.md)

A short writeup about lessons learned throughout the project may be accessed [here](https://github.com/dnwissel/sib-ai/lessons_learned.md).

## Questions and issues
In case you have any questions, please get in touch via email (dwissel@ethz.ch) or feel free to open an issue. Thanks!