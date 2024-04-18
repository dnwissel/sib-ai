# Cell-type classification tools for the SIB AI project

## Execution
The only dependencies are snakemake (>= 7.30.1) and singularity. Afterwards, you can run both the tissue annotator and the cell type classification using:

```sh
snakemake --use-conda --use-singularity
```

All configuration happens in each respective workflow's config.yaml file. Most of the parameters are hopefully self-documenting, although the command line args in each Python script should usually also have somewhat useful help messages.

Just let me know if you would like more documentation in any part of this repo.

## Limitations
- The tissue has only been tested in a very limited manner on FCA data. It is also rather slow so if this is a common task, we should rediscuss the model that we're using and/or how we're applying it.
