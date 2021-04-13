# **SAMBA: Standardized and Automated MetaBarcoding Analyses workflow**.

[![SAMBA version](https://img.shields.io/badge/samba%20version-v3.0.1-red?labelColor=000000)](https://www.nextflow.io/)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![Run with with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![Run with with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![SeBiMER Docker](https://img.shields.io/badge/docker%20build-SeBiMER-yellow?labelColor=000000)](https://hub.docker.com/u/sebimer)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io) and [nf-core best practices](https://nf-co.re/developers/guidelines), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility (please only use [`Conda`](https://conda.io/miniconda.html) as a last resort)

iii. Download the pipeline and test it on a minimal dataset with a single command

* for short reads test :
```bash
nextflow run main.nf -profile shortreadstest,<docker/singularity/conda>
```

* for long reads test :
```bash
nextflow run main.nf -profile longreadstest,<docker/singularity/conda>
```

> To use samba on a computing cluster, it is necessary to provide a configuration file for your system. For some institutes, this one already exists and is referenced on [nf-core/configs](https://github.com/nf-core/configs#documentation). If so, you can simply download your institute custom config file and simply use `-c <institute_config_file>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

```bash
nextflow run main.nf -profile <docker/singularity/conda>,custom [-c <institute_config_file>]
```

See [usage docs](docs/usage.md) for a complete description of all of the options available when running the pipeline.

## Documentation

The samba workflow comes with documentation about the pipeline, found in the `docs/` directory:

1. [Introduction](docs/usage.md#introduction)
2. [Pipeline installation](docs/usage.md#install-the-pipeline)
    * [Local installation](docs/usage.md#local-installation)
    * [Adding your own system config](docs/usage.md#your-own-config)
3. [Running the pipeline](docs/usage.md#running-the-pipeline)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

Here is an overview of the many steps available in samba pipeline:

![SAMBA Workflow](./docs/images/samba-v2.0.0.png)

At the end of samba pipeline execution, you get an interactive HTML report that's look like this one:

![SAMBA report](docs/images/samba-report.gif)

Full report description is available in [samba pipeline documentation](docs/output.md).

## Credits

samba is written by [SeBiMER](https://ifremer-bioinformatics.github.io/), the Bioinformatics Core Facility of [IFREMER](https://wwz.ifremer.fr/en/).

## Contributions

We welcome contributions to the pipeline. If such case you can do one of the following:
* Use issues to submit your questions 
* Fork the project, do your developments and submit a pull request
* Contact us (see email below) 

## Support

For further information or help, don't hesitate to get in touch with the samba developpers: 

![samba email](assets/samba-email-address-image.png)

## Citation

<!-- If you use  samba for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

### References 

References databases (SILVA v132, PR2, UNITE) are available on IFREMER FTP at [ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/2019.10](ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/SAMBA/2019.10).

Training dataset used from [Qiime2 Tutorial] (https://docs.qiime2.org/2019.7/tutorials/atacama-soils), [associated publication](https://msystems.asm.org/content/2/3/e00195-16).
