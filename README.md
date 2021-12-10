# GATTACA

[![License](https://img.shields.io/github/license/Feat-FeAR/GATTACA)](https://choosealicense.com/licenses/mit/)
[![Issues](https://img.shields.io/github/issues/Feat-FeAR/GATTACA)](https://github.com/Feat-FeAR/GATTACA/issues)
[![R Version](https://img.shields.io/badge/R%20version-4.1.1-informational)](https://www.r-project.org/)
[![Latest docker image size](https://img.shields.io/docker/image-size/mrhedmad/gattaca/latest)](https://hub.docker.com/repository/docker/mrhedmad/gattaca)


<b>GATTACA</b> is a short for <b>General Algorithm for The Transcriptional Analysis by one-Channel Arrays</b>.
It consists of a set of R scripts for the analysis of Transportome Transcription in Cancer Datasets.
As the name suggests, GATTACA was originally written with (one-color/high-density) microarray technology in mind, however it can be easily adapted to analyze RNA-Seq data as well.

GATTACA assumes data to be already background-subtracted, log2-transformed, and interarray-normalized.
To this purpose, use the platform-specific <i>*_to_Expression.R</i> script.

The available docker containers live [on Docker Hub](https://hub.docker.com/repository/docker/mrhedmad/gattaca#).

> **Disclaimer**: These instructions are given with a unix-like environment in mind. The program was developed and tested in an Ubuntu installation. See issue #3.

## Installation

These scripts are bundled with their dependencies in a Docker container to make usage as painless and reproducible as possible. To install, follow these instructions:

1. Install docker. A tutorial can be found [at the official docker website](https://docs.docker.com/get-docker/).
2. Get the GATTACA script. For example, you could run:
    ```bash
    curl https://raw.githubusercontent.com/Feat-FeAR/GATTACA/main/GATTACA --create-dirs -o ~/bin/GATTACA && chmod +x ~/bin/GATTACA
    ```
3. Run `GATTACA -h` to begin. The `GATTACA` script is the access point to all the other scripts in the repo. Read the usage section in this readme to learn more. The script takes care of downloading the containers and running them.

### Troubleshooting
If the `GATTACA` script cannot be run due to permissions faults, run `chmod +x GATTACA` to give it run permissions, or give it to bash directly using `/bin/bash GATTACA`.

For more help, [open an issue](https://github.com/Feat-FeAR/GATTACA/issues/new).

## Usage
Using the scripts is as simple as starting the `GATTACA` bash script. The various help messages (accessed using `-h`) detail the usage of each subcommand.

These are the currently available subcommands:
- `GATTACA init`: Used to create configuration files for `GATTACA run`.
- `GATTACA run`: Runs differential expression analysis on expression data.
- `GATTACA annotate`: Annotates expression data created by other commands with annotations from a variety of bioconductor packages.
- `GATTACA prepaffy`: Generates expression matrices from a collection of affymetrix .CEL files.

The full GATTACA manual, which covers installation, usage and result interpretation can be downloaded [here](https://github.com/Feat-FeAR/GATTACA/tree/main/docs/GATTACA_manual.pdf), from this repository, inside the `/docs` folder.

## Building the image
If you wish to rebuild the Docker container, simply clone the repo, and run `docker build .`. More information on the [docker build documentation page](https://docs.docker.com/engine/reference/commandline/build/).
Note that rebuilding the image can take a long, long time due to the need to install required packages (that need to be compiled).
