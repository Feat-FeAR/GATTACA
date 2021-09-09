# GATTACA

[![License](https://img.shields.io/github/license/Feat-FeAR/GATTACA)](https://choosealicense.com/licenses/mit/)
[![Issues](https://img.shields.io/github/issues/Feat-FeAR/GATTACA)](https://github.com/Feat-FeAR/GATTACA/issues)
[![R Version](https://img.shields.io/badge/R%20version-4.1.1-informational)](https://www.r-project.org/)
[![Latest docker image size](https://img.shields.io/docker/image-size/mrhedmad/gattaca/latest?label=Docker%20image%20size)](https://hub.docker.com/repository/docker/mrhedmad/gattaca#)


<b>GATTACA</b> is a short for <b>General Algorithm for The Transcriptional Analysis by one-Channel Arrays</b>.
It consists of a set of R scripts for the analysis of Transportome Transcription in Cancer Datasets.
As the name suggests, GATTACA was originally written with (one-color/high-density) microarray technology in mind, however it can be easily adapted to analyze RNA-Seq data as well.

GATTACA assumes data to be already background-subtracted, log2-transformed, and interarray-normalized.
To this purpose, use the platform-specific <i>*_to_Expression.R</i> script.

The available docker containers live [here on Docker Hub](https://hub.docker.com/repository/docker/mrhedmad/gattaca#).

> **Disclaimer**: These instructions are given with a unix-like environment in mind. The program was developed and tested in an Ubuntu installation. See issue #3.

## Installation

These scripts are bundled with their dependencies in a Docker container to make usage as painless and reproducible as possible. To install, follow these instructions:

1. Install docker. A tutorial can be found [at the official docker website](https://docs.docker.com/get-docker/).
2. Install [`git`](https://git-scm.com/).
3. Pull down the repository using:
    ```bash
    git clone https://github.com/Feat-FeAR/GATTACA.git
    ```
    This clones the repository by default in the current working directory under a new folder `./GATTACA`.
4. Run `.GATTACA/GATTACA -h` to begin. The `GATTACA/GATTACA` script is the access point to all the other scripts in the repo. Read the usage section in this readme to learn more. The script takes care of downloading the containers and running them.

### Troubleshooting
If the `GATTACA` script cannot be run due to permissions faults, run `chmod +x GATTACA` to give it run permissions, or give it to bash directly using `/bin/bash GATTACA`.

For more help, [open an issue](https://github.com/Feat-FeAR/GATTACA/issues/new).

## Usage
Using the scripts is as simple as starting the `GATTACA` bash script. The various help messages (accessed using `-h`) detail the usage of each subcommand.

These are the currently available subcommands:
- `GATTACA init`: Used to make configuration files for `GATTACA run`.
- `GATTACA run`: Runs the `GATTACA.R` script. 

### Script: GATTACA.R
The `GATTACA.R` script is the main program of this collection. It performs differential gene expression analysis from single-channel microarray data.

To access the script, use the `GATTACA/GATTACA run` function. Append `-h` for usage information and available runtime options (not listed here).

#### Options
To run GATTACA, you must provide an options file encoded in [yaml](https://yaml.org/). You can create a default version of the options by running `GATTACA init <path>`. Here are the supported options:

- General options:
    - `slowmode` (bool): If set to `true`, stops the execution of the script at various points of the analysis, allowing you to stop the analysis if you detect something has gone terribly wrong. If set to `false`, these checkpoints are silently ignored.
    - `show_data_snippets` (bool): If set to `true`, prints out small snippets of data at and around `slowmode` checkpoints. This is useful to check if something has gone terribly wrong. Pretty useless when `slowmode` is `false`, as these data snippets are not logged.
    - `log_path` (str | `null`): If set to `null`, spawn a logfile in the output folder named `gattaca_(current time).log`. Otherwise, a string with the filename (that will be overwritten) of the log file, always created in the output folder.
    - `save_pdf` (bool): If `true`, saves plots to `pdf` files in a subfolder in the output directory. If neither `save_pdf` nor `save_png` are set to `true`, the script will not save any plots. Usually only one of the two are set.
    - `save_png` (bool): If `true`, saves plots to `png` files in a subfolder in the output directory. If neither `save_pdf` nor `save_png` are set to `true`, the script will not save any plots. Usually only one of the two are set.
    - `plot_width` (int) and `plot_height` (int): Plot sizes in inches for `save_pdf` and pixels for `save_png`. 
    - `append_annotation` (bool): ?????
- Switches:
    - `dryrun` (bool): Run the analysis, but don't save any output file. Useful to test the analysis before committing to it. Currently not implemented. See issue #4.
    - `log2_transform` (bool): If the input isn't log-transformed, setting this to true will perform log2 transformation for you.
    - `renormalize` (bool): If set to true, runs quantile normalization after performing Robust Multichip Average (RMA) normalization.
    - Limma options:
        - `run_DEA` (bool): If true, detects Differentially Expressed Genes (DEGs) using the `limma` package (Performs Differential Expression Analysis). This performs better than `rankProd` for well-behaved data.
        - `run_paired_DEA` (bool): If false, runs the DEA in paired mode. Paired mode confronts specifically two samples. Note that to run in paired mode, you must select just a single contrast (see the design section).
    - Rankproduct options:
        - `run_DEA` (bool): If true, detects Differentially Expressed Genes (DEGs) using the `rankproduct` package (Performs Differential Expression Analysis). This performs better than `rankProd` for well-behaved data.
        - `run_paired_DEA` (bool): If false, runs the DEA in paired mode. Paired mode confronts specifically two samples. Note that to run in paired mode, you must select just a single contrast (see the design section).
- Chip options:
    - `type` (str): Define the version of the chip. The name must be given as one of these strings **exactly**. The chip type is used to determine which annotation database to use. If your chip isn't listed here, [open an issue](https://github.com/Feat-FeAR/GATTACA/issues/new) to request adding it. The currently supported chips are as follows:
      - Affymetrix Human Genome U133 Set (A)
      - Affymetrix Human Genome U133 Set (B)
      - Affymetrix Human Genome HG-U133 Plus 2.0 Array
      - Affymetrix Human Gene 1.0-ST Array
      - Agilent-026652 Whole Human Genome Microarray 4x44K v2
    - `use_remote` (bool): If true, the tool will use the remote databases to fetch annotation data. If false, the tool will use a local database instead, as defined in the `local_database_path` option.
    - `local_database_path` (str): See `use_remote`.
- Experimental design:
    - `experimental_groups` (list of str): A list with the (arbitrary) names of the experimental groups of the experiment. Note that the first group will be treated as the control/background group (regardless of its name).
    - `experimental_design` (list of int): A list defining the experimental matrix. This defines at what experimental groups the various samples belong to. Assign a number representing at which experimental group each sample (in the order it appears in the input file) belongs to.
    For example: in my experiment I have 9 samples: three are controls, three are treated with treatment 1 and three are treated with treatment 2, and they appear in the input file in the same order. I would set the `experimental_groups` to be `["control", "treat1", "treat2"]` and the `experimental_design` to be `[1, 1, 1, 2, 2, 2, 3, 3, 3]`.
    If the samples are many (more in each category than the number of experimental groups themselves), one can give as `experimental_design` just the number of samples in each group (again, in order). For example: in my experiment I have 75 samples, 25 in each of the categories of the previous example, my `experimental_design` could be `[25, 25, 25]`.
    - `group_colors` (list of str): A list of strings representing names of colors to be used in some plots. The list must be as long or longer than the number of `experimental_groups`.
    - `contrasts` (list of int): This list defines which group pairs are tested against each other to find differentially expressed genes. The pairs are generated from the experimental groups as non-palindrome ordered combinations. If the `experimental_groups` are `["control", "treat1", "treat2"]` then the generated pairs are `["treat1-control", "treat2-control", "treat2-treat1"]` A `slowmode` `dryrun` is especially useful to check if this option was set correctly.
    - Filters applied to the dataset:
        - `log2_expression` (num): Filter out (mark as non-differentally-expressed) any genes that have lower log2 expression than this value. This gets rid of genes that might be detected as differentially expressed just because their low expression fluctuates a lot between samples. This value depends a lot on the experiment design. For Agilent, check the un-hybridized `(-)3xSLv1 NegativeControl` probe as a plausible value (a way to check this will be implemented in the future (maybe)). For  Affymetrix, 4 is a good baseline.
        - `fold_change` (num): Fold Change Threshold. Usually, `|log2FC| > 0.5` or `|log2FC| > 1`. The reasoning behind this filter is the same as `log2_expression`.
        - `min_groupwise_presence` (num): The percentage (rouded up) of samples in each group that a gene has to have over `fold_change` expression to have the possibility to be included as as differentially expressed gene. This is done to be sure that the differentially expressed genes are represented in most analysed samples and therefore robust from the biological point of view.

### Script: Affy2Expression
This script converts a collection of `.CEL` files from Affymetrix chips to an expression matrix.
It uses the `oligo` package that can automatically detect and download a wide variety of annotation databases for chips.

The script is called with `GATTACA prepaffy` command. The options are straightforward to understand (use `GATTACA prepaffy -h` to see them).

The `.CEL` files in the input folder are found by looking at their file extension, and are not found recursively.

## Building the image
If you wish to rebuild the Docker container, simply clone the repo, and run `docker build .`. More information on the [docker build documentation page](https://docs.docker.com/engine/reference/commandline/build/).
Note that rebuilding the image can take a long, long time due to the need to install required packages (that need to be compiled).
