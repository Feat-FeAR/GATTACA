# GATTACA

[![License](https://img.shields.io/github/license/Feat-FeAR/GATTACA)](https://choosealicense.com/licenses/mit/)
[![Issues](https://img.shields.io/github/issues/Feat-FeAR/GATTACA)](https://github.com/Feat-FeAR/GATTACA/issues)
[![R Version](https://img.shields.io/badge/R%20version-4.1.1-informational)](https://www.r-project.org/)
[![Latest docker image size](https://img.shields.io/docker/image-size/mrhedmad/gattaca/v0.0.0-beta.2?label=Docker%20image%20size)](https://hub.docker.com/repository/docker/mrhedmad/gattaca#)


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

### Command: `GATTACA run`
The `GATTACA run` command performs differential gene expression analysis from single-channel microarray data.

Append `-h` for usage information and available runtime options (not listed here).

#### YAML Options
To run GATTACA, you must provide an options file encoded in [yaml](https://yaml.org/). You can create a default version of the options by running `GATTACA init <path>`. Here are the supported options:

- General options:
    - `slowmode` (bool): If set to `true`, stops the execution of the script at various points of the analysis, allowing you to stop the analysis if you detect something has gone terribly wrong. If set to `false`, these checkpoints are silently ignored.
    - `show_data_snippets` (bool): If set to `true`, prints out small snippets of data at and around `slowmode` checkpoints. This is useful to check if something has gone terribly wrong. Pretty useless when `slowmode` is `false`, as these data snippets are not logged.
    - `save_pdf` (bool): If `true`, saves plots to `pdf` files in a subfolder in the output directory. If neither `save_pdf` nor `save_png` are set to `true`, the script will not save any plots. Usually only one of the two are set.
    - `save_png` (bool): If `true`, saves plots to `png` files in a subfolder in the output directory. If neither `save_pdf` nor `save_png` are set to `true`, the script will not save any plots. Usually only one of the two are set.
    - `plot_width` (int) and `plot_height` (int): Plot sizes in inches.
    - `png_ppi` (int): The pixels per inch resolution for `png` plots.
    - `enumerate_plots` (bool): If `true`, the plots created by GATTACA will be enumerated following the order by which they are created. This allows easier inspection of them.
    - `annotation_chip_id` (str or null): If annotation are added, some plots (like volcano plots) are more legible as they use HUGO symbols instead of probe IDs. Specify here one available annotation ID based on your chip. Leave to `null` to not annotate. Available annotations (chip_name -> `id`):
        - Affymetrix Human Genome U133 Set (A) -> `hgu133a`
        - Affymetrix Human Genome U133 Set (B) -> `hgu133b`
        - Affymetrix Human Genome HG-U133 Plus 2.0 Array -> `hgu133plus2`
        - Agilent-026652 Whole Human Genome Microarray 4x44K v2 -> `HsAgilentDesign026652`
        - Affymetrix Human Gene 1.0-ST Array -> `hugene10st`
- Switches:
    - `dryrun` (bool): Run the analysis, but don't save any output file. Useful to test the analysis before committing to it.
    - `log2_transform` (bool): If the input isn't log-transformed, setting this to true will perform log2 transformation for you.
    - `renormalize` (bool): If set to true, runs quantile normalization after performing Robust Multichip Average (RMA) normalization.
    - `limma` (bool): If true, detects Differentially Expressed Genes (DEGs) using the `limma` package (Performs Differential Expression Analysis). This performs better than `rankProd` for well-behaved data.
    - `rankproduct` options: If true, detects Differentially Expressed Genes (DEGs) using the `rankproduct` package (Performs Differential Expression Analysis). This performs better than `rankProd` for well-behaved data.
- Design options:
    - `experimental_design` (str): A string representing the experimental design. See the special section at the end of this list.
    - `group_colors` (list of str): A list of strings representing names of colors to be used in some plots. The list must be as long or longer than the number of `experimental_groups`.
    - `contrasts` (list of int): This is a list of strings. The strings represent the contrasts to check, in the format "difference-baseline". Each group specified must also be a group in the experimental design.
    - Filters applied to the dataset:
        - `log2_expression` (num): Filter out (mark as non-differentally-expressed) any genes that have lower log2 expression than this value. This gets rid of genes that might be detected as differentially expressed just because their low expression fluctuates a lot between samples. This value depends a lot on the experiment design. For Agilent, check the un-hybridized `(-)3xSLv1 NegativeControl` probe as a plausible value (a way to check this will be implemented in the future (maybe)). For  Affymetrix, 4 is a good baseline.
        - `fold_change` (num): Fold Change Threshold. Usually, `|log2FC| > 0.5` or `|log2FC| > 1`. The reasoning behind this filter is the same as `log2_expression`.
        - `min_groupwise_presence` (num): The percentage (rounded up) of samples in each group that a gene has to have over `fold_change` expression to have the possibility to be included as as differentially expressed gene. This is done to be sure that the differentially expressed genes are represented in most analysed samples and therefore robust from the biological point of view.
        
##### Setting the experimental design string
Each sample is marked with a letter (or series of letters) representing conditions. They can also be marked with arbitrary numbers representing sample pairs for paired designs. The letter and optional numbers are known as labels, and are separated by commas (all spaces are ignored).

The parser respects already expanded labels and will retain them as-is:
```r
> design_parser("a, b, c") # Unpaired
  [1] "a, b, c"
> design_parser("a1, b2, c3") # Paired
  [1] "a1, b2, c3"
```
To avoid high repetition, parsing supports two types of pattern expansion. Round brackets `()` represent intercalated expansions, where labels inside the brackets are repeated in the order inside the brackets for the number of times specified:
```r
> design_parser("(a, b):3") # Unpaired 'intercalated'
  [1] "a, b, a, b, a, b"
```
Square brackets represent ordered expansions, where each label in the brackets is repeated on its own the number of times specified:
```r
> design_parser("[a, b]:3") # Unpaired 'ordered'
  [1] "a, a, a, b, b, b"
> design_parser("[a, b, c]:3") # Unpaired 'ordered'
  [1] "a, a, a, b, b, b, c, c, c"
```
The * wildcard inside the brackets will be replaced with unique numbers, so that expansion is actually useful for paired designs:
```r
> design_parser("(a*, b*):3") # Paired 'intercalated'
  [1] "a1, b1, a2, b2, a3, b3"
> design_parser("[a*, b*]:3") # Paired 'ordered'
  [1] "a1, a2, a3, b1, b2, b3"
> design_parser("[a*, b*, c*]:2") # Paired 'ordered'
  [1] "a1, a2, b1, b2, c1, c2"
```
The wildcard is assured to not collide with any already used numbers in the pattern, starting one number after the largest number found:
```r
> design_parser("a3, b7, (a*, b*):3") # Paired 'intercalated'
  [1] "a3, b7, a8, b8, a9, b9, a10, b10"
> design_parser("a3, b7, [a*, b*]:3") # Paired 'ordered'
  [1] "a3, b7, a8, a9, a10, b8, b9, b10"
```
The various patterns can be mixed and matched:
```r
> design_parser("(a*, b*):3, a3, b6, [a*]:2")
  [1] "a7, b7, a8, b8, a9, b9, a3, b6, a10, a11"
```  
Note: Both types of expansion work identically if only one label is specified in the brackets.


### Command: `prepaffy`
Convert a collection of `.CEL` files from Affymetrix chips to an expression matrix. The command uses the `oligo` package that can automatically detect and download a wide variety of annotation databases for chips.

The options are straightforward to understand (use `GATTACA prepaffy -h` to see them).

The `.CEL` files in the input folder are found by looking at their file extension, and are not found recursively.


### Command: `prepagil`
Convert a collection of files from a variety of scanners for Agilent microarrays to an expression matrix. The command uses `limma` to parse the input files. As the format of these files is varied, they are found by `grep`-ping files in the target folder. See `./GATTACA prepaffy --help` for more information.


### Command: `annotate`
Annotate data created by these commands with annotations from a variety of databases present in bioconductor. The databases are downloaded at runtime to have the latest annotations. All bioconductor databases support the following fields: `ACCNUM`, `CHR`, `CHRLOC`, `CHRLOCEND`, `ENSEMBL`, `ENTREZID`, `ENZYME`, `GENENAME`, `GO`, `MAP`, `OMIM`, `PATH`, `PMID`, `REFSEQ`, `SYMBOL` and `UNIPROT`.

The options are straightforward to understand (use `GATTACA annotate -h` to see them).

Available annotations (chip_name -> `id`):
  - Affymetrix Human Genome U133 Set (A) -> `hgu133a`
  - Affymetrix Human Genome U133 Set (B) -> `hgu133b`
  - Affymetrix Human Genome HG-U133 Plus 2.0 Array -> `hgu133plus2`
  - Agilent-026652 Whole Human Genome Microarray 4x44K v2 -> `HsAgilentDesign026652`
  - Affymetrix Human Gene 1.0-ST Array -> `hugene10st`

## Building the image
If you wish to rebuild the Docker container, simply clone the repo, and run `docker build .`. More information on the [docker build documentation page](https://docs.docker.com/engine/reference/commandline/build/).
Note that rebuilding the image can take a long, long time due to the need to install required packages (that need to be compiled).
