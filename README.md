# GATTACA

<b>GATTACA</b> is a short for <b>General Algorithm for The Transcriptional Analysis by one-Channel Arrays</b>.
It consists of a set of R scripts for the analysis of Transportome Transcription in Cancer Datasets.
As the name suggests, GATTACA was originally written with (one-color/high-density) microarray technology in mind, however it can be easily adapted to analyze RNA-Seq data as well.

GATTACA assumes data to be already background-subtracted, log2-transformed, and interarray-normalized.
To this purpose, use the platform-specific <i>*_to_Expression.R</i> script.

## Installation

To run these scripts, you will need to install R version `4.1.0` and the linked dependencies. To increase reproducibility, the project uses the `renv` package to manage its dependencies. 

First, install the required external dependencies. In ubuntu `21.04` use:
```bash
sudo apt install libgmp3-dev libmpfr-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev proj-bin libproj-dev libcairo2-dev libxt-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```
For other distros, google the names of the respective libraries for your specific OS.

Then, launch the `install_dep` script:
```bash
./install_dep
```
This installs `remotes`, `here` and `renv` packages in your global environment, then installs the needed dependencies. Note that
the installation needs to compile the packages, so it could take a very long time.

Once the installation is done, the script can be run any numbers of times without reinstalling.

## GATTACA Options

To run GATTACA, you must provide an options file encoded in (yaml)[https://yaml.org/]. Here are the supported options:

- General options:
    - `slowmode` (bool): If set to `true`, stops the execution of the script at various points of the analysis, allowing you to stop the analysis if you detect something has gone terribly wrong. If set to `false`, these checkpoints are silently ignored.
    - `show_data_snippets` (bool): If set to `true`, prints out small snippets of data at and around `slowmode` checkpoints. This is useful to check if something has gone terribly wrong. Pretty useless when `slowmode` is `false`, as these data snippets are not logged.
    - `log_path` (str | `null`): If set to `null`, spawn a logfile in the output folder named `gattaca_(current time).log`. Otherwise, a string with the absolute path to a file (that will be overwritten) to save the log to.
    - `save_pdf` (bool): If `true`, saves plots to `pdf` files in a subfolder in the output directory. If neither `save_pdf` nor `save_png` are set to `true`, the script will not save any plots. Usually only one of the two are set.
    - `save_png` (bool): If `true`, saves plots to `png` files in a subfolder in the output directory. If neither `save_pdf` nor `save_png` are set to `true`, the script will not save any plots. Usually only one of the two are set.
    - `plot_width` (int) and `plot_height` (int): Plot sizes in iches for `save_pdf` and pixels for `save_png`. 
    - `append_annotation` (bool): ?????
- Switches:
    - `dryrun` (bool): Run the analysis, but don't save any output file. Useful to test the analysis before committing to it.
    - `log2_transform` (bool): If the input isn't log-transformed, setting this to true will perform log2 transformation for you.
    - `renormalize` (bool): If set to true, runs quantile normalization after performing Robust Multichip Average (RMA) normalization.
    - Limma options:
        - `run_DEA`: If true, detects Differentially Expressed Genes (DEGs) using the `limma` package (Performs Differential Expression Analysis). This performs better than `rankProd` for well-behaved data.
        - `run_paired_DEA`: If false, runs the DEA in paired mode. Paired mode  


## SWITCHES - Various switches to turn parts of the analysis on or off
switches:
  # Run the analysis, but don't save any output file. Useful to test the
  # analysis before committing.
  dryrun: false # bool
  # Log2 transform the input intensities before running the analysis?
  log2_transform: false # bool
  # Should re-run normalization after performing RMA?
  renormalize: true # bool
  # Options related to the `limma` package
  limma:
    # Run Differential gene expression detection with `limma`
    # This uses parametric tests to find DEGs.
    run_DEA: true # bool
    # The package can run both paired and unpaired differental expression
    # analises. If `true`, run the paired analysis. If `false`, run the unpaired
    # analysis. Note:: Paired analysis can only be run when only one contrast is
    # selected. Check the README if you intend to use this option.
    run_paired_DEA: false # bool
  # Options related to the `rankproduct` package
  rankproduct:
    # Run Differential gene expression detection with `rankproduct`
    # This uses non-parametric tests to find DEGs.
    run_DEA: true # bool
    # The package can run both paired and unpaired differental expression
    # analises. If `true`, run the paired analysis. If `false`, run the unpaired
    # analysis. Note:: Paired analysis can only be run when only one contrast is
    # selected. Check the README if you intend to use this option.
    run_paired_DEA: false # bool


## CHIP - Information about the chips.
chip:
  # Type of the chip. This is used to source annotation data
  # Currently supported chips:
  # -- Affymetrix Human Genome U133 Set (A)
  # -- Affymetrix Human Genome U133 Set (B)
  # -- Affymetrix Human Genome HG-U133 Plus 2.0 Array
  # -- Affymetrix Human Gene 1.0-ST Array
  # -- Agilent-026652 Whole Human Genome Microarray 4x44K v2
  type: null # str
  # Should the tool use remote or local databases?
  # Note: setting this to false makes the `type` option obsolete.
  use_remote: true # bool
  # If a local database is used, where is it stored (absolute paths only)?
  # Read the readme for more info on the expected format of the database file
  local_database_path: null # str pathlike


# DESIGN - Experimental design
design:
  # Define any experimental groups here.
  # NOTE!! : The FIRST group MUST be the control/background/wild-type!
  experimental_groups: ["Control", "Treated"] # list of str
  # Define the experimental matrix.
  # The experimental matrix describes which samples are of which group, as
  # defined above. Read the README for a guide on how to define this.
  experimental_design: [1, 1, 2, 2] # list of int
  # Specify the colours that represent the groups above. Must cover all groups (or more)
  group_colors: ["cornflowerblue","firebrick3","olivedrab3","darkgoldenrod1","purple","magenta3"] # list of str
  # Select which experimental group pairs should be analised. Group pairs are
  # generated from the group names defined above. Use slowmode if uncertain.
  contrasts: [1] # list of int
  # Filters:
  filters:
    # Filter out (mark as non-differentally-expressed) any genes that have lower
    # log2 expression than this value. This value depends a lot on the experiment.
    # The readme has pointers on how to choose the correct value.
    log2_expression: 4 # num
    # Filter out any genes that have lower absolute Fold Change than this value.
    # Usually 0.5 or 1.
    fold_change: 0.5 # num
    # Minimum gene presence per group
    min_groupwise_presence: 0.8