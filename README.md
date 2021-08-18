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
    libfribidi-dev
```
For other distros, google the names of the respective libraries for your specific OS.

Then, launch the `install_dep` script:
```bash
./install_dep
```
This installs `remotes`, `here` and `renv` packages in your global environment, then installs the needed dependencies. Note that
the installation needs to compile the packages, so it could take a very long time.

Once the installation is done, the script can be run any numbers of times without reinstalling.

