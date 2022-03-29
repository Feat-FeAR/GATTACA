# This file, when sourced, must do the following:
# - Load necessary packages for the module
# - Handle the incoming command line argument(s)
# - Run the actual command.
#
# The environment is shaped as such:
# GATTACA/
#   src/
#       shared/
#           ...
#       modules/
#           ...
#   input/  < Mounted to the input folder
#   target/ < Mounted to the output folder
#   logs/   < Mounted to the logs folder
#
# This script has access to the following:
# > The tidyverse;
# > The `tools` library in src/shared/tools.R;
#
# The logging module is already sourced, so the log object is present in the
# environment. 
#
# The command line options passed by the user (the wrapper), are in the
# `module.args` global option.

# Load the arguments -----------------------------------------------------------
args <- getOption("module.args")

# Test that the arguments are valid --------------------------------------------
# I expect these arguments, in order:
# options.path, input.file
# Pass "NULL" or NULL to use the defaults. Setting NULL in the defaults
# signifies a required argument.

defaults = list(
    # Plot options
    use_pdf = TRUE,
    plot_width = 16,
    plot_height = 9,
    png_ppi = 250,
    enumerate_plots = TRUE
)

fun_args <- validate_arguments(args, defaults)

# Add the hardcoded arguments
fun_args$input.file <- paste0("/GATTACA/input/", fun_args$input.file)
fun_args$options.path <- paste0("/GATTACA/input/", fun_args$options.path)
fun_args$output.dir <- paste0("/GATTACA/target/")

# Load required libraries.
module.packages <- c(
    "limma", "rankprod", "reshape2", "EnhancedVolcano", "gplots", "UpSetR",
    "yaml", "statmod"
)
graceful_load(module.packages)

# Load the functions for this module 
source("/GATTACA/modules/analize/analysis_core.R")

# Run the module
# Set options for printPlots
options(
    scriptName = "...",
    save.PNG.plot = !fun_args$use_pdf,
    save.PDF.plot = fun_args$use_pdf,
    plot.width = fun_args$plot_width,
    plot.height = fun_args$plot_height,
    png_ppi = fun_args$png_resolution,
    enumerate.plots = fun_args$enumerate_plots
)
# Remove plot-related options
fun_args[c("use_pdf", "plot_width", "plot_heigth", "png_ppi", "enumerate_plots")] <- NULL
do.call("GATTACA", args = fun_args)
