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
# output_file, grep_pattern, remove_controls, width, heigth, use.pdf, n_plots
# Pass "NULL" or NULL to use the defaults

defaults = list(
    output_file = NULL,
    grep_pattern = "*.(txt|TXT)",
    remove_controls = TRUE,
    plot.width = 16,
    plot.height = 9,
    use.pdf = TRUE,
    n_plots = Inf
)

fun_args <- validate_arguments(args, defaults)
print(fun_args)

# Add the hardcoded arguments
fun_args$input_dir <- "/GATTACA/input/"
fun_args$output_file <- paste0("/GATTACA/target/", fun_args$output_file)

# Load required libraries
module.packages <- c("limma", "oligo", "reshape2")
graceful_load(module.packages)

# Load the functions for this module 
source("/GATTACA/modules/prepagil/Agilent_TXT_to_Expression.R")

# Run the module
do.call("agil2expression", args = fun_args)
