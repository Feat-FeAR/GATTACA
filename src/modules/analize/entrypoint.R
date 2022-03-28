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

# Load the arguments
args <- getOption("module.args")

# Test that the arguments are valid

# Load required libraries
module.packages <- c()
graceful_load(module.packages)

# Load the functions for this module 
source("/GATTACA/modules/thismodule/...")

# Run the module
