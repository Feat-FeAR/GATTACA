#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

print(args)

# Args here are in this order:
# > UUID GUID command logname loglevel_console loglevel_file ...
# The ... are the args to pass on the command/module.

# Restore UUIDs and GUIDs on exit
# TODO: Test if this is possible to do with the docker run -u option,
# and an internal user.
.restore_uids <- function(UUID, GUID) {
    cat(paste0("Restoring files to UUID/GUID ", UUID, ":", GUID, "\n"))
    exclude_files <- c(
        list.files("/GATTACA/target", all.files = TRUE, recursive = TRUE),
        list.files("/GATTACA/log", all.files = TRUE, recursive = TRUE),
        list.files("/GATTACA/input", all.files = TRUE, recursive = TRUE)
    )

    restorer <- function(...) {
        current_files <- c(
            list.files("/GATTACA/target", all.files = TRUE, recursive = TRUE),
            list.files("/GATTACA/log", all.files = TRUE, recursive = TRUE),
            list.files("/GATTACA/input", all.files = TRUE, recursive = TRUE)
        )
        
        new_files <- current_files[! current_files %in% exclude_files]
        if (length(new_files) == 0) {cat("No files to restore permissions to.\n"); return()}

        owner <- paste0(UUID, ":", GUID)
        cat("Restoring file permissions...\n")
        system2("xargs -0 chown", args = owner, stdin = new_files)
    }

    return(restorer)
}

UIDREST <- .restore_uids(args[1], args[2])

args <- args[-c(1, 2)]

COMMAND <- args[1]

if (COMMAND == "test") {
    # This is a testrun.
    # This means a few things:
    #  - The usual paths are not mounted (e.g. /GATTACA/target...)
    #  - All logging has to be redirected to stdout (e.g., just CAT, and not to file.)
    #  - No initial packages are loaded by the tests themselves. Use the
    #  setup/teardown files to load things instead.
    library(testthat)

    cat("Starting shared tests...\n")
    test_file("/GATTACA/shared/tests/test-shared.R")

    cat("Starting module tests...\n")
    modules <- c("prepagil", "prepaffy", "annotation", "analize")
    # Test all modules
    for (path in paste0("/GATTACA/modules/", modules, "/tests/test-module.R")) {
        test_file(path)
    }

    cat("Done executing tests.\n")
    quit(save = "no", status = 0)
}

invisible(reg.finalizer(environment(), UIDREST, onexit = TRUE))

args <- args[-1]

# Args here are in this order:
# > logname loglevel_console loglevel_file ...
# The ... are the args to pass on the command/module.

# Load the base packages
suppressPackageStartupMessages({
    library(tidyverse)
    library(progress)
})

source("/GATTACA/shared/tools.R")

if (length(args) < 4) {
    stop(paste0("Not enough arguments: ", length(args)))
}

# Setup logging facilities
log$set_name(args[1])
log$set_console_level(args[2])
log$set_file_level(args[3])

log$debug("Logfile initialized as '", args[1], "'")

COMMAND_ARGS <- args[-c(1:3)]

run_module <- function(module_name, module_args, exit_immediately = FALSE) {
    log$info("Saving environment...")
    .PREVIOUS_LS <- ls(all.names = TRUE)
    .PREVIOUS_NAMESPACES <- loadedNamespaces()
    .PREVOUS_OPTIONS <- options()

    log$info("Setting module arguments (", paste(module_args, collapse = ", "), ")")
    options(module.args = module_args)

    # Run the module
    log$info("Running module '", module_name, "'...")
    source(paste0("/GATTACA/modules/", module_name, "/entrypoint.R"))

    if (exit_immediately) {
        quit(save = "no", status = 0)
    }

    # Restore the environment
    log$info("Restoring namespaces...")
    .NEW_NAMESPACES <- loadedNamespaces()
    .TO_UNLOAD <- .NEW_NAMESPACES[! .NEW_NAMESPACES %in% .PREVIOUS_NAMESPACES]
    for (namespace in .TO_UNLOAD) {
        detach(namespace, character.only = TRUE)
    }

    log$info("Cleaning memory...")
    .NEW_LS <- ls(all.names = TRUE)
    .TO_REMOVE <- .NEW_LS[! .NEW_LS %in% .PREVIOUS_LS]
    rm(list = .TO_REMOVE)

    log$info("Resetting options...")
    .NEW_OPTIONS <- modifyList(options(), NULL, keep.null = TRUE)
    options(modifyList(.NEW_OPTIONS, .PREVOUS_OPTIONS, keep.null = TRUE))
}

# Run the desired module
switch(
    COMMAND,
    prepaffy = {
        run_module("prepaffy", COMMAND_ARGS, exit_immediately = TRUE)
    },
    prepagil = {
        run_module("prepagil", COMMAND_ARGS, exit_immediately = TRUE)
    },
    annotation = {
        run_module("annotation", COMMAND_ARGS, exit_immediately = TRUE)
    },
    analize = {
        run_module("analize", COMMAND_ARGS, exit_immediately = TRUE)
    }
)
