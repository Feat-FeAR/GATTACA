#!/usr/bin/env bash

echo # This random echo is to give a bit more space to the help commands.

function _gattaca_print_general_help {
    echo "Usage: GATTACA [-v | --version <version>] [-h | --help] [ -l | --log_name <name> ]"
    echo "       [--verbose] <command> [<args>]"
    echo
    echo "Options:"
    echo "    -v | --version <version>      Specify the version of the docker to run <command> in. MUST be specified first."
    echo "    -l | --log_name <name>        Specify the filename of the log, created alongside the"
    echo "                                  output of many commands."
    echo "    --verbose                     Passing this flag includes debug information in the"
    echo "                                  logfiles."
    echo
    echo "Available commands:"
    echo "    init              Create a configuration file for 'GATTACA run'."
    echo "    run               Perform a differential expression analysis."
    echo "    prepaffy          Create an expression matrix from Affymetrix .CEL files."
    echo "    prepagil          Create an expression matrix from Agilent .txt files."
    echo "    annotate          Annotate an expression file with remote annotations."
    echo "    versions          Get a list of R packages installed in the docker and their versions."
    echo "    test              Run the tests associated with the code."
    echo "    parse             Return the parsed-and-expanded preview of an experimental design."
    echo
    echo "Use '-h' or '--help' in any subcommand to get additional help and possible options."
}


function _gattaca_help_init {
    echo "Usage: GATTACA init [-h | --help] <file_path>"
    echo
    echo "Create a configuration file for 'GATTACA run' at <file_path>."
    echo
    echo "NOTE: This command does not produce a log file."
}


function _gattaca_help_run {
    echo "Usage: GATTACA run [-h | --help] <output_dir> <input_file> <options_file>"
    echo
    echo "Perform differential gene analysis on an expression dataset '<input_file>', saving"
    echo "the output plots and files in '<output_dir>'. The runtime options for this command are"
    echo "specified in the <options_file>, encoded in yaml. See the 'init' command and the"
    echo "GATTACA repository README for more information about a specific container version."
    echo
    echo "NOTE: The options file will be copied in the same folder as the input file, due"
    echo "to (self-imposed) limitations in the docker engine."
}


function _gattaca_help_prepaffy {
    echo "Usage: GATTACA prepaffy [-h | --help] [-r | --remove-controls] [--plot-number <int>]"
    echo "       [--plot-size <x,y>] [--png] <input_dir> <output_file>"
    echo
    echo "Run preprocessing on all .CEL files present in <input_dir> to get a CSV expression matrix,"
    echo "saved as <output_file>."
    echo "Diagnostic plots will be saved alongside the output file, in a folder named 'prepaffy Figures'."
    echo
    echo "Options:"
    echo "  -r | --remove-controls      If set, removes control probes from the expression data."
    echo "  --plot-number <int>         Print only top <int> plots, sorted by distortion score."
    echo "  --plot-size <x,y>           Plot sizes (in inches) to use for the diagnostic PDF"
    echo "                              plots. Defaults to '12,5'. Note: Using the '--png'"
    echo "                              flag produces PNG files at a resolution of 300 ppi."
    echo "  --png                       Produce PNG files instead of PDF files."
}


function _gattaca_help_annotate {
    echo "Usage: GATTACA annotate [-h | --help] [-d | --database <database_name>] <input_file>"
    echo "       <output_file> <chip_id>"
    echo
    echo "Annotate an expression dataset '<input_file>' with extra annotations."
    echo
    echo "The default annotations cover most human chips. If, instead, you need a"
    echo "more specific annotation package, use the flag '-d | --database <database_name>'"
    echo "with the name of the package with the annotations you wish to use from"
    echo "bioconductor."
    echo
    echo "The annotations include the Gene Symbol (SYMBOL), the gene name (GENENAME),"
    echo "the collection of Ensembl IDs (ENSEMBL), the package(s) used to source the"
    echo "annotations (package_name) and its version(s) (version)."
    echo "It is currently impossible to change which information is annotated."
    echo
    echo "Save the output in <output_file>."
    echo
    echo "Options:"
    echo "  -d | --database <database_name>         Specify the name of the package to"
    echo "                                          annotate with, overriding the internal"
    echo "                                          annotation data. The package will be"
    echo "                                          installed on-the-fly."
}


function _gattaca_help_prepagil {
    echo "Usage: GATTACA prepagil [-h | --help] [-r | --remove-controls] [-p | --grep-pattern <pattern>]"
    echo "       [--plot-number <int>] [--plot-size <x,y>] [--png] <input_dir> <output_file>"
    echo
    echo "This command prepares normalized, log2 transformed expression matrices from Agilent files."
    echo "Runs on all files in <input_dir> that end in a '.txt' (or <pattern>) to get a CSV expression"
    echo "matrix, saved as <output_file>."
    echo "Diagnostic plots will be saved alongside the output file, in a folder named 'prepagil Figures'."
    echo
    echo "Options:"
    echo "  -r | --remove-controls            If set, removes control probes from the expression data."
    echo "  -p | --grep-pattern <pattern>     The RegEx pattern to match input files with."
    echo "  --plot-number <int>               Print only top <int> plots, sorted by distortion score."
    echo "  --plot-size <x,y>                 Plot sizes (in inches) to use for the diagnostic PDF"
    echo "                                    plots. Defaults to '12,5'. Note: Using the '--png'"
    echo "                                    flag produces PNG files at a resolution of 300 ppi."
    echo "  --png                             Produce PNG files instead of PDF files."
}

function _gattaca_help_parse {
    echo "Usage: GATTACA parse [-h | --help] \"<unparsed>\""
    echo
    echo "Return the parsed-and-expanded preview of <unparsed> experimental design."
    echo
    echo "NOTE: The <unparsed> expression needs to be enclosed in quotation marks."
    echo "This command does not produce a log file."
}

function _gattaca_assert_is_file {
    # Assert if a certain path is pointing to a file that exists
    if test ! -f "$1"; then
        >&2 echo "Error: '$1' is not a file."
        exit 1
    fi
}


function _gattaca_assert_is_directory {
    # Assert if a certain path is pointing to a directory that exists
    if test ! -d "$1"; then
        >&2 echo "Error: '$1' is not a directory."
        exit 1
    fi
}


function _gattaca_run_init {
    target_path="."

    while test $# -gt 0; do
        case "$1" in
            -h | --help)
                _gattaca_help_init
                exit 0
            ;;
            * )
                break
            ;;
        esac
    done

    if test $# -eq 1; then
        target_path="$1"
    elif test $# -eq 0; then
        >&2 echo "Missing required argument '<file_path>'"
        exit 1
    else
        >&2 echo "Wrong number of arguments: ${#} (expected 1)"
        exit 1
    fi

    target_mountpoint="$(dirname "$target_path")"
    target_name=$(basename -- "$target_path")

    mkdir -p "$target_mountpoint"

    echo "Spooling up docker instance..."

    docker run \
        --rm \
        --mount type=bind,source="$target_mountpoint",target=/GATTACA/target \
        cmalabscience/gattaca:"$version" \
        "unnamed", "INFO" \
        "init" "$target_name"

    echo "Created options file in '$target_path'"
}

function _gattaca_run_run {

    while test $# -gt 0; do
        case "$1" in
            -h | --help)
                _gattaca_help_run
                exit 0
            ;;
            * )
                break
            ;;
        esac
    done

    # <output_dir> <input_file> <options_file>
    if test $# -eq 3; then
        target_mountpoint=$(realpath "$1")
        input_path=$(realpath "$2")
        options_path=$(realpath "$3")
    else
        >&2 echo "Wrong number of arguments: ${#} (expected 3)"
        exit 1
    fi

    _gattaca_assert_is_file "$input_path"
    _gattaca_assert_is_file "$options_path"

    mkdir -p "$target_mountpoint"

    input_filename=$(basename -- "$input_path")
    input_mountpoint="$(dirname "$input_path")"
    options_filename=$(basename -- "$options_path")
    options_mountpoint="$(dirname "$options_path")"

    if [ "${input_filename}" = "${options_filename}" ]; then
        >&2 echo "Input and option files cannot have the same name. Sorry!"
        exit 1
    fi

    if [ ! "$input_mountpoint" = "$options_mountpoint" ]; then
        echo "Copying options file to input dir..."
        cp "$options_path" "$input_mountpoint"
    fi

    echo "Spooling up docker instance..."

    docker run \
        -it --rm \
        --mount type=bind,source="$target_mountpoint",target=/GATTACA/target \
        --mount type=bind,source="$input_mountpoint",target=/GATTACA/input,readonly \
        cmalabscience/gattaca:"$version" \
        "$log_name" "$log_level" \
        "gattaca" "$input_filename" "$options_filename"

    if [ ! "$input_mountpoint" = "$options_mountpoint" ]; then
        echo "Removing copy of options file from input dir..."
        rm "$input_mountpoint/$options_filename"
    fi
}

function _gattaca_run_prepaffy {
    remove_controls="FALSE"
    plot_size="12,5"
    use_pdf="TRUE"
    plot_number="Inf"

    while test $# -gt 0; do
        case "$1" in
            -h | --help)
                _gattaca_help_prepaffy
                exit 0
            ;;
            -r | --remove-controls)
                shift
                remove_controls="TRUE"
            ;;
            --plot-number)
              shift
              plot_number="$1"
              shift
            ;;
            --plot-size)
              shift
              plot_size=$1
              shift
            ;;
            --png)
              shift
              use_pdf="FALSE"
              ;;
            * )
                break
            ;;
        esac
    done

    # <input_dir> <output_file>

    if test $# -eq 2; then
        input_mountpoint=$(realpath "$1")
        output_path=$(realpath "$2")
    else
        >&2 echo "Wrong number of arguments: ${#} (expected 2)"
        exit 1
    fi

    output_mountpoint=$(dirname "$output_path")
    output_filename=$(basename -- "$output_path")

    mkdir -p "$output_mountpoint"

    _gattaca_assert_is_directory "$input_mountpoint"

    echo "Spooling up docker instance..."
    docker run \
        -it --rm \
        --mount type=bind,source="$output_mountpoint",target=/GATTACA/target \
        --mount type=bind,source="$input_mountpoint",target=/GATTACA/input,readonly \
        cmalabscience/gattaca:"$version" \
        "$log_name" "$log_level" \
        "prepaffy" "$output_filename" "$remove_controls" "$plot_size" "$use_pdf" "$plot_number"
}

function _gattaca_run_annotate {
    database="NA"
    while test $# -gt 0; do
        case "$1" in
            -h | --help)
                _gattaca_help_annotate
                exit 0
            ;;
            -d | --database)
                shift
                database="$1"
                shift
            ;;
            * )
                break
            ;;
        esac
    done

    # <input_file> <output_file> <database_name>

    if test $# -eq 2; then
        input_path=$(realpath "$1")
        target_path="$2"
    else
        >&2 echo "Wrong number of arguments: ${#} (expected 2)"
        exit 1
    fi

    input_mountpoint="$(dirname "$input_path")"
    target_mountpoint="$(dirname "$target_path")"

    output_filename=$(basename "$target_path")
    input_filename=$(basename "$input_path")

    _gattaca_assert_is_file "$input_path"

    mkdir -p "$target_mountpoint"

    echo "Spooling up docker instance..."
    docker run \
        -it --rm \
        --mount type=bind,source="$target_mountpoint",target=/GATTACA/target \
        --mount type=bind,source="$input_mountpoint",target=/GATTACA/input,readonly \
        cmalabscience/gattaca:"$version" \
        "$log_name" "$log_level" \
        "annotate" "$input_filename" "$output_filename" "$database"
}

function _gattaca_run_prepagil {
    grep_pattern="*.(txt|TXT)"
    remove_controls="FALSE"
    plot_size="12,5"
    use_pdf="TRUE"
    plot_number="Inf"

     while test $# -gt 0; do
        case "$1" in
            -h | --help)
                _gattaca_help_prepagil
                exit 0
            ;;
            -r | --remove-controls)
                shift
                remove_controls="TRUE"
            ;;
            -p | --grep-pattern)
                shift
                grep_pattern="$1"
                shift
            ;;
            --plot-number)
                shift
                plot_number="$1"
                shift
            ;;
            --plot-size)
                shift
                plot_size="$1"
                shift
            ;;
            --png)
                shift
                use_pdf="FALSE"
            ;;
            * )
                break
            ;;
        esac
    done

    if test $# -eq 2; then
        input_dir=$(realpath "$1")
        target_path="$2"
    else
        >&2 echo "Wrong number of arguments: ${#} (expected 2)"
        exit 1
    fi

    target_mountpoint="$(dirname "$target_path")"
    output_filename=$(basename -- "$target_path")

    _gattaca_assert_is_directory "$input_dir"

    mkdir -p "$target_mountpoint"

    echo "Spooling up docker instance..."
    docker run \
        -it --rm \
        --mount type=bind,source="$target_mountpoint",target=/GATTACA/target \
        --mount type=bind,source="$input_dir",target=/GATTACA/input,readonly \
        cmalabscience/gattaca:"$version" \
        "$log_name" "$log_level" \
        "prepagil" "$output_filename" "$grep_pattern" "$remove_controls" "$plot_size" "$use_pdf" "$plot_number"
}

function _gattaca_run_tests {
  docker run \
      --rm \
      cmalabscience/gattaca:"$version" \
      "NULL" "NULL" \
      "testall"
}

function _gattaca_run_versions {
  docker run \
      --rm \
      cmalabscience/gattaca:"$version" \
      "NULL" "NULL" \
      "versions"
}

function _gattaca_run_parser {

    while test $# -gt 0; do
        case "$1" in
            -h | --help)
                _gattaca_help_parse
                exit 0
            ;;
            * )
                break
            ;;
        esac
    done

    if test $# -eq 1; then
        design="$1"
    elif test $# -eq 0; then
        >&2 echo "Missing required argument '<unparsed>'"
        exit 1
    else
        >&2 echo "Wrong number of arguments: ${#} (expected 1)"
        exit 1
    fi

    docker run \
        --rm \
        cmalabscience/gattaca:"$version" \
        "NULL" "NULL" \
        "parser" "$design"
}

if [ $# == 0 ]; then
    _gattaca_print_general_help
fi

version="latest"
log_level="INFO"

while test $# -gt 0; do
    case "$1" in
        -h | --help)
            _gattaca_print_general_help
            exit 0
        ;;
        -v | --version)
            shift # Get rid of `-v`
            version="$1"
            ## TODO: Print out warnings for very old docker images that do
            # not support some commands.
            shift
        ;;
        -l | --log-name)
            shift
            log_name="$1"
            shift
        ;;
        --verbose)
            shift
            log_level="DEBUG"
        ;;
        init)
            shift
            _gattaca_run_init "$@"
            exit 0
        ;;
        run)
            shift
            _gattaca_run_run "$@"
            exit 0
        ;;
        prepaffy)
            shift
            _gattaca_run_prepaffy "$@"
            exit 0
        ;;
        annotate)
            shift
            _gattaca_run_annotate "$@"
            exit 0
        ;;
        prepagil)
            shift
            _gattaca_run_prepagil "$@"
            exit 0
        ;;
        test)
            shift
            _gattaca_run_tests
            exit 0
        ;;
        versions)
            shift
            _gattaca_run_versions
            exit 0
        ;;
        parse)
            shift
            _gattaca_run_parser "$@"
            exit 0
        ;;
        *)
            echo "Invalid parameter '{$1}'."
            _gattaca_print_general_help
            exit 0
    esac
done
