#!/usr/bin/env bash

order=$3

if [ $order != "testall" ] && [ $order != "versions" ] && [ $order != "parser" ]; then
    # We need the host UID/GID to give permissions back when we're done.
    # This is a bit of a hack, but I cannot find a better way to do it.
    host_uid=$(stat -c "%u" /GATTACA/target)
    host_gid=$(stat -c "%g" /GATTACA/target)

    # We take out the log name to setup logging.
    log_name=$1
    log_level=$2
fi

shift 3 # This gets rid of the order too

case "$order" in
    init)
        cp "/GATTACA/src/resources/GATTACA_default_options.yaml" "/GATTACA/target/$1"
    ;;
    gattaca)
        # As this should only be run from the cli.sh script, we can be sure that
        # $1 is input_filename and $2 is options_filename, both in /GATTACA/input
        # We are also sure that we need to write to /GATTACA/target
        inputfile="/GATTACA/input/${1}"
        optionfile="/GATTACA/input/${2}"

        Rscript \
            --vanilla \
            -e "source('/GATTACA/src/__init__.R')" \
            -e "source('/GATTACA/src/GATTACA.R')" \
            -e "setup_file_logging('/GATTACA/target', '${log_name}')" \
            -e "logger::log_threshold(${log_level})" \
            -e "GATTACA('${optionfile}', '${inputfile}', '/GATTACA/target')"
        ;;
    prepaffy)
        # Again, we know the order of the inputs from GATTACA cli.sh:
        # "$output_filename" "$remove_controls" "$plot_size" "$use_pdf" "$plot_number"
        outputfile="/GATTACA/target/${1}"

        Rscript \
            --vanilla \
            -e "source('/GATTACA/src/__init__.R')" \
            -e "source('/GATTACA/src/toExpression/Affy_CEL_to_Expression.R')" \
            -e "setup_file_logging('/GATTACA/target', '${log_name}')" \
            -e "logger::log_threshold(${log_level})" \
            -e "strsplit('${3}', ',')[[1]] |> as.numeric() -> sizes" \
            -e "affy2expression('/GATTACA/input', '${outputfile}', remove.controls = ${2}, plot.width = sizes[1], plot.height = sizes[2], use.pdf = ${4}, n_plots = ${5})"
        ;;
    annotate)
        # "$input_filename" "$output_filename" "$database"
        outputfile="/GATTACA/target/${2}"
        inputfile="/GATTACA/input/${1}"

        Rscript \
            --vanilla \
            -e "source('/GATTACA/src/__init__.R')" \
            -e "source('/GATTACA/src/tools/annotations.R')" \
            -e "setup_file_logging('/GATTACA/target', '${log_name}')" \
            -e "logger::log_threshold(${log_level})" \
            -e "selections <- strsplit('${4}', ',')[[1]]" \
            -e "annotate_to_file(expression_data_path = '${inputfile}', output_path = '${outputfile}', database_name = '${3}')"
        ;;
    prepagil)
        # "$output_filename" "$grep_pattern" "$remove_controls" "$plot_size" "$use_pdf" "$plot_number"
        outputfile="/GATTACA/target/${1}"

        Rscript \
            --vanilla \
            -e "source('/GATTACA/src/__init__.R')" \
            -e "source('/GATTACA/src/toExpression/Agilent_TXT_to_Expression.R')" \
            -e "setup_file_logging('/GATTACA/target', '${log_name}')" \
            -e "logger::log_threshold(${log_level})" \
            -e "strsplit('${4}', ',')[[1]] |> as.numeric() -> sizes" \
            -e "agil2expression('/GATTACA/input', '${outputfile}', grep_pattern='${2}', remove_controls=${3}, plot.width = sizes[1], plot.height = sizes[2], use.pdf = ${5}, n_plots = ${6})"
        ;;
    testall)
        Rscript --vanilla \
            -e "source('/GATTACA/src/__init__.R')" \
            -e "source('/GATTACA/tests/test_all.R')"
            exit 0
        ;;
    versions)
        Rscript --vanilla \
            -e "ip = as.data.frame(installed.packages())[,c('Package', 'Version')]" \
            -e "rownames(ip) <- NULL" \
            -e "cat('List of installed packages:\n')" \
            -e "print(ip)"
            exit 0
        ;;
    parser)
        # In order to speed it up, here we Rscript just the essential elements to
        # run the parser, rather than sourcing the entire __init__.R as before. 
        # Here the only input from GATTACA cli.sh is "$design"
        Rscript --vanilla \
            -e "suppressMessages(library(logger))" \
            -e "ROOT <- '/GATTACA'" \
            -e "source('/GATTACA/src/tools/tools.R')" \
            -e "design_parser('${1}') -> design" \
            -e "cat('Parsed experimental design:\n', design, '\n')"
            exit 0
        ;;
    *)
        echo "Unrecognized order: '${order}'"
        exit 1
        ;;
esac

# If we get to here we have to do this:
chown -R $host_uid:$host_gid /GATTACA/target
