#!/usr/bin/env bash
set -e

echo "Rebuilding annotations."
echo ""
echo "Assuming that we are inside /GATTACA/scripts/"

echo "STEP 1/6 - Scrape Bioconductor"
python ./scrape_bioc.py ./biocpackages.csv

echo "STEP 2/6 - Generate Annotations"
Rscript --vanilla ./generate_raw_annotations.R

echo "STEP 3/6 - Process annotations"
python ./process_annotation_data.py ./annotdataraw.csv ./annotdata.csv

echo "STEP 4/6 - Generate .RData"
Rscript --vanilla ./format_annotations.R

echo "STEP 5/6 - Move new RData to repo"
mv ./full_annotations.RData ../src/resources/full_annotations.RData

echo "STEP 6/6 - Cleanup temp files"
rm ./annotdataraw.csv
rm ./annotdata.csv
rm ./biocpackages.csv
rm ./sortedinput.csv

echo "DONE"
