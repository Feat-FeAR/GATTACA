# This is a remider on how to rebuild the GATTACA_manual.pdf

set -e
# First, run Xelatex to get a preliminary build
echo -n "XeLaTeX - Pass 1/3..."
xelatex main.tex

# 2 Then, run biber on the resulting files, to build biber build files:
echo -n "Generating bibliograpy with biber..."
biber main > /dev/null
echo "  Done."

# 3 Rerun Xelatex
echo -n "XeLaTeX - Pass 2/3..."
xelatex main.tex > /dev/null
echo "  Done."

# 4 Rerun it again to finalize bibliographic refs
echo -n "XeLaTeX - Pass 3/3..."
xelatex main.tex > /dev/null
echo "  Done."

# 5 Finally, move the result over
echo -n "Move rebuilt file..."
mv ./main.pdf ../GATTACA_Manual.pdf
echo "  Done."

echo "Removing junk..."
rm main.acn main.aux main.bbl main.bcf main.blg main.glo main.ist main.out main.run.xml main.toc
