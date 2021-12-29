# This is a remider on how to rebuild the GATTACA_manual.pdf

# First, run Xelatex to get a preliminary build
xelatex main.tex

# 2 Then, run biber on the resulting files, to build biber build files:
biber main

# 3 Rerun Xelatex
xelatex main.tex

# 4 Rerun it again to finalize bibliographic refs
xelatex main.tex

# 5 Finally, move the result over
mv ./main.pdf ../GATTACA_Manual.pdf