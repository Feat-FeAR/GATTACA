FROM rocker/r-ver:4.1.1

# Make the hardcoded mountpoints for the outputs and inputs
RUN mkdir ./target && mkdir ./input

# Install necessary libs to compile R packages
# Some packages are not available from rocker, so I install as many deps as
# possible here.
RUN apt update && \
  apt install -y --no-install-recommends \
    libgmp3-dev libmpfr-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev proj-bin libproj-dev libcairo2-dev libxt-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libbz2-dev

WORKDIR /GATTACA

# Disable threading for preprocessCore since it is bugged rigth now.
RUN echo "options(configure.args = list(preprocessCore = '--disable-threading'))" >> /usr/local/lib/R/etc/Rprofile.site

# Package layer
RUN install2.r -n $(nproc --all) \
    tidyverse logger progress yaml statmod VennDiagram gplots UpSetR \
    testthat BiocManager && \
    /usr/local/lib/R/site-library/littler/examples/installBioc.r \
    preprocessCore PCAtools limma RankProd oligo affycoretools EnhancedVolcano org.Hs.eg.db

# Annotation layer - the most likely to be changed
RUN /usr/local/lib/R/site-library/littler/examples/installBioc.r \
    hgu133a.db hgu133b.db hgu133plus2.db HsAgilentDesign026652.db hugene10sttranscriptcluster.db \
    pd.hg.u133a pd.hg.u133b pd.hg.u133.plus.2 pd.hugene.1.0.st.v1 pd.drosophila.2

# Copy the rest of the files
COPY ./ ./

# Setup the entrypoint
ENTRYPOINT [ "./src/entryswitch.sh" ]
