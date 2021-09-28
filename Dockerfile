FROM rocker/r-ver:4.1.0 AS dependabox

# Install necessary libs to compile R packages
RUN apt update && \ 
  apt install -y --no-install-recommends \
    libgmp3-dev libmpfr-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev proj-bin libproj-dev libcairo2-dev libxt-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libbz2-dev

WORKDIR /GATTACA

# Copy just the installation scripts 
COPY ./dock_install.R ./
COPY ./renv.lock ./

# Install the dependencies of the R scripts...
RUN Rscript --vanilla \
  -e "source('./dock_install.R')" \
  -e "renv::isolate()"

# Copy the rest of the files
COPY ./ ./ 

# Make the hardcoded mountpoints for the outputs and inputs
RUN mkdir ./target && mkdir ./input

# Setup the entrypoint
ENTRYPOINT [ "./src/entry" ]
