FROM rocker/r-ver:4.1.0

# Install necessary libs to compile R packages
RUN apt update
RUN apt install -y --no-install-recommends \
    libgmp3-dev libmpfr-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev proj-bin libproj-dev libcairo2-dev libxt-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

WORKDIR /GATTACA

# Copy just the installation scripts 
COPY ./src/install.R .
COPY ./src/install_dep .
COPY ./src/renv.lock .

# Install the dependencies of the R scripts...
RUN ./install_dep

# Copy the rest of the files
COPY ./src . 

# Make the hardcoded mountpoints for the outputs and inputs
RUN mkdir ./target
RUN mkdir ./input

# Setup the entrypoint
ENTRYPOINT [ "./entry" ]
