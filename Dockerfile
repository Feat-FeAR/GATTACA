FROM rocker/r-ver:4.1.0

RUN apt update

RUN apt install -y libgmp3-dev libmpfr-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev proj-bin libproj-dev libcairo2-dev libxt-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
  
COPY ./renv.lock .
COPY ./install_dep .
COPY ./install.R .

RUN ./install_dep

COPY ./GATTACA2.R .
COPY ./STALKER_Functions.R .

