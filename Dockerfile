FROM cmalabscience/gattaca-base:0.1.0

# Make the hardcoded mountpoints for the outputs and inputs
RUN mkdir /GATTACA && \
  mkdir /GATTACA/target && mkdir /GATTACA/input && mkdir /GATTACA/logs

WORKDIR /GATTACA

# Copy the source code
COPY ./src/ /GATTACA/

# Setup the entrypoint
ENTRYPOINT [ "/GATTACA/entrypoint.R" ]
