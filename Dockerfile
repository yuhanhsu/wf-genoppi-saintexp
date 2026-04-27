### Dockerfile for Genoppi and SAINTexpress installation

# base image from Rocker: https://github.com/rocker-org/rocker-versioned2/wiki
# use shiny-verse for preinstalled shiny and tidyverse packages
# pick R version compatible with Bioconductor: https://www.bioconductor.org/install/
FROM rocker/shiny-verse:4.4.0
LABEL maintainer="Yu-Han Hsu <yuhanhsu@broadinstitute.org>"

# prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# install build-essential (make and gcc), wget, ca-certificates (for secure download)
RUN apt-get update && apt-get install -y \
	build-essential \
	wget \
	ca-certificates \
	&& rm -rf /var/lib/apt/lists/*

# set build directory
WORKDIR /tmp/build

# download, extract, and compile SAINTexpress source code
RUN wget -O SAINTexpress_v3.6.3__2018-03-09.tar.gz \
	https://sourceforge.net/projects/saint-apms/files/SAINTexpress_v3.6.3__2018-03-09.tar.gz/download \
	&& tar -xzvf SAINTexpress_v3.6.3__2018-03-09.tar.gz \
	&& cd SAINTexpress_v3.6.3__2018-03-09 \
	&& make all

# copy SAINTexpress binaries to /usr/local/bin (included in default system PATH)
RUN cp SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress* /usr/local/bin/

# install genoppi from GitHub development branch
RUN Rscript -e "options(timeout=1000); \
	remotes::install_github('lagelab/Genoppi',ref='development')"

# copy Genoppi and SAINTexpress R scripts to /usr/local/src
WORKDIR /usr/local/src
COPY runGenoppi.r runSAINTexpressInt.r runSAINTexpressSpc.r ./

# clean up build files
RUN rm -rf /tmp/build

