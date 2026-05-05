### Dockerfile for Genoppi, SAINTexpress, and gcloud CLI installation

# base image from Rocker: https://github.com/rocker-org/rocker-versioned2/wiki
# use shiny-verse for preinstalled shiny and tidyverse packages
# pick R version compatible with Bioconductor: https://www.bioconductor.org/install/
FROM rocker/shiny-verse:4.4.0
LABEL maintainer="Yu-Han Hsu <yuhanhsu@broadinstitute.org>"

# prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# install dependencies and gcloud CLI
RUN apt-get update && apt-get install -y \
	curl \
	ca-certificates \
	gnupg \
	build-essential \
	wget \
	&& echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
	&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
	&& apt-get update && apt-get install -y google-cloud-cli \
	&& rm -rf /root/.cache/pip/ \
	&& find /usr/lib/google-cloud-sdk -name "*.pyc" -delete \
	&& find /usr/lib/google-cloud-sdk -name "*__pycache__*" -delete \
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

