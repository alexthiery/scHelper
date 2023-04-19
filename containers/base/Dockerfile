FROM  rocker/tidyverse:4.0.4

LABEL authors="alex.thiery@crick.ac.uk" \
        description="Docker image containing requirements for schelpr"

# Install basic linux libraries and extras required to build R packages
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        apt-utils \
        build-essential \
        libxtst6 \
        vim

# Install linux libraries listed in Satija Seurat docker image
RUN apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libpng-dev \
        libboost-all-dev \
        libxml2-dev \
        openjdk-8-jdk \
        python3-dev \
        python3-pip \
        wget \
        git \
        libfftw3-dev \
        libgsl-dev \
        libgeos-dev \
        libglpk-dev

# Set CRAN and BioC repositories for package installations
RUN   echo 'local({library(BiocManager); r <- repositories(); r["CRAN"] <- "https://cran.rstudio.com/"; options(repos = r)}); \
            options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' >> /usr/local/lib/R/etc/Rprofile.site

# RUN   R -e "install.packages('future')" && \
#       R -e "install.packages('getopt')" && \
#       R -e "install.packages('optparse')" && \
#       R -e "remotes::install_github('mojaveazure/seurat-disk')" && \
#       R -e "remotes::install_github('carmonalab/STACAS@1.1.0')"
      
# RUN   R -e "devtools::install_github('alexthiery/scHelper@65c756cb5b4bd14bf6a3959aecd0fc106f5183f4', dependencies = TRUE)"



