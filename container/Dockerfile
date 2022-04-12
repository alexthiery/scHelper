FROM  rocker/tidyverse:4.0.4

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat"

# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      build-essential \
      libhdf5-dev \
      libboost-filesystem1.67-dev \
      libxtst6 \
      libxt-dev \
      libfontconfig1-dev \
      libfribidi-dev

# Add conda to $PATH
ENV PATH /opt/conda/bin:$PATH

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Update conda and install phate
RUN   conda update -n base -c defaults conda

# Set default python for reticulate
RUN   echo 'RETICULATE_PYTHON = /opt/conda/bin/python' >> /usr/local/lib/R/etc/Renviron

# Set CRAN and BioC repositories for package installations
RUN   echo 'local({library(BiocManager); r <- repositories(); r["CRAN"] <- "https://cran.rstudio.com/"; options(repos = r)}); \
            options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' >> /usr/local/lib/R/etc/Rprofile.site

RUN   R -e "install.packages('future')" && \
      R -e "install.packages('getopt')" && \
      R -e "install.packages('optparse')" && \
      R -e "remotes::install_github('mojaveazure/seurat-disk')" && \
      R -e "remotes::install_github('carmonalab/STACAS@1.1.0')"
      
RUN   R -e "devtools::install_github('alexthiery/scHelper@65c756cb5b4bd14bf6a3959aecd0fc106f5183f4', dependencies = TRUE)"


