FROM  alexthiery/10x-npb-seurat:base-2.8

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with 10x seurat packages, schelpr and shiny"

# Install scHelper
RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.3.4', dependencies = TRUE)"

# Install shiny
RUN   R -e "install.packages('shiny')" && \
      R -e "install.packages('shinydashboard')" && \
      R -e "install.packages('bs4Dash')" && \
      R -e "install.packages('shinythemes')"