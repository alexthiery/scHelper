FROM  alexthiery/archr_macs2_seurat:3.0

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with R4.1.2 and ArchR and macs2 and seurat and schelpr"

# Install scHelper
RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.3.5', dependencies = TRUE)"

# Install shiny
RUN   R -e "install.packages('shiny')" && \
      R -e "install.packages('shinydashboard')" && \
      R -e "install.packages('bs4Dash')" && \
      R -e "install.packages('shinythemes')"