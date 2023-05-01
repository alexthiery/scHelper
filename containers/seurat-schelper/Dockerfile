FROM  alexthiery/10x-npb-seurat:base-2.8

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with 10x seurat packages and schelpr"

# Install scHelper
RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.3.3', dependencies = TRUE)"