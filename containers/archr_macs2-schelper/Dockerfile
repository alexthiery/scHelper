FROM  alexthiery/archr_macs2:3.0

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with R4.1.2 and ArchR and macs2 and schelpr"

# Install scHelper
RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.3.1', dependencies = TRUE)"