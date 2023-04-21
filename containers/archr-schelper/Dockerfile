FROM  alexthiery/archr:2.0

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with R4.1.2 and ArchR and schelpr"

# Install scHelper 
RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.2.4', dependencies = TRUE)"