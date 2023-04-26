FROM  alexthiery/archr:3.0

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with R4.1.2 and ArchR and schelpr"

# Install scHelper 
RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.2.8', dependencies = TRUE)"