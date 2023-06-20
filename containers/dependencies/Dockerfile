FROM  alexthiery/schelper:base-0.3.0

LABEL authors="alex.thiery@crick.ac.uk" \
        description="Docker image containing requirements for schelpr"

RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.3.5', dependencies = TRUE)"
