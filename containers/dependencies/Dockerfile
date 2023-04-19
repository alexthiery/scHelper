FROM  alexthiery/schelper:base-0.2.7

LABEL authors="alex.thiery@crick.ac.uk" \
        description="Docker image containing requirements for schelpr"

RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.2.7', dependencies = TRUE)"
