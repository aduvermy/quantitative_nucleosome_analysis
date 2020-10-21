FROM lbmc/r-base:4.0.0

## copy files
COPY install_packages.R /usr/bin/
COPY parsingLog.R /usr/bin/
COPY diagramm_Venn.R /usr/bin/
COPY get_normalization_factor.R /usr/bin/
COPY normalized_occupancy.R /usr/bin/
COPY occupancy_statistics.R /usr/bin/

RUN apk add  --update --no-cache bash file
RUN chmod +x /usr/bin/*.R &&\
Rscript /usr/bin/install_packages.R  ## install R-packages

CMD ["bash"]


