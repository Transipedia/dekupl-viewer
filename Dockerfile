# Basic secured-centos image
FROM r-base:latest

LABEL maintainer="Abdoulaye DIALLO abdoulaye.diallo@seq.one"
LABEL maintainer="Dimitri Larue dimitri.larue@seq.one"
LABEL maintainer="Jerome AUDOUX jerome.audoux@seq.one"

RUN apt-get update && apt-get install -y curl libcurl4-openssl-dev libv8-3.14-dev

# Install dependencies
RUN apt-get update && apt-get install -y curl libcurl4-openssl-dev libv8-3.14-dev

# setup R lib for shiny
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
Rscript -e "install.packages('ggplot2')" && \
Rscript -e "install.packages('shiny')" && \
Rscript -e "install.packages('DT')" && \
Rscript -e "install.packages('shinydashboard')" && \
Rscript -e "install.packages('optparse')" && \
Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('ComplexHeatmap')" && \
Rscript -e "install.packages('factoextra')"
RUN Rscript -e "install.packages('calibrate')"


# Expose 8080 for Docker
EXPOSE 8080

LABEL service_type="dekupl" service_name="viewer"

COPY files/health /usr/local/bin/health
HEALTHCHECK --start-period=4s --interval=10s --timeout=10s --retries=3 CMD /usr/local/bin/health

# Build application
WORKDIR /root/app/
COPY src /root/app/

ENTRYPOINT ["Rscript", "app.R"]
