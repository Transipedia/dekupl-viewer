# Basic secured-centos image
FROM r-base:3.5.1

LABEL maintainer="Abdoulaye DIALLO abdoulaye.diallo@seq.one"
LABEL maintainer="Dimitri Larue dimitri.larue@seq.one"
LABEL maintainer="Jerome AUDOUX jerome.audoux@seq.one"

# Build application
WORKDIR /dekupl
COPY install_r_packages.R .
# Install dependencies
RUN apt-get update && apt-get install -y curl libcurl4-openssl-dev libv8-3.14-dev libssl-dev libxml2-dev

RUN Rscript install_r_packages.R

EXPOSE 8080

LABEL service_type="dekupl" service_name="viewer"

COPY files/health /usr/local/bin/health
HEALTHCHECK --start-period=4s --interval=10s --timeout=10s --retries=3 CMD /usr/local/bin/health

COPY . .
WORKDIR /dekupl/src

ENTRYPOINT ["Rscript", "app.R"]
