# Basic secured-centos image
FROM debian:stretch

LABEL maintainer="Abdoulaye DIALLO abdoulaye.diallo@seq.one"
LABEL maintainer="Dimitri Larue dimitri.larue@seq.one"
LABEL maintainer="Jerome AUDOUX jerome.audoux@seq.one"

# R
ENV R_BASE_VERSION 3.3.3

## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
## Also install stringr to make dococt install (from source) easier
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    littler \
    r-cran-littler \
    r-cran-stringr \
    r-base=${R_BASE_VERSION}-* \
    r-base-dev=${R_BASE_VERSION}-* \
    r-recommended=${R_BASE_VERSION}-* \
      && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> /etc/R/Rprofile.site \
      && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
  && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
  && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
  && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
  && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
  && install.r docopt \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
&& rm -rf /var/lib/apt/lists/*

# Build application
WORKDIR /dekupl
COPY install_r_packages.R .
# Install dependencies
RUN apt-get update && apt-get install -y curl libcurl4-openssl-dev libv8-3.14-dev libssl-dev libxml2-dev

RUN Rscript install_r_packages.R

EXPOSE 8080

COPY . .

ENTRYPOINT ["src/dekupl-viewer.sh"]
