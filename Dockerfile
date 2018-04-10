FROM continuumio/miniconda

RUN conda env create -f install/rnaseq1.5.yml

# Install R
RUN curl -fsSL https://cran.r-project.org/src/base/R-3/R-3.4.2.tar.gz -o /opt/R-3.4.2.tar.gz && \
    tar xvzf /opt/R-3.4.2.tar.gz -C /opt/ && \
    cd /opt/R-3.4.2 && \
    ./configure && \
    make && \
    make install && \
    rm /opt/R-3.4.2.tar.gz

# Install R Packages v2
RUN echo 'source("https://bioconductor.org/biocLite.R")' > /opt/packages.r && \
    echo 'biocLite()' >> /opt/packages.r && \
    echo 'biocLite(c("Rsubread", "dupRadar", "limma", "lattice", "locfit", "edgeR", "chron", "data.table", "gtools", "gdata", "bitops", "caTools", "gplots", "markdown"))' >> /opt/packages.r && \
    Rscript /opt/packages.r && \
    mkdir -p  /usr/local/lib/R/site-library

# Install iRods
RUN apt-get update \
    apt-get install -y wget \
    wget ftp://ftp.renci.org/pub/irods/releases/4.1.10/ubuntu14/irods-icommands-4.1.10-ubuntu14-x86_64.deb \
    dpkg -i irods-icommands-4.1.10-ubuntu14-x86_64.deb \
    wget ftp://ftp.renci.org/pub/irods/plugins/irods_auth_plugin_krb/1.4/irods-auth-plugin-krb-1.4-ubuntu14-x86_64.deb \
    dpkg -i irods-auth-plugin-krb-1.4-ubuntu14-x86_64.deb \
    apt-get install -y krb5-auth-dialog krb5-locales krb5-user \
    DEFAULT_USER="ubuntu" \
    echo search internal.sanger.ac.uk >> /etc/resolvconf/resolv.conf.d/base \
    resolvconf -u

ADD .irods .
