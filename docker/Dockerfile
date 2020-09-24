FROM debian:10

ENV DEBIAN_FRONTEND=noninteractive

# Install system libraries required for python and R installations
RUN apt-get update && apt-get install -y --no-install-recommends build-essential apt-utils ca-certificates zlib1g-dev gfortran locales libxml2-dev libcurl4-openssl-dev libssl-dev libreadline6-dev xorg-dev libcairo-dev libpango1.0-dev libbz2-dev liblzma-dev libffi-dev libsqlite3-dev

# Install common linux tools
RUN apt-get update && apt-get install -y --no-install-recommends wget curl htop less nano vim emacs git

# Configure default locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8



# Download and compile R from source
WORKDIR /app/R36
RUN wget https://cran.rstudio.com/src/base/R-3/R-3.6.2.tar.gz
RUN tar xvfz R-3.6.2.tar.gz && rm R-3.6.2.tar.gz

WORKDIR /app/R36/R-3.6.2
RUN ./configure --enable-R-shlib --with-cairo --with-libpng --prefix=/app/R36/
RUN make && make install
ENV PATH="/app/R36/bin:${PATH}"
ENV LD_LIBRARY_PATH="/app/R36/lib/R/lib:${LD_LIBRARY_PATH}"

RUN Rscript -e "update.packages(ask=FALSE, repos='https://ftp.gwdg.de/pub/misc/cran/')"
RUN Rscript -e "install.packages(c('devtools', 'gam', 'RColorBrewer', 'BiocManager'), repos='https://ftp.gwdg.de/pub/misc/cran/')"
RUN Rscript -e "devtools::install_github(c('patzaw/neo2R', 'patzaw/BED'))"
RUN Rscript -e "BiocManager::install(c('scran', 'MAST','monocle','ComplexHeatmap','limma','slingshot','clusterExperiment'), version = '3.10')"

RUN Rscript -e 'writeLines(capture.output(sessionInfo()), "../../package_versions_r.txt")' --default-packages=scran,RColorBrewer,slingshot,monocle,gam,clusterExperiment,ggplot2,plyr,MAST



# Download and compile python from source
WORKDIR /app/python37
RUN wget https://www.python.org/ftp/python/3.7.5/Python-3.7.5.tgz
RUN tar zxfv Python-3.7.5.tgz && rm Python-3.7.5.tgz

WORKDIR /app/python37/Python-3.7.5
RUN ./configure --enable-optimizations --with-lto --prefix=/app/python37/
RUN make && make install

WORKDIR /app/python37
RUN rm -rf /opt/python37/Python-3.7.5
RUN ln -s /app/python37/bin/python3 /app/python37/bin/python
RUN ln -s /app/python37/bin/pip3 /app/python37/bin/pip
ENV PATH="/app/python37/bin:${PATH}"

COPY python-packages.txt /app/python37/python-packages.txt
RUN pip install --no-cache-dir -U pip wheel setuptools cmake
RUN pip install --no-cache-dir -r /app/python37/python-packages.txt
RUN jupyter contrib nbextension install --system
RUN jupyter nbextension enable --py widgetsnbextension

RUN pip freeze > ../../package_versions_py.txt

COPY .bashrc_docker /root/.bashrc
COPY .profile_docker /root/.profile

RUN apt-get clean -y && apt-get autoremove -y
