# TR_Bear_Connectivity pipeline — canonical computational environment.
#
# This Dockerfile IS the archived environment for this study: it is the single,
# authoritative recipe from which the analysis container is rebuilt. No prebuilt
# image is distributed; build it from this file.
#
# R packages install from a dated CRAN snapshot and the two source-installed
# tools (ENMTML, UNICOR) are pinned through build arguments, so a rebuild years
# from now resolves the same versions as the original run.
#
# Build:
#   docker build -t trbear:latest -f Dockerfile .
#
# Reproduce the exact study environment by pinning the two upstream tools to the
# commits used (see README §2):
#   docker build -t trbear:latest -f Dockerfile . \
#       --build-arg ENMTML_REF=<commit-sha> --build-arg UNICOR_REF=<commit-sha>
#
# Optional — run on an HPC with Apptainer/Singularity instead of Docker:
#   docker save trbear:latest | gzip > trbear.tar.gz
#   # transfer trbear.tar.gz to the cluster, then:
#   apptainer build trbear.sif docker-archive://trbear.tar.gz
# All cluster scripts then run via `apptainer exec trbear.sif ...`.

FROM rocker/geospatial:4.3.2

## Upstream refs for the two source-installed tools. Default to the latest
## upstream state; override at build time to pin the exact commits used.
## This study used ENMTML 1.0.0.
ARG ENMTML_REF=main
ARG UNICOR_REF=main

LABEL project="TR_Bear_Connectivity" \
      description="Turkey-wide brown bear SDM (ENMTML) + UNICOR connectivity" \
      r_version="4.3.2" \
      python="3.10" \
      sdm_pkg="ENMTML"

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    TZ=Europe/Istanbul \
    R_LIBS_USER=/usr/local/lib/R/site-library \
    UNICOR_HOME=/opt/programs/UNICOR \
    CONEFOR_HOME=/opt/programs/conefor_cli \
    OMP_NUM_THREADS=1 \
    PROJ_NETWORK=OFF \
    MAKEFLAGS="-j2" \
    PATH="/opt/programs/UNICOR:/opt/programs/conefor_cli:${PATH}"

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        build-essential cmake git wget curl unzip ca-certificates \
        software-properties-common \
        python3 python3-pip python3-dev python3-venv \
        libgdal-dev gdal-bin libgeos-dev libproj-dev libudunits2-dev \
        libssl-dev libxml2-dev libcurl4-openssl-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libcairo2-dev libxt-dev libmagick++-dev \
        openjdk-17-jdk-headless \
        libnetcdf-dev libhdf5-dev \
        libv8-dev libnode-dev \
        bc parallel && \
    ln -sf /usr/bin/python3 /usr/local/bin/python && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

## ---- R base + spatial + tidy + plotting ----
RUN R -e "options(Ncpus = parallel::detectCores()); \
      install.packages(c( \
        'remotes','pak','devtools', \
        'terra','sf','stars','raster','sp','geosphere', \
        'tidyverse','data.table','readxl','writexl','openxlsx', \
        'ggplot2','ggpubr','ggspatial','scales','patchwork','cowplot', \
        'tidyterra','rasterVis','gridExtra','viridis','RColorBrewer','colorspace', \
        'rnaturalearth','rnaturalearthdata', \
        'gtools','foreach','doParallel','future','future.apply','furrr', \
        'logger','glue','fs','here','jsonlite','yaml', \
        'fasterize' \
      ), repos = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-12-01')"

RUN R -e "remotes::install_github('ropensci/rnaturalearthhires')" || true

## ---- LEGACY rgdal + rgeos for ENMTML (archived on CRAN, Oct 2023 snapshot) ----
RUN R -e "install.packages('rgdal', repos='https://packagemanager.posit.co/cran/__linux__/jammy/2023-10-01'); \
      if (!requireNamespace('rgdal', quietly=TRUE)) remotes::install_github('cran/rgdal', upgrade='never'); \
      if (!requireNamespace('rgdal', quietly=TRUE)) stop('rgdal install FAILED'); \
      cat('rgdal OK\n')"
RUN R -e "install.packages('rgeos', repos='https://packagemanager.posit.co/cran/__linux__/jammy/2023-10-01'); \
      if (!requireNamespace('rgeos', quietly=TRUE)) remotes::install_github('cran/rgeos', upgrade='never'); \
      if (!requireNamespace('rgeos', quietly=TRUE)) stop('rgeos install FAILED'); \
      cat('rgeos OK\n')"

## ---- ENMTML Imports (full DESCRIPTION list) ----
RUN R -e "install.packages(c( \
        'ade4','caret','flexclust','foreach','glmnet','igraph', \
        'pgirmess','doParallel','dplyr','plyr','tools','usdm', \
        'dismo','gam','kernlab','randomForest','maxnet','adehabitatHS', \
        'maxlike','gbm','spThin','mgcv','sp','raster' \
      ), repos = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-12-01'); \
      ok <- sapply(c('ade4','caret','flexclust','glmnet','igraph','pgirmess', \
                     'dismo','gam','kernlab','randomForest','maxnet', \
                     'adehabitatHS','maxlike','gbm','spThin','usdm'), \
                   requireNamespace, quietly = TRUE); \
      if (!all(ok)) stop('Missing ENMTML Imports: ', paste(names(ok)[!ok], collapse=', '))"

## ---- SDM extras (downstream use) ----
RUN R -e "install.packages(c( \
        'PresenceAbsence','ecospat','blockCV','landscapemetrics', \
        'SDMtune','rJava','earth','nnet','MASS','e1071','reshape2' \
      ), repos = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-12-01')"

## ---- GRaF (GAU algorithm; optional) ----
RUN R -e "remotes::install_github('goldingn/GRaF')" || true

## ---- ENMTML itself — git clone + R CMD INSTALL (verbose) ----
##   REPRODUCIBILITY: this recipe is the canonical environment. ENMTML is
##   installed from source at the ref given by ENMTML_REF (default: upstream
##   main). This study used ENMTML 1.0.0; pass the corresponding commit via
##   --build-arg ENMTML_REF=<sha> for a bit-for-bit rebuild.
RUN cd /tmp && \
    git clone https://github.com/andrefaa/ENMTML.git enmtml_src && \
    git -C enmtml_src checkout "${ENMTML_REF}" && \
    R CMD INSTALL --no-multiarch --with-keep.source enmtml_src 2>&1 | tee /tmp/enmtml_install.log && \
    R --vanilla -e "if (!requireNamespace('ENMTML', quietly=TRUE)) stop('ENMTML namespace not loadable'); \
                     cat('ENMTML', as.character(packageVersion('ENMTML')), 'OK\n')" && \
    rm -rf /tmp/enmtml_src

## ---- Python (UNICOR + general geospatial) ----
RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel && \
    pip3 install --no-cache-dir \
        numpy scipy pandas \
        gdal==$(gdal-config --version) \
        rasterio fiona shapely pyproj \
        geopandas \
        networkx \
        matplotlib seaborn \
        openpyxl xlrd \
        psutil tqdm pyyaml

## ---- UNICOR connectivity simulator ----
##   REPRODUCIBILITY: no '|| true' -- a failed clone aborts the build instead
##   of silently producing an image without UNICOR. Pin the exact commit with
##   --build-arg UNICOR_REF=<sha>.
RUN mkdir -p /opt/programs && cd /opt/programs && \
    git clone https://github.com/ComputationalEcologyLab/UNICOR.git && \
    git -C UNICOR checkout "${UNICOR_REF}" && \
    chmod -R a+rX UNICOR

## ---- Conefor CLI (OPTIONAL) ----
##   Not required by this pipeline: every connectivity index (PC, IIC, dPC,
##   dIIC, ECA) is computed in R following the published Conefor definitions,
##   and no analysis script invokes this binary. It is installed only as a
##   convenience for independent cross-checking. The upstream download URL is
##   unstable, so a failed fetch is tolerated and simply leaves it absent.
RUN cd /opt/programs && \
    wget -q http://www.conefor.org/files/usuarios/Conefor_command_line_Linux64.zip -O conefor.zip || true && \
    if [ -s conefor.zip ]; then \
        unzip -o conefor.zip -d conefor_cli && chmod -R a+rx conefor_cli; \
    fi && rm -f conefor.zip

## ---- fasterize + lwgeom, installed and verified individually ----
## (Both silently failed when bundled into the big install.packages() call
##  above -- install.packages() does not stop on a single package failure in
##  a vector call, so the miss was only caught much later by the final
##  verification step, after the whole image had already been built once.
##  Installing + checking them in their own RUN layer makes the failure loud
##  immediately, and keeps the fix cheap to rebuild: everything above this
##  line is unchanged, so Docker reuses it from cache.)
RUN R -e "install.packages('fasterize', repos = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-12-01'); \
          if (!requireNamespace('fasterize', quietly=TRUE)) stop('fasterize install FAILED'); \
          cat('fasterize', as.character(packageVersion('fasterize')), 'OK\n')"
RUN R -e "install.packages('lwgeom', repos = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-12-01'); \
          if (!requireNamespace('lwgeom', quietly=TRUE)) stop('lwgeom install FAILED'); \
          cat('lwgeom', as.character(packageVersion('lwgeom')), 'OK\n')"

## ---- Final verification (fails the build if anything required is missing) ----
RUN R --no-save -e "stopifnot( \
        requireNamespace('ENMTML'), \
        requireNamespace('terra'), \
        requireNamespace('sf'), \
        requireNamespace('dismo'), \
        requireNamespace('gam'), \
        requireNamespace('kernlab'), \
        requireNamespace('randomForest'), \
        requireNamespace('maxnet'), \
        requireNamespace('gbm'), \
        requireNamespace('spThin'), \
        requireNamespace('landscapemetrics'), \
        requireNamespace('lwgeom'), \
        requireNamespace('fasterize') \
    ); cat('ENMTML version:', as.character(packageVersion('ENMTML')), '\n'); \
       cat('fasterize version:', as.character(packageVersion('fasterize')), '\n')" && \
    python3 -c "import numpy, scipy, pandas, networkx, rasterio, geopandas; print('python ok')" && \
    test -d /opt/programs/UNICOR && echo "UNICOR present"

WORKDIR /workspace
CMD ["/bin/bash"]
