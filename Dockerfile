# syntax=docker/dockerfile:1
#
# IMPORTANT: base image is a VERSION-PINNED R build (rocker/r-ver:4.3.3),
# not a generic Debian/Ubuntu image with "whatever r-base apt currently
# ships". This is deliberate and fixes a real problem we hit repeatedly
# testing this pipeline interactively: R itself drifts forward over time
# just like CRAN does, so "apt-get install r-base" today gives a DIFFERENT
# R version (and therefore a different Bioconductor release) than it would
# have a year ago. That mismatch between an old, frozen Bioconductor release
# and constantly-moving-forward CRAN packages caused multiple real build
# failures (treeio, ggtree, fgsea all broke against newer CRAN dependency
# versions). Pinning R to the exact version this pipeline was validated
# against (4.3.3 -> Bioconductor 3.18) makes this reproducible forever,
# rather than a moving target that can silently break again on a future
# rebuild.
FROM rocker/r-ver:4.3.3

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv \
    build-essential gfortran \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libjpeg-dev \
    libglpk-dev libgmp-dev libuv1-dev \
    libcairo2-dev libxt-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    pandoc git cmake \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY install_cran_packages.R /app/install_cran_packages.R
RUN Rscript /app/install_cran_packages.R

COPY install_bioc_packages.R /app/install_bioc_packages.R
RUN Rscript /app/install_bioc_packages.R

COPY requirements.txt /app/requirements.txt
RUN pip3 install --no-cache-dir --break-system-packages -r /app/requirements.txt

COPY oscc_pipeline.py /app/oscc_pipeline.py
COPY oscc_web.py /app/oscc_web.py
COPY scripts/ /app/scripts/
# COPY data/ /app/data/

RUN mkdir -p /app/data /app/results /app/plots /app/logs

ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
EXPOSE 8080

CMD streamlit run oscc_web.py --server.port=${PORT:-8080} --server.address=0.0.0.0
