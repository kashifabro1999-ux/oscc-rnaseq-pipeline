# syntax=docker/dockerfile:1
FROM python:3.11-slim

# ---------------------------------------------------------------------------
# 1. System libraries: R itself, build tools, and every header/lib the
#    CRAN + Bioconductor packages below compile against. This is the same
#    list used in the Ubuntu tutorial (Section 0), adapted to Debian package
#    names for this base image.
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-base-dev \
    build-essential gfortran \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libjpeg-dev \
    libglpk-dev libgmp-dev \
    libcairo2-dev libxt-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    pandoc git cmake \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---------------------------------------------------------------------------
# 2. R packages — installed as two separate layers (CRAN, then Bioconductor)
#    BEFORE the application code is copied in. This means editing
#    oscc_web.py or an R script later and rebuilding will reuse these two
#    layers from Docker's cache instead of reinstalling ~30 R packages
#    every time.
# ---------------------------------------------------------------------------
COPY install_cran_packages.R /app/install_cran_packages.R
RUN Rscript /app/install_cran_packages.R

COPY install_bioc_packages.R /app/install_bioc_packages.R
RUN Rscript /app/install_bioc_packages.R

# ---------------------------------------------------------------------------
# 3. Python side: just Streamlit (the pipeline itself is orchestrated R,
#    not reimplemented in Python).
# ---------------------------------------------------------------------------
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# ---------------------------------------------------------------------------
# 4. Application code (changes most often -> copied last so it invalidates
#    the fewest cached layers on rebuild).
# ---------------------------------------------------------------------------
COPY oscc_pipeline.py /app/oscc_pipeline.py
COPY oscc_web.py /app/oscc_web.py
COPY scripts/ /app/scripts/
# If you want the Space to be able to actually RUN the analysis (not just
# browse pre-computed results), also add your public data files, e.g.:
# COPY data/ /app/data/

# Pre-create the working folders the pipeline expects.
RUN mkdir -p /app/data /app/results /app/plots /app/logs

# ---------------------------------------------------------------------------
# 5. Cloud Run injects a PORT environment variable at runtime (default 8080)
#    and requires the container to listen on whatever value it provides —
#    unlike Hugging Face Spaces, which fixes this at 7860. Using the shell
#    form of CMD (not the ["exec","form"] array) lets $PORT expand at
#    container start time.
# ---------------------------------------------------------------------------
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
EXPOSE 8080

CMD streamlit run oscc_web.py --server.port=${PORT:-8080} --server.address=0.0.0.0
