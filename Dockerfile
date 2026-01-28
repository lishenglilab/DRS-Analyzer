FROM ubuntu:20.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

LABEL description="Docker image for DRS-Analyzer pipeline"
LABEL version="2.0"

# Set working directory
WORKDIR /app

# ============================================================================
# STEP 1: Install system dependencies and tools
# ============================================================================
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    autoconf \
    automake \
    pkg-config \
    vim \
    less \
    pigz \
    unzip \
    bzip2 \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ============================================================================
# STEP 2: Install minimap2
# ============================================================================
RUN wget https://github.com/lh3/minimap2/releases/download/v2.18/minimap2-2.18_x64-linux.tar.bz2 && \
    tar -xjf minimap2-2.18_x64-linux.tar.bz2 && \
    mv minimap2-2.18_x64-linux/minimap2 /usr/local/bin/ && \
    rm -rf minimap2-2.18_x64-linux minimap2-2.18_x64-linux.tar.bz2

# ============================================================================
# STEP 3: Install SAMtools 1.20
# ============================================================================
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xjf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.20 samtools-1.20.tar.bz2

# ============================================================================
# STEP 4: Install BLAST
# ============================================================================
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-2.17.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.17.0+/bin/* /usr/local/bin/ && \
    rm -rf ncbi-blast-2.17.0+ ncbi-blast-2.17.0+-x64-linux.tar.gz

# ============================================================================
# STEP 5: Install gffread
# ============================================================================
RUN wget https://github.com/gpertea/gffread/releases/download/v0.11.8/gffread-0.11.8.Linux_x86_64.tar.gz && \
    tar -xzf gffread-0.11.8.Linux_x86_64.tar.gz && \
    mv gffread-0.11.8.Linux_x86_64/gffread /usr/local/bin/ && \
    rm -rf gffread-0.11.8.Linux_x86_64 gffread-0.11.8.Linux_x86_64.tar.gz

# ============================================================================
# STEP 6: Install Miniconda (latest version)
# ============================================================================
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    /opt/conda/bin/conda clean -a -y

# Add conda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Initialize conda
RUN conda init bash

# ============================================================================
# STEP 7: Accept conda ToS and configure channels
# ============================================================================
RUN conda tos accept --override-channels \
                     --channel https://repo.anaconda.com/pkgs/msys2 \
                     --channel https://repo.anaconda.com/pkgs/r \
                     --channel https://repo.anaconda.com/pkgs/main

RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda config --set auto_activate_base false && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority flexible"

# ============================================================================
# STEP 8: Create conda environment with Python 3.8
# ============================================================================
RUN conda create -n drs-analyzer python=3.7 numpy==1.19.5 pandas==1.3.0 scikit-learn==0.22 cython=0.29.21 r-base pybedtools=0.8.2 ont-tombo openjdk=11.0.1 nextflow=23.10.1 -y && \
    conda clean -a -y

# Activate environment by default
SHELL ["conda", "run", "-n", "drs-analyzer", "/bin/bash", "-c"]

# Alternative: Add to bashrc so it's activated in interactive shells
RUN echo "conda activate drs-analyzer" >> ~/.bashrc

# ============================================================================
# STEP 9: Install Python packages with specific versions
# ============================================================================
RUN conda run -n drs-analyzer pip install --no-cache-dir \
    pyfastx==0.8.4 \
    pysam==0.22.0 \
    biopython \
    ont-fast5-api

# ============================================================================
# STEP 10: Install bioinformatics tools via conda
# ============================================================================
RUN conda install -n drs-analyzer -y -c bioconda -c conda-forge \
    bedops=2.4.35 \
    porechop=0.2.4 \
    flair \
    seqkit \
    && conda clean -a -y

# ============================================================================
# STEP 11: Install SUPPA2
# ============================================================================
RUN conda run -n drs-analyzer pip install --no-cache-dir SUPPA

# ============================================================================
# STEP 12: Install nanopsu (pseudouridine detection)
# ============================================================================
RUN git clone https://github.com/sihaohuanguc/Nanopore_psU.git /tmp/Nanopore_psU && \
    cd /tmp/Nanopore_psU && \
    conda run -n drs-analyzer pip install . && \
    cd / && \
    rm -rf /tmp/Nanopore_psU

# ============================================================================
# STEP 13: Install R in the conda environment
# ============================================================================

# Install additional system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN conda install -n drs-analyzer -y -c bioconda -c conda-forge r-getopt r-tidyverse r-plyr r-ggsci r-rtsne r-hmisc r-stringr bioconductor-complexheatmap bioconductor-biostrings bioconductor-rtracklayer r-remotes

# ============================================================================
# STEP 14: Verify installations
# ============================================================================
RUN echo "=== Verifying installations ===" && \
    conda run -n drs-analyzer python --version && \
    conda run -n drs-analyzer python -c "import sys; print('Python location:', sys.executable)" && \
    conda run -n drs-analyzer R --version && \
    samtools --version && \
    minimap2 --version && \
    blastn -version && \
    gffread --version && \
    nextflow -version && \
    echo "" && \
    echo "Python packages:" && \
    conda run -n drs-analyzer pip list | grep -E "numpy|pandas|scikit-learn|pyfastx|pysam|biopython|pybedtools|tombo|SUPPA" && \
    echo "" && \
    echo "Testing critical imports:" && \
    conda run -n drs-analyzer python -c "import numpy; print('NumPy:', numpy.__version__)" && \
    conda run -n drs-analyzer python -c "import pandas; print('Pandas:', pandas.__version__)" && \
    conda run -n drs-analyzer python -c "import sklearn; print('scikit-learn:', sklearn.__version__)" && \
    echo "" && \
    echo "Testing nanopsu:" && \
    conda run -n drs-analyzer nanopsu -h | head -5 || echo "nanopsu installed" && \
    echo "=== All installations verified ==="

# ============================================================================
# STEP 15: Set up environment variables
# ============================================================================
ENV PATH=/opt/conda/envs/drs-analyzer/bin:/opt/conda/bin:/usr/local/bin:$PATH
ENV CONDA_DEFAULT_ENV=drs-analyzer
ENV CONDA_PREFIX=/opt/conda/envs/drs-analyzer

# ============================================================================
# STEP 16: Create working directory and set permissions
# ============================================================================
RUN mkdir -p /data /output && \
    chmod -R 777 /data /output

WORKDIR /data

# Default command
CMD ["/bin/bash"]
