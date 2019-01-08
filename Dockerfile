FROM continuumio/miniconda3:4.5.11
LABEL maintainer="Ricardo Lebrón <rlebron@go.ugr.es>" \
      authors="Ricardo Lebrón <rlebron@go.ugr.es>" \
      description="Container image containing all requirements for the methflow pipeline" \
      version='0.0.0'

COPY env /env
COPY jar/GenomeAnalysisTK.jar /root/
COPY bin /usr/local/bin
COPY include /usr/local/include
COPY lib /usr/local/lib

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update \
    && apt-get dist-upgrade -y \
    && apt-get install -y build-essential gfortran apt-utils procps \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean -y \
    && conda update -y --all \
    && conda env update -n root -f /env/requirements.yml \
    && conda env update -n root -f /env/main.yml \
    && conda env update -n root -f /env/two_ref.yml \
    && conda env update -n root -f /env/diff_meth.yml \
    && conda env update -n root -f /env/data_dump.yml \
    && conda env update -n root -f /env/tools.yml \
    && conda install --yes -c conda-forge ncurses=6.1 \
    && conda install --yes -c anaconda gcc_linux-64=7.3.0 gxx_linux-64=7.3.0 gfortran_linux-64=7.3.0 \
    && conda install --yes -c bioconda perl-gd=2.69 perl-gdgraph=1.54 \
    && conda install --yes -c conda-forge xorg-libx11=1.6.6 \
    && conda clean -y --all \
    && /opt/conda/opt/gatk-3.8/gatk3-register.sh /root/GenomeAnalysisTK.jar \
    && /opt/conda/bin/cpanm inc::latest \
    && /opt/conda/bin/cpanm GD \
    && /opt/conda/bin/cpanm GD::Graph
