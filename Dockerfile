FROM continuumio/miniconda3:4.5.11
LABEL maintainer="Ricardo Lebrón <rlebron@go.ugr.es>" \
      authors="Ricardo Lebrón <rlebron@go.ugr.es>" \
      description="Container image containing all requirements for the methflow pipeline" \
      version='0.0.0'

COPY environment.yml bin/GenomeAnalysisTK.jar /
COPY bin/M-IndelRealigner bin/software_versions /usr/local/bin/

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update \
    && apt-get install -y procps \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean -y \
    && conda update -y --all \
    && conda env update -n root -f /requirements.yml \
    && conda env update -n root -f /environment.yml \
    && conda install --yes -c conda-forge ncurses=6.1 \
    && conda clean -y --all \
    && /opt/conda/opt/gatk-3.8/gatk3-register.sh /GenomeAnalysisTK.jar
