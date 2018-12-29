From:continuumio/miniconda3:4.5.11
Bootstrap:docker

%labels
    MAINTAINER Ricardo Lebrón <rlebron@go.ugr.es>
    AUTHORS Ricardo Lebrón <rlebron@go.ugr.es>
    DESCRIPTION Container image containing all requirements for the methflow pipeline
    VERSION 0.0.0

%files
    env /env
    bin/GenomeAnalysisTK.jar /root/
    bin/M-IndelRealigner /usr/local/bin/
    bin/software_versions /usr/local/bin/

%post
    /usr/bin/apt-get update
    /usr/bin/apt-get install -y apt-utils
    /usr/bin/apt-get install -y procps
    /bin/rm -rf /var/lib/apt/lists/*
    /usr/bin/apt-get clean -y
    /opt/conda/bin/conda update -y --all 
    /opt/conda/bin/conda env update -n root -f /env/requirements.yml 
    /opt/conda/bin/conda env update -n root -f /env/main.yml 
    /opt/conda/bin/conda env update -n root -f /env/two_ref.yml 
    /opt/conda/bin/conda env update -n root -f /env/diff_meth.yml 
    /opt/conda/bin/conda env update -n root -f /env/data_dump.yml 
    /opt/conda/bin/conda env update -n root -f /env/tools.yml 
    /opt/conda/bin/conda install --yes -c conda-forge ncurses=6.1 
    /opt/conda/bin/conda clean -y --all 
    /opt/conda/opt/gatk-3.8/gatk3-register.sh /root/GenomeAnalysisTK.jar
