From:continuumio/miniconda3:4.5.11
Bootstrap:docker

%labels
    MAINTAINER Ricardo Lebrón <rlebron@go.ugr.es>
    AUTHORS Ricardo Lebrón <rlebron@go.ugr.es>
    DESCRIPTION Container image containing all requirements for the methflow pipeline
    VERSION 0.0.0

%files
    env /env
    jar/GenomeAnalysisTK.jar /root/
    bin /usr/local/bin
    include /usr/local/include
    lib /usr/local/lib

%post
    /usr/bin/apt-get update
    /usr/bin/apt-get dist-upgrade -y
    /usr/bin/apt-get install -y build-essential gfortran apt-utils procps 
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
    /opt/conda/bin/conda install --yes -c anaconda gcc_linux-64=7.3.0 gxx_linux-64=7.3.0 gfortran_linux-64=7.3.0 
    /opt/conda/bin/conda clean -y --all 
    /opt/conda/opt/gatk-3.8/gatk3-register.sh /root/GenomeAnalysisTK.jar
    /opt/conda/bin/cpanm inc::latest
    /opt/conda/bin/cpanm GD
    /opt/conda/bin/cpanm GD::Graph
