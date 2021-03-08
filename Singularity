Bootstrap: docker
From: ubuntu:20.04

%post
    apt-get update -y && apt-get install -y build-essential bc wget git unzip

    rm -rf /var/lib/apt/lists/*


    CONDA_INSTALL_PATH="/opt/miniconda"
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $CONDA_INSTALL_PATH

    rm -r miniconda.sh

    . $CONDA_INSTALL_PATH/etc/profile.d/conda.sh

    mkdir /bart_dir
    cd /bart_dir
    mkdir run
    git clone --recursive https://github.com/exosports/BART BART/
    conda env create -f /bart_dir/BART/environment.yml
    conda activate bart
    echo ". $CONDA_INSTALL_PATH/etc/profile.d/conda.sh" >> /.singularity_bashrc
    echo "conda activate bart" >> /.singularity_bashrc



    cd /bart_dir/BART/modules/transit/
    make

    cd /bart_dir/BART/modules/MCcubed/
    make


%environment
    CONDA_INSTALL_PATH="/opt/miniconda"
    CONDA_BIN_PATH="/opt/minconda/bin"
    SINGULARITY_SHELL=/bin/bash
    topdir="/bart_dir"
    action="${0##*/}"
    if [ "$action" == "shell" ]; then
        if [ "${SINGULARITY_SHELL:-}" == "/bin/bash" ]; then
            set -- --noprofile --init-file /.singularity_bashrc
        elif test -z "${SINGULARITY_SHELL:-}"; then
            export SINGULARITY_SHELL=/bin/bash
            set -- --noprofile --init-file /.singularity_bashrc
        fi
    fi
    export PATH="$CONDA_BIN_PATH:$PATH"
    export topdir

%runscript
    exec /bin/bash --noprofile --init-file /.singularity_bashrc "$@"
%labels
    Maintainer David Wright <davecwright@knights.ucf.edu>
    License See README or /bart_dir/BART/license in container
    Version v1.0
