# syntax=docker/dockerfile:1

# https://biohpc.cornell.edu/doc/CondaInContainer.html
FROM continuumio/miniconda3:latest

SHELL [ "/bin/bash", "-c" ]
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update && apt install -y \
        build-essential \
        cmake \
        vim \
        curl \
        git \
        unzip \
        openssh-client

RUN conda config --add channels conda-forge --add channels defaults --add channels r

RUN conda install --yes --name base \
        python=3.10 \
        r-base \
        r-reticulate \
        r-truncnorm \
        r-dplyr \
        r-magrittr \
        r-invgamma \
        r-tidyr \
        r-devtools \
        r-openxlsx \
        r-hmisc \
        r-mclust \
        r-mnormt \
        r-rcpp \
        r-rcpparmadillo

RUN Rscript -e "install.packages('AMISforInfectiousDiseases', repos='https://cloud.r-project.org/', lib='/opt/conda/lib/R/library/')"
RUN conda clean -a -y

RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="/root/.local/bin:$PATH"

# Verify installations
RUN R --version && python --version && poetry --version

ARG TRACHOMA_AMIS_DIR=/ntdmc/trachoma-amis-integration
ARG TRACHOMA_MODEL_DIR=${TRACHOMA_AMIS_DIR}/model/ntd-model-trachoma
ARG FITTING_PREP_DIR=${TRACHOMA_AMIS_DIR}/fitting-prep
ARG FITTING_DIR=${TRACHOMA_AMIS_DIR}/fitting
ARG PROJECTIONS_PREP_DIR=${TRACHOMA_AMIS_DIR}/projections-prep
ARG PROJECTIONS_DIR=${TRACHOMA_AMIS_DIR}/projections

RUN mkdir -p ${FITTING_PREP_DIR}/{inputs,artefacts,scripts} && \
        mkdir -p ${FITTING_PREP_DIR}/artefacts/{Maps,model_output} && \
        mkdir -p ${FITTING_DIR}/{inputs,artefacts,scripts} && \
        mkdir -p ${PROJECTIONS_PREP_DIR}/{inputs,artefacts,scripts} && \
        mkdir -p ${PROJECTIONS_DIR}/{inputs,artefacts,scripts}

# https://medium.com/datamindedbe/how-to-access-private-data-and-git-repositories-with-ssh-while-building-a-docker-image-a-ea283c0b4272
RUN mkdir -p -m 0600 ~/.ssh && \
        ssh-keyscan github.com >> ~/.ssh/known_hosts

# Get trachoma model
ADD git@github.com:NTD-Modelling-Consortium/ntd-model-trachoma.git ${TRACHOMA_AMIS_DIR}/model/ntd-model-trachoma

WORKDIR ${TRACHOMA_AMIS_DIR}

# Install the trachoma model
RUN cd ${TRACHOMA_MODEL_DIR} && \
        python -m venv .venv && \
        source .venv/bin/activate && \
        pip install ${TRACHOMA_MODEL_DIR}



ENTRYPOINT [ "bash" ]
