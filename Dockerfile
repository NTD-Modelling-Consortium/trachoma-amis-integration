# syntax=docker/dockerfile:1

# https://hub.docker.com/r/condaforge/miniforge3
FROM condaforge/miniforge3

SHELL [ "/bin/bash", "-c" ]
ARG DEBIAN_FRONTEND=noninteractive

ARG TRACHOMA_AMIS_DIR=/ntdmc/trachoma-amis-integration
ARG TRACHOMA_MODEL_DIR=${TRACHOMA_AMIS_DIR}/model/ntd-model-trachoma

ENV TRACHOMA_AMIS_DIR=${TRACHOMA_AMIS_DIR}
ENV TRACHOMA_MODEL_DIR=${TRACHOMA_MODEL_DIR}

RUN apt update && apt install -y \
    build-essential \
    cmake \
    vim \
    curl \
    git \
    unzip \
    openssh-client \
    libssl-dev \
    libgdal-dev \
    libudunits2-dev

RUN conda install --yes --name base \
    python=3.10 \
    pandas \
    joblib \
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
    r-optparse \
    r-sf \
    r-renv && \
    conda clean -a -y

# Cannot activate the conda environment easily
# So instead adjust shell to run everything inside Conda
# from here on out
# https://pythonspeed.com/articles/activate-conda-dockerfile/
SHELL ["conda", "run", "--no-capture-output", "/bin/bash", "-c"]

RUN Rscript -e "install.packages('AMISforInfectiousDiseases', repos='https://cloud.r-project.org/')"

# Verify installations
RUN R --version && python --version

# https://medium.com/datamindedbe/how-to-access-private-data-and-git-repositories-with-ssh-while-building-a-docker-image-a-ea283c0b4272
RUN mkdir -p -m 0600 ~/.ssh && \
    ssh-keyscan github.com >> ~/.ssh/known_hosts

ADD ESPEN_IU_2021 ${TRACHOMA_AMIS_DIR}/ESPEN_IU_2021
ADD Maps/prepare_histories_and_maps.R \
    Maps/prepare_histories_projections.R \
    Maps/trachoma_IU_Match_PB.csv \
    Maps/trachomaComb_IU.csv ${TRACHOMA_AMIS_DIR}/Maps/
ADD post_AMIS_analysis ${TRACHOMA_AMIS_DIR}/post_AMIS_analysis
ADD trachoma_amis ${TRACHOMA_AMIS_DIR}/trachoma_amis
ADD preprocess_for_projections.R \
    realocate_InputPars_MTP.R \
    trachoma_fitting.R \
    RunProjectionsTo2026.py \
    run_fit_and_process.sh \
    run_proj.sh \
    run_historic_simulations.sh ${TRACHOMA_AMIS_DIR}

# Get trachoma model
ADD --keep-git-dir git@github.com:NTD-Modelling-Consortium/ntd-model-trachoma.git ${TRACHOMA_MODEL_DIR}
RUN cd ${TRACHOMA_MODEL_DIR} && git checkout acf7d8b

WORKDIR ${TRACHOMA_AMIS_DIR}

# Install the trachoma model
RUN --mount=type=cache,target=/root/.cache/pip cd ${TRACHOMA_MODEL_DIR} && pip install .

ENV RETICULATE_PYTHON=/opt/conda/bin/python
ENV RETICULATE_PYTHON_FALLBACK=FALSE
ENV TARGET_ESS=500
ENV ESS_THRESHOLD=200
ENV AMIS_N_SAMPLES=1000

VOLUME [${TRACHOMA_MODEL_DIR}/projections/trachoma]
ENTRYPOINT [ "bash", "run_historic_simulations.sh" ]
