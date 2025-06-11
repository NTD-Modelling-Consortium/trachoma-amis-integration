# syntax=docker/dockerfile:1

# https://hub.docker.com/r/condaforge/miniforge3
FROM condaforge/miniforge3

SHELL [ "/bin/bash", "-c" ]
ARG DEBIAN_FRONTEND=noninteractive

ARG TRACHOMA_AMIS_DIR=/ntdmc/trachoma-amis-integration
ARG TRACHOMA_MODEL_DIR=${TRACHOMA_AMIS_DIR}/model/ntd-model-trachoma
ARG FITTING_PREP_DIR=${TRACHOMA_AMIS_DIR}/fitting-prep
ARG FITTING_DIR=${TRACHOMA_AMIS_DIR}/fitting
ARG PROJECTIONS_PREP_DIR=${TRACHOMA_AMIS_DIR}/projections-prep
ARG PROJECTIONS_DIR=${TRACHOMA_AMIS_DIR}/projections

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
    run_pipeline.py ${TRACHOMA_AMIS_DIR}

# Get trachoma model
ADD --keep-git-dir git@github.com:NTD-Modelling-Consortium/ntd-model-trachoma.git ${TRACHOMA_MODEL_DIR}
RUN cd ${TRACHOMA_MODEL_DIR} && git checkout acf7d8b

WORKDIR ${TRACHOMA_AMIS_DIR}

# Install the trachoma model
RUN --mount=type=cache,target=/root/.cache/pip cd ${TRACHOMA_MODEL_DIR} && pip install .

ENV TRACHOMA_AMIS_DIR=${TRACHOMA_AMIS_DIR}
ENV TRACHOMA_MODEL_DIR=${TRACHOMA_MODEL_DIR}
ENV PATH_TO_FITTING_PREP_ARTEFACTS="$TRACHOMA_AMIS_DIR/fitting-prep/artefacts"
ENV PATH_TO_FITTING_PREP_INPUTS="$TRACHOMA_AMIS_DIR/fitting-prep/inputs"
ENV PATH_TO_FITTING_PREP_SCRIPTS="$TRACHOMA_AMIS_DIR/fitting-prep/scripts"

ENV PATH_TO_FITTING_ARTEFACTS="$TRACHOMA_AMIS_DIR/fitting/artefacts"
ENV PATH_TO_FITTING_INPUTS="$TRACHOMA_AMIS_DIR/fitting/inputs"
ENV PATH_TO_FITTING_SCRIPTS="$TRACHOMA_AMIS_DIR/fitting/scripts"

ENV PATH_TO_PROJECTIONS_PREP_ARTEFACTS="$TRACHOMA_AMIS_DIR/projections-prep/artefacts"
ENV PATH_TO_PROJECTIONS_PREP_INPUTS="$TRACHOMA_AMIS_DIR/projections-prep/inputs"
ENV PATH_TO_PROJECTIONS_PREP_SCRIPTS="$TRACHOMA_AMIS_DIR/projections-prep/scripts"

ENV PATH_TO_PROJECTIONS_ARTEFACTS="$TRACHOMA_AMIS_DIR/projections/artefacts"
ENV PATH_TO_PROJECTIONS_INPUTS="$TRACHOMA_AMIS_DIR/projections/inputs"
ENV PATH_TO_PROJECTIONS_SCRIPTS="$TRACHOMA_AMIS_DIR/projections/scripts"

ENV RETICULATE_PYTHON=/opt/conda/bin/python
ENV RETICULATE_PYTHON_FALLBACK=FALSE

VOLUME [${TRACHOMA_MODEL_DIR}/projections/trachoma]
ENTRYPOINT [ "python", "run_pipeline.py" ]
