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

RUN conda install --override-channels -c conda-forge -c r --yes --name base \
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
    r-renv

# Cannot activate the conda environment easily
# So instead adjust shell to run everything inside Conda
# from here on out
# https://pythonspeed.com/articles/activate-conda-dockerfile/
SHELL ["conda", "run", "--no-capture-output", "/bin/bash", "-c"]

RUN Rscript -e "install.packages('weights', repos='https://cran.r-project.org/', lib='/opt/conda/lib/R/library')"
RUN Rscript -e "install.packages('AMISforInfectiousDiseases', repos='https://cran.r-project.org/', lib='/opt/conda/lib/R/library')"
RUN conda clean -a -y

# Verify installations
RUN R --version && python --version

# https://medium.com/datamindedbe/how-to-access-private-data-and-git-repositories-with-ssh-while-building-a-docker-image-a-ea283c0b4272
RUN mkdir -p -m 0600 ~/.ssh && \
    ssh-keyscan github.com >> ~/.ssh/known_hosts

# Copy the scripts into the container
# Order of the scripts follows the stages of the pipeline
ADD fitting-prep ${FITTING_PREP_DIR}
ADD fitting ${FITTING_DIR}
ADD projections-prep ${PROJECTIONS_PREP_DIR}
ADD projections ${PROJECTIONS_DIR}
ADD post_AMIS_analysis ${TRACHOMA_AMIS_DIR}/post_AMIS_analysis
ADD run_pipeline.py ${TRACHOMA_AMIS_DIR}

ADD https://storage.googleapis.com/ntd-data-storage/pipeline/trachoma/Maps.tar.gz ${FITTING_PREP_DIR}/inputs/
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/trachoma/ESPEN_IU_2021.tar.gz ${FITTING_PREP_DIR}/inputs/
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/trachoma/fitting-prep-artefacts.tar.gz ${FITTING_PREP_DIR}/artefacts/

RUN tar -xzf ${FITTING_PREP_DIR}/inputs/Maps.tar.gz -C ${FITTING_PREP_DIR}/inputs/ && \
    rm ${FITTING_PREP_DIR}/inputs/Maps.tar.gz
RUN tar -xzf ${FITTING_PREP_DIR}/artefacts/fitting-prep-artefacts.tar.gz -C ${FITTING_PREP_DIR}/artefacts/ && \
    rm ${FITTING_PREP_DIR}/artefacts/fitting-prep-artefacts.tar.gz
RUN tar -xzf ${FITTING_PREP_DIR}/inputs/ESPEN_IU_2021.tar.gz -C ${FITTING_PREP_DIR}/inputs/ && \
    rm ${FITTING_PREP_DIR}/inputs/ESPEN_IU_2021.tar.gz

# Get trachoma model
ADD --keep-git-dir git@github.com:NTD-Modelling-Consortium/ntd-model-trachoma.git#v2.0.0 ${TRACHOMA_MODEL_DIR}
RUN cd ${TRACHOMA_MODEL_DIR}

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

VOLUME [${PROJECTIONS_DIR}/artefacts/projections/trachoma]
ENTRYPOINT [ "python", "run_pipeline.py" ]
