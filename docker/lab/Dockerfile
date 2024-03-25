FROM python:3.10-slim

ARG LAB_USER_ID=10001
ARG LAB_GROUP_ID=10001
ARG WORK_HOME=/opt/mycgan
ARG VENV_DIR=/opt/venv
ARG LAB_STORE_DIR=${WORK_HOME}/store

# install system libraries
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends git ca-certificates pkg-config \
    libgdal-dev libgeos-dev libproj-dev gdal-bin libcgal-dev libxml2-dev libsqlite3-dev  \
    gcc g++ texlive openssl dvipng texlive-latex-extra texlive-fonts-recommended cm-super \
    libfreetype-dev libfontconfig-dev libjpeg-dev libspng-dev libx11-dev libgbm-dev 

# service runtime filesystem directories
RUN mkdir -p ${WORK_HOME} ${LAB_STORE_DIR} ${VENV_DIR}

RUN groupadd --gid ${LAB_GROUP_ID} mycgan && \ 
    useradd --home-dir ${WORK_HOME} --uid ${LAB_USER_ID} --gid ${LAB_GROUP_ID} mycgan && \
    chown -Rf mycgan:mycgan ${WORK_HOME} ${LAB_STORE_DIR} ${VENV_DIR}

USER mycgan
WORKDIR ${WORK_HOME}

ENV PATH=${VENV_DIR}/bin:$PATH WORK_HOME=${WORK_HOME} \
    LAB_STORE_DIR=${LAB_STORE_DIR} VENV_DIR=${VENV_DIR}

RUN git clone --depth 1 https://github.com/jaysnm/ensemble-cgan.git -b Jurre_brishti /tmp/cgan-code && \
    python -m venv ${VENV_DIR} && . ${VENV_DIR}/bin/activate && \
    cd /tmp/cgan-code && pip install --no-cache -e . && pip install jupyterlab

CMD ["jupyter-lab", "--IdentityProvider.token=''", "--no-browser"]
