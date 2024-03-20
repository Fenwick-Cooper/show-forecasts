FROM python:3.11-slim

ARG WORK_HOME=/opt/app
ARG APP_DIR=/code
ARG APP_STORE_DIR=${APP_DIR}/store

# install system libraries
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends ca-certificates curl wget pkg-config libgdal-dev \
    libgeos-dev libproj-dev gdal-bin libcgal-dev libxml2-dev libsqlite3-dev gcc g++ texlive  \
    openssl dvipng texlive-latex-extra texlive-fonts-recommended cm-super libfreetype-dev \
    libfontconfig-dev libjpeg-dev libspng-dev libx11-dev libgbm-dev

# service runtime filesystem directories
RUN mkdir -p ${WORK_HOME} ${APP_DIR} ${APP_STORE_DIR}

COPY poetry.lock pyproject.toml README.md ${APP_DIR}/
COPY ./webapps ${APP_DIR}/webapps
COPY ./cgan_ui ${APP_DIR}/cgan_ui

ARG APP_USER_ID=1000
ARG APP_GROUP_ID

ENV PATH=${WORK_HOME}/bin:/usr/local/bin:/usr/local/chrome-linux64:/usr/local/chromedriver-linux64:$PATH \
    POETRY_HOME=${WORK_HOME} USER_ID=${APP_USER_ID} GROUP_ID=${APP_GROUP_ID} APP_HOME=${APP_DIR} \
    APP_STORE_DIR=${APP_STORE_DIR}

WORKDIR ${APP_HOME}

RUN curl -sSL https://install.python-poetry.org | python3 - && \
    poetry install --all-extras


