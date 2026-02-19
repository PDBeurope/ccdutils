FROM python:3.10-slim-bookworm AS builder

RUN apt-get update && \
    apt-get install -y procps libexpat1 libxrender1 libxtst6 libxi6

COPY pyproject.toml poetry.lock README.md /ccdutils/
COPY pdbeccdutils /ccdutils/pdbeccdutils

RUN cd /ccdutils && pip install .
