# Licensed under a 3-clause BSD style license - see LICENSE.rst
FROM python:slim

# Install build dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential git

# Install test dependencies
RUN pip install hypothesis flake8

# Copy source to image and install development version
COPY . .
RUN pip install -e .[test]
