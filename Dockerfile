# Use an official Python runtime as a parent image
# Working from Python base instead of Julia base, because we had conflicts
# with cloudvolume and some of the already installed python packages.
FROM python:2.7

# Attach author
MAINTAINER macrintr

# Download git & other packages
RUN apt-get update
RUN apt-get install -y curl
RUN apt-get install -y libblas-dev liblapack-dev liblapacke-dev gfortran git
RUN apt-get install -y vim

# Install cloud-volume
RUN pip install pip --upgrade
RUN pip install cloud-volume==0.5.5
RUN mkdir -p /root/.cloudvolume/
RUN echo $GOOGLE_STORAGE_PROJECT > /root/.cloudvolume/project_name

# Install Julia
ENV JULIA_PATH /usr/local/julia
RUN curl -fL -o julia.tar.gz https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.1-linux-x86_64.tar.gz
RUN mkdir "$JULIA_PATH"
RUN tar -xzf julia.tar.gz -C "$JULIA_PATH" --strip-components 1;
RUN rm julia.tar.gz
ENV PATH $JULIA_PATH/bin:$PATH

# Install Alembic
RUN julia -e 'Pkg.clone("https://github.com/seung-lab/Alembic.git")'
RUN julia /root/.julia/v0.6/Alembic/UNREGISTERED_REQUIRE.jl
RUN julia -e 'using Alembic'

RUN ln -s /root/.julia/v0.6/Alembic/src/clients /
# Create secrets
# COPY aws-secret.json /root/.cloudvolume/secrets/aws-secret.json
# COPY google-secret.json /root/.cloudvolume/secrets/google-secret.json
