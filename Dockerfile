# Use official Julia runtime as parent image
FROM julia

# Attach author
MAINTAINER macrintr

RUN apt-get update
RUN apt-get install -y libblas-dev liblapack-dev liblapacke-dev gfortran git python-dev python-pip

# Setup cloud-volume dirs & vars
RUN pip install pip --upgrade
RUN pip install python --upgrade
RUN pip install cloud-volume==0.5.5
RUN mkdir -p /root/.cloudvolume/
RUN echo $GOOGLE_STORAGE_PROJECT > /root/.cloudvolume/project_name

RUN julia -e 'Pkg.clone("https://github.com/seung-lab/Alembic.git")'
RUN julia /root/.julia/v0.6/Alembic/UNREGISTERED_REQUIRE.jl
# WORKDIR /root/.julia/v0.6/Alembic

#
# RUN pip install -r requirements.txt


# Create secrets
# COPY aws-secret.json /root/.cloudvolume/secrets/aws-secret.json
# COPY google-secret.json /root/.cloudvolume/secrets/google-secret.json