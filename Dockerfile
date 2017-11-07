# Use an official Python runtime as a parent image
# Working from Python base instead of Julia base, because we had conflicts
# with cloudvolume and some of the already installed python packages.
FROM rscohn2/julia-mkl 

# Attach author
MAINTAINER macrintr

# Download git & other packages
RUN yum update -y
RUN yum install -y curl
RUN yum install -y libblas-dev liblapack-dev liblapacke-dev gfortran git
RUN yum install -y vim
RUN yum install -y python-dev
RUN yum install -y python-pip 

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

# Fix requests to make Alembic work
RUN pip uninstall -y requests
RUN pip install requests

# Install Alembic
RUN julia -e 'Pkg.clone("https://github.com/seung-lab/Alembic.git")'
RUN julia /root/.julia/v0.6/Alembic/UNREGISTERED_REQUIRE.jl
RUN julia -e 'using Alembic'

RUN ln -s /root/.julia/v0.6/Alembic/src/clients /
# Create secrets
# COPY aws-secret.json /root/.cloudvolume/secrets/aws-secret.json
# COPY google-secret.json /root/.cloudvolume/secrets/google-secret.json
