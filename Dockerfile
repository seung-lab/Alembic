FROM centos:7

# install julia
ENV JULIA_PATH /usr/local/julia
RUN curl -fL -o julia.tar.gz https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.1-linux-x86_64.tar.gz
RUN mkdir "$JULIA_PATH"
RUN tar -xzf julia.tar.gz -C "$JULIA_PATH" --strip-components 1;
RUN rm julia.tar.gz
ENV PATH $JULIA_PATH/bin:$PATH

RUN yum install -y \
        bzip2 \
        cmake \
        gcc \
        gcc-c++ \
        gcc-gfortran \
        git \
        make \
        m4 \
        openssl \
        openssl-devel \
        patch \
        perl \
        pkgconfig \
        python \
        wget \
        which \
        vim \
        python-devel \
        libblas liblapack liblapacke gfortran git \
        curl

# install cloudvolume
RUN    wget https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm  \
    && yum install -y epel-release-latest-7.noarch.rpm  \
    && yum install -y python34-setuptools \
    && yum install -y python34-devel \
    && easy_install-3.4 pip
RUN pip3 install pip --upgrade 
RUN pip3 install cloud-volume
RUN mkdir -p /root/.cloudvolume/ 
RUN echo $GOOGLE_STORAGE_PROJECT > /root/.cloudvolume/project_name 
RUN pip3 uninstall -y requests 
RUN pip3 install requests
RUN pip3 install protobuf==3.1.0

RUN pwd
RUN pwd

# Install Alembic
RUN ls 
ENV PYTHON /usr/bin/python3
RUN julia -e 'Pkg.clone("https://github.com/seung-lab/Alembic.git")' \
    && julia /root/.julia/v0.6/Alembic/UNREGISTERED_REQUIRE.jl \
    && julia -e 'using Alembic' \
    && ln -s /root/.julia/v0.6/Alembic/src/tasks
