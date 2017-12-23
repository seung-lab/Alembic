FROM centos:7

#install MKL julia prereq's
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
        curl \
    && mkdir tmp_install \
    && pushd tmp_install \
    && wget -q https://cmake.org/files/v3.7/cmake-3.7.1.tar.gz \
    && tar zxf cmake-3.7.1.tar.gz \
    && cd cmake-3.7.1 \
    && cmake . \
    && make \
    && make install \
    && popd \
    && rm -rf tmp_install

RUN rpm --import http://yum.repos.intel.com/2017/setup/RPM-GPG-KEY-intel-psxe-runtime-2017 \
    && rpm -Uhv http://yum.repos.intel.com/2017/setup/intel-psxe-runtime-2017-reposetup-1-0.noarch.rpm \
    && yum install -y intel-ifort-runtime-64 \
                      intel-mkl-runtime-64   \
    && cp -r /opt/intel/psxe_runtime_2017.5.239/linux/compiler/lib/intel64_lin/* /lib64 \
    && cp -r /opt/intel/psxe_runtime/linux/mkl/lib/intel64/* /lib64 \
    && rm -rf /opt/intel

# build Julia
RUN git clone https://github.com/JuliaLang/julia.git \
    && cd julia \
    && git checkout v0.6.1 \
    && echo LLVM_VER=4.0.0 > Make.user \
    && echo USE_INTEL_MKL=1 >> Make.user \
    && echo USE_INTEL_LIBM=1 >> Make.user \
    && make -j8 && make install \
    && yes | cp -rf /julia/usr / \ 
    && rm -rf /julia

#install cloudvolume
RUN wget https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm \
    && yum install -y epel-release-latest-7.noarch.rpm \
    && yum install -y python36u \
    && yum install -y python36u-pip \
    && yum install -y python36u-devel \
    && python3.6 -m venv venv \
    && . venv/bin/activate \
    && pip install pip --upgrade \
    && pip install cloud-volume==0.6.8 \
    && mkdir -p /root/.cloudvolume/ \
    && echo $GOOGLE_STORAGE_PROJECT > /root/.cloudvolume/project_name \
    && pip uninstall -y requests \
    && pip install requests

# Install Alembic
RUN julia -e 'Pkg.clone("https://github.com/seung-lab/Alembic.git")' \
 && julia /root/.julia/v0.6/Alembic/UNREGISTERED_REQUIRE.jl \
 && julia -e 'rm(Pkg.dir("PyCall","deps","PYTHON")); Pkg.build("PyCall")' \
 && julia -e 'using Alembic' \
 && ln -s /root/.julia/v0.6/Alembic/src/tasks
