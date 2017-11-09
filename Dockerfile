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

RUN yum install -y vim
RUN rpm --import http://yum.repos.intel.com/2017/setup/RPM-GPG-KEY-intel-psxe-runtime-2017 \
    && rpm -Uhv http://yum.repos.intel.com/2017/setup/intel-psxe-runtime-2017-reposetup-1-0.noarch.rpm
RUN yum install -y intel-ifort-runtime-64 \
                   intel-mkl-runtime-64 

RUN cp -r /opt/intel/psxe_runtime_2017.5.239/linux/compiler/lib/intel64_lin/* /lib64
RUN cp -r /opt/intel/psxe_runtime/linux/mkl/lib/intel64/* /lib64

# build Julia
RUN git clone https://github.com/JuliaLang/julia.git \
    && cd julia \
    && git checkout v0.6.1 \
    && echo LLVM_VER=4.0.0 > Make.user \
    && echo USE_INTEL_MKL=1 >> Make.user \
    && echo USE_INTEL_LIBM=1 >> Make.user
RUN cd julia && make -j8 && make install -j9
ENV PATH=/julia/usr/bin:$PATH

#install cloudvolume
RUN yum install -y vim 
RUN yum install -y python-devel 
RUN yum install -y libblas-dev liblapack-dev liblapacke-dev gfortran git 
RUN yum install -y curl
RUN wget https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm && yum install -y epel-release-latest-7.noarch.rpm
RUN yum install -y python-pip 
RUN pip install pip --upgrad
RUN pip install cloud-volume==0.5.5
RUN mkdir -p /root/.cloudvolume/
RUN echo $GOOGLE_STORAGE_PROJECT > /root/.cloudvolume/project_name

# Fix cloudvolume requests conflict
RUN pip uninstall -y requests
RUN pip install requests

# Install Alembic
RUN julia -e 'Pkg.clone("https://github.com/seung-lab/Alembic.git")'
RUN julia /root/.julia/v0.6/Alembic/UNREGISTERED_REQUIRE.jl
RUN julia -e 'using Alembic'

# Link the tasks forlder to root for ez access
RUN ln -s /root/.julia/v0.6/Alembic/src/tasks /
