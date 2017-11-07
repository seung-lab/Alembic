FROM centos:7

RUN rpm --import http://yum.repos.intel.com/2017/setup/RPM-GPG-KEY-intel-psxe-runtime-2017 
RUN rpm -Uhv http://yum.repos.intel.com/2017/setup/intel-psxe-runtime-2017-reposetup-1-0.noarch.rpm 
RUN yum install -y \
           bzip2 \
           cmake \
           epel-release \
           gcc \
           gcc-c++ \
           gcc-gfortran \
           git \
           intel-icc-runtime-64 \
           intel-ifort-runtime-64 \
           intel-mkl-runtime-64 \
           make \
           m4 \
           patch \
           perl \
           pkgconfig \
           python \
           which \
           wget 
RUN mkdir /tmp_install 
RUN cd /tmp_install \
    && wget -q https://cmake.org/files/v3.7/cmake-3.7.1.tar.gz \
    && tar zxf cmake-3.7.1.tar.gz \
    && cd cmake-3.7.1 \
    && cmake . \
    && make -j 4 \
    && make install \
    && cd /tmp_install \
    && wget -q http://registrationcenter-download.intel.com/akdlm/irc_nas/11113/xppsl-1.4.3-rhel7.2.tar \
    && tar xf xppsl-1.4.3-rhel7.2.tar \
    && yum install -y xppsl-1.4.3/rhel7.2/x86_64/*-memkind-* \
    && cd / \
    && rm -rf /tmp_install \
    && git clone https://github.com/JuliaLang/julia.git \
    && cd /julia \
    && echo USE_ICC = 1 >> Make.user \
    && echo USE_IFC = 1 >> Make.user \
    && echo USE_INTEL_MKL = 1 >> Make.user \
    && echo USE_INTEL_MKL_FFT = 1 >> Make.user \
    && echo USE_INTEL_LIBM = 1 >> Make.user \
    && echo override LLVM_VER=svn >> Make.user \
    && echo override USE_INTEL_MKL=1 >> Make.user \
    && make -C deps get-llvm

RUN yum install -y curl
RUN yum install -y libblas-dev liblapack-dev liblapacke-dev gfortran git  
RUN yum install -y vim 
RUN yum install -y python-devel 
RUN yum install -y python-pip 

RUN pip install pip --upgrad
RUN pip install cloud-volume==0.5.5
RUN mkdir -p /root/.cloudvolume/
RUN echo $GOOGLE_STORAGE_PROJECT > /root/.cloudvolume/project_name

# Install Julia
#ENV JULIA_PATH /usr/local/julia
#RUN curl -fL -o julia.tar.gz https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.1-linux-x86_64.tar.gz
#RUN mkdir "$JULIA_PATH"
#RUN tar -xzf julia.tar.gz -C "$JULIA_PATH" --strip-components 1;
#RUN rm julia.tar.gz
#ENV PATH $JULIA_PATH/bin:$PATH

# Fix requests to make Alembic work
RUN pip uninstall -y requests
RUN pip install requests

# Install Alembic
RUN julia -e 'Pkg.clone("https://github.com/seung-lab/Alembic.git")'
RUN julia /root/.julia/v0.6/Alembic/UNREGISTERED_REQUIRE.jl
RUN julia -e 'using Alembic'

RUN ln -s /root/.julia/v0.6/Alembic/src/clients /

