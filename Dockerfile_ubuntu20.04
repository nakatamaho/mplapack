FROM ubuntu:20.04

ARG MPLAPACK_VER="2.0.1"

ARG DEBIAN_FRONTEND=noninteractive
RUN apt -y update
RUN apt -y upgrade
RUN apt install -y sudo 
RUN apt install -y tzdata
# set your timezone
ENV TZ Asia/Tokyo
RUN echo "${TZ}" > /etc/timezone \
  && rm /etc/localtime \
  && ln -s /usr/share/zoneinfo/Asia/Tokyo /etc/localtime \
  && dpkg-reconfigure -f noninteractive tzdata

RUN apt install -y build-essential gfortran python3
RUN apt install -y autotools-dev automake libtool gnuplot
RUN apt install -y git wget time parallel
RUN apt install -y ng-common ng-cjk emacs-nox
RUN sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1

ARG DOCKER_UID=1000
ARG DOCKER_USER=docker
ARG DOCKER_PASSWORD=docker
RUN useradd -u $DOCKER_UID -m $DOCKER_USER --shell /bin/bash && echo "$DOCKER_USER:$DOCKER_PASSWORD" | chpasswd && echo "$DOCKER_USER ALL=(ALL) ALL" >> /etc/sudoers
USER ${DOCKER_USER}

RUN cd /home/$DOCKER_USER && echo "cd /home/$DOCKER_USER" >> .bashrc
RUN cd /home/$DOCKER_USER && wget https://github.com/nakatamaho/mplapack/releases/download/v${MPLAPACK_VER}/mplapack-${MPLAPACK_VER}.tar.xz
RUN cd /home/$DOCKER_USER && tar xvfJ mplapack-${MPLAPACK_VER}.tar.xz
RUN cd /home/$DOCKER_USER/mplapack-${MPLAPACK_VER}
ARG CXX="g++"
ARG CC="gcc"
ARG FC="gfortran"
RUN /bin/bash -c 'set -ex && \
    ARCH=`uname -m` && \
    if [ "$ARCH" == "x86_64" ]; then \
       cd /home/$DOCKER_USER/mplapack-${MPLAPACK_VER} && ./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes --enable-benchmark=yes; \
    else \
       cd /home/$DOCKER_USER/mplapack-${MPLAPACK_VER} && ./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes --enable-benchmark=yes ; \
    fi'
RUN cd /home/$DOCKER_USER/mplapack-${MPLAPACK_VER} && make -j`getconf _NPROCESSORS_ONLN`
RUN cd /home/$DOCKER_USER/mplapack-${MPLAPACK_VER} && make install
