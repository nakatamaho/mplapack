FROM ubuntu:20.04
RUN apt update
RUN apt -y upgrade
RUN apt install -y sudo 
RUN apt install -y tzdata
# set your timezone
ENV TZ Asia/Tokyo
RUN echo "${TZ}" > /etc/timezone \
  && rm /etc/localtime \
  && ln -s /usr/share/zoneinfo/Asia/Tokyo /etc/localtime \
  && dpkg-reconfigure -f noninteractive tzdata

RUN apt install -y git
RUN apt install -y gcc-9 g++-9 gfortran-9
RUN apt install -y autotools-dev automake libtool make
RUN apt install -y ng-common ng-cjk

ARG DOCKER_UID=1000
ARG DOCKER_USER=docker
ARG DOCKER_PASSWORD=docker
　RUN useradd -m --uid ${DOCKER_UID} --groups sudo ${DOCKER_USER} \
  && echo ${DOCKER_USER}:${DOCKER_PASSWORD} | chpasswd
USER ${DOCKER_USER}
RUN cd /home/$DOCKER_USER && git clone https://github.com/nakatamaho/mplapack.git
RUN cd /home/$DOCKER_USER/mplapack && bash -x misc/reconfig.ubuntu20.04.sh ; make
