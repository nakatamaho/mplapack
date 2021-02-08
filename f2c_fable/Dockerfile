FROM ubuntu:18.04

LABEL maintainer="maho.nakata@gmail.com"

RUN apt -y update && apt -y upgrade && apt -y dist-upgrade
RUN apt install -y tzdata
ENV TZ Asia/Tokyo
RUN apt -y install python3 ng-cjk wget sudo git expect build-essential

ARG DOCKER_UID=1000
ARG DOCKER_USER=docker
ARG DOCKER_PASSWORD=docker
RUN useradd -u $DOCKER_UID -m $DOCKER_USER --shell /bin/bash && echo "$DOCKER_USER:$DOCKER_PASSWORD" | chpasswd && echo "$DOCKER_USER ALL=(ALL) ALL" >> /etc/sudoers

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.6 1
ARG WORK=/home/$DOCKER_USER

USER docker     
RUN cd ${WORK} && who am i && wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py --no-check-certificate
RUN cd ${WORK} && who am i && python bootstrap.py hot update --builder=cctbx
RUN cd ${WORK} && who am i && wget https://repo.anaconda.com/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh && chmod +x ./Miniconda2-4.5.11-Linux-x86_64.sh 
SHELL ["/bin/bash", "-c"]
RUN echo -e "#!/usr/bin/expect -f\n\
spawn ./Miniconda2-4.5.11-Linux-x86_64.sh\n\
set timeout -1\n\
expect \">>>\"\n\
send \"\\\n\"\n\
expect \"More\"\n\
send \"   \"\n\
expect \">>>\"\n\
send \"yes\\\n\"\n\
expect \">>>\"\n\
send \"/home/${DOCKER_USER}/miniconda2\\\n\"\n\
expect \">>>\"\n\
send \"no\\\n\"" >> ${WORK}/answerfile_miniconda
RUN cat ${WORK}/answerfile_miniconda
RUN chown $DOCKER_USER ${WORK}/answerfile_miniconda
RUN cd ${WORK} && chmod +x answerfile_miniconda && ./answerfile_miniconda

RUN cd ${WORK} && source miniconda2/etc/profile.d/conda.sh && conda create -y --name fable36 python=3.6 && conda activate fable36 && conda install six future && mkdir build36 && cd ${WORK}/build36 && python ../modules/cctbx_project/libtbx/configure.py fable
RUN source ${WORK}/build36/setpaths.sh && cd ${WORK}/build36 && make && cd ${WORK} && cd ${WORK}/build36 && make ; cd ${WORK}
RUN cd ${WORK} && echo "source /home/docker/miniconda2/etc/profile.d/conda.sh" >> .bashrc
RUN cd ${WORK} && echo "conda activate fable36" >> .bashrc
RUN cd ${WORK} && echo "source /home/docker/build36/setpaths.sh" >> .bashrc

