docker build  -f Dockerfile_CentOS7                -t mplapack:centos7           . 2>&1 | tee log.centos7
docker build  -f Dockerfile_CentOS8                -t mplapack:centos8           . 2>&1 | tee log.centos8
docker build  -f Dockerfile_ubuntu18.04            -t mplapack:ubuntu1804        . 2>&1 | tee log.ubuntu1804
docker build  -f Dockerfile_ubuntu20.04            -t mplapack:ubuntu2004        . 2>&1 | tee log.ubuntu2004
docker build  -f Dockerfile_ubuntu20.04_intel      -t mplapack:ubuntu2004intel   . 2>&1 | tee log.ubuntu2004intel
docker build  -f Dockerfile_ubuntu20.04_mingw64    -t mplapack:ubuntu2004mingw64 . 2>&1 | tee log.ubuntu2004mingw64
#docker buildx -f Dockerfile_ubuntu20.04            -t mplapack:ubuntu2004ppc64le . 2>&1 | tee log.ubuntu2004ppc64le
