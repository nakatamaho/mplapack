TAGS="centos8 centos7 ubuntu2004 ubuntu2004intel ubuntu2004mingw64 ubuntu1804 centos8aarch64 centos7aarch64 ubuntu2004aarch64 ubuntu2004ppc64le ubuntu2004riscv64"
for _tag in $TAGS; do
    rm -f a
    docker run --cidfile="a" -u docker -it mplapack:${_tag} bash -x /home/docker/a.sh 2>&1 | tee log.${_tag}
    container_id=`cat a`
    docker commit $container_id mplapack:${_tag}
    rm -f a
done
