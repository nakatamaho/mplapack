TAGS="fable centos8 centos7 ubuntu2004 ubuntu2004intel ubuntu2004mingw64 ubuntu1804 centos8aarch64 centos7aarch64 ubuntu2004aarch64 ubuntu2004ppc64le ubuntu2004riscv64 ubuntu2004s390x debianmips64le"

GITCOMMAND="git checkout v1.0.0"
HEAD="TRUE"

for _tag in $TAGS; do
if [ $HEAD != "TRUE" ]; then
    rm -f a
    docker run --cidfile="a" -u docker -it mplapack:${_tag} bash -c -x "rm -rf /home/docker/.ccache /home/docker/ccache.tar.gz ; cp /home/docker/aorg.sh /home/docker/a.sh ; sed -i 's/%%GIT%%/$GITCOMMAND/g' /home/docker/a.sh ; cat /home/docker/a.sh"
    container_id=`cat a`
    docker commit $container_id mplapack:${_tag}
else
    rm -f a
    docker run --cidfile="a" -u docker -it mplapack:${_tag} bash -c -x "cp /home/docker/aorg.sh /home/docker/a.sh ; sed -i 's/%%GIT%%/echo/g' /home/docker/a.sh ; cat /home/docker/a.sh"
    container_id=`cat a`
    docker commit $container_id mplapack:${_tag}
fi
    rm -f a
    docker run --cidfile="a" -u docker -it mplapack:${_tag} bash -x /home/docker/a.sh 2>&1 | tee log.${_tag}
    container_id=`cat a`
    docker commit $container_id mplapack:${_tag}

    docker run --cidfile="a" mplapack:$_tag echo
    container_id=`cat a`
    echo $container_id
    docker export $container_id | gzip > mplapack:${_tag}.tar.gz ; cat mplapack:${_tag}.tar.gz | gunzip | docker import - mplapack:${_tag}
    rm -f a
    docker images | grep none | awk '{print "docker rmi -f " $3}' | sh
done
