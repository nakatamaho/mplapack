TAGS="centos8 centos7 ubuntu2004 ubuntu2004intel ubuntu2004mingw64 ubuntu1804 ubuntu2004ppc64le centos8aarch64 centos7aarch64 ubuntu2004aarch64"
for _tag in $TAGS; do 
    rm -f a
    docker run --cidfile="a" mplapack:$_tag echo
    container_id=`cat a`
    echo $container_id
    docker export $container_id | gzip > mplapack:${_tag}.tar.gz ; cat mplapack:${_tag}.tar.gz | gunzip | docker import - mplapack:${_tag}
    rm -f a
done

