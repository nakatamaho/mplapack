TAGS="mplapack:ubuntu2004 mplapack:centos8 mplapack:centos7 mplapack:ubuntu2004ppc64le mplapack:ubuntu2004intel mplapack:ubuntu2004mingw64 mplapack:ubuntu1804"
for _tag in $TAGS; do 
    rm -f a
    docker run --cidfile="a" $_tag ccache -C
    container_id=`cat a`
    echo $container_id
    docker export $container_id | gzip > ${_tag}.tar.gz ; cat ${_tag}.tar.gz | gunzip | docker import - ${_tag}
done

