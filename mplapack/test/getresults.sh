##########
cd /home/docker/mplapack/mplapack/test/eig/
cd results/aarch64-unknown-linux-gnu
rsync -arv maho@172.27.109.65:/data/mnt/MPLAPACK/lib/aarch64-unknown-linux-gnu/mplapack/test/eig/log* .

##########
cd /home/docker/mplapack/mplapack/test/eig/
cd results/x86_64-apple-darwin20.6.0
rsync -arv maho@172.27.109.87:~/MPLAPACK/lib/x86_64-apple-darwin20.6.0/mplapack/test/eig/log* .

##########
cd /home/docker/mplapack/mplapack/test/eig/
cd results/x86_64-pc-linux-gnu
rsync -arv ~/MPLAPACK/lib/x86_64-pc-linux-gnu/mplapack/test/eig/log* .

##########
cd /home/docker/mplapack/mplapack/test/eig/
cd results/x86_64-pc-linux-gnu_inteloneapi
rsync -arv ~/MPLAPACK_INTELONEAPI/lib/x86_64-pc-linux-gnu/mplapack/test/eig/log* .

##########
cd /home/docker/mplapack/mplapack/test/eig/
cd results/x86_64-w64-mingw32
rsync -arv ~/MPLAPACK_MINGW/lib/x86_64-w64-mingw32/mplapack/test/eig/log* .

##########
cd /home/docker/mplapack/mplapack/test/lin/
cd results/aarch64-unknown-linux-gnu
rsync -arv maho@172.27.109.65:/data/mnt/MPLAPACK/lib/aarch64-unknown-linux-gnu/mplapack/test/lin/log* .

##########
cd /home/docker/mplapack/mplapack/test/lin/
cd results/x86_64-apple-darwin20.6.0
rsync -arv maho@172.27.109.87:~/MPLAPACK/lib/x86_64-apple-darwin20.6.0/mplapack/test/lin/log* .

##########
cd /home/docker/mplapack/mplapack/test/lin/
cd results/x86_64-pc-linux-gnu
rsync -arv ~/MPLAPACK/lib/x86_64-pc-linux-gnu/mplapack/test/lin/log* .

##########
cd /home/docker/mplapack/mplapack/test/lin/
cd results/x86_64-pc-linux-gnu_inteloneapi
rsync -arv ~/MPLAPACK_INTELONEAPI/lib/x86_64-pc-linux-gnu/mplapack/test/lin/log* .

##########
cd /home/docker/mplapack/mplapack/test/lin/
cd results/x86_64-w64-mingw32
rsync -arv ~/MPLAPACK_MINGW/lib/x86_64-w64-mingw32/mplapack/test/lin/log* .
