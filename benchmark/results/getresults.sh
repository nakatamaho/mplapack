cd /home/docker/mplapack/benchmark/results
mkdir -p 2022/x86_64-apple-darwin20.6.0-Corei5-8500B
mkdir -p 2022/x86_64-apple-darwin20.6.0-Corei5-8500B/logs
cd 2022/x86_64-apple-darwin20.6.0-Corei5-8500B
rsync -arv maho@172.27.109.87:~/MPLAPACK/lib/*/mplapack/benchmark/*.pdf .
rsync -arv maho@172.27.109.87:~/MPLAPACK/lib/*/mplapack/benchmark/*.plt logs/
rsync -arv maho@172.27.109.87:~/MPLAPACK/lib/*/mplapack/benchmark/log.* logs/
 
cd /home/docker/mplapack/benchmark/results
mkdir -p 2022/x86_64-pc-linux-gnu-Xeon-E5-2623-TeslaV100
mkdir -p 2022/x86_64-pc-linux-gnu-Xeon-E5-2623-TeslaV100/logs
cd 2022/x86_64-pc-linux-gnu-Xeon-E5-2623-TeslaV100
rsync -arv maho@dadeba:~/MPLAPACK/lib/*/mplapack/benchmark/*.pdf .
rsync -arv maho@dadeba:~/MPLAPACK/lib/*/mplapack/benchmark/*.plt logs/
rsync -arv maho@dadeba:~/MPLAPACK/lib/*/mplapack/benchmark/log.* logs/

cd /home/docker/mplapack/benchmark/results
mkdir -p 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-RTX3080
mkdir -p 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-RTX3080/logs
cd 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-RTX3080
rsync -arv ~/MPLAPACK/lib/*/mplapack/benchmark/*.pdf .
rsync -arv ~/MPLAPACK/lib/*/mplapack/benchmark/*.plt logs/
rsync -arv ~/MPLAPACK/lib/*/mplapack/benchmark/log.* logs/

cd /home/docker/mplapack/benchmark/results
mkdir -p 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-inteloneapi
mkdir -p 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-inteloneapi/logs
cd 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-inteloneapi
rsync -arv ~/MPLAPACK_INTELONEAPI/lib/*/mplapack/benchmark/*.pdf .
rsync -arv ~/MPLAPACK_INTELONEAPI/lib/*/mplapack/benchmark/*.plt logs/
rsync -arv ~/MPLAPACK_INTELONEAPI/lib/*/mplapack/benchmark/log.* logs/

cd /home/docker/mplapack/benchmark/results
mkdir -p 2022/aarch64-unknown-linux-gnu-raspberry-pi4
mkdir -p 2022/aarch64-unknown-linux-gnu-raspberry-pi4/logs
cd 2022/aarch64-unknown-linux-gnu-raspberry-pi4
rsync -arv maho@172.27.109.65:/data/mnt/MPLAPACK/lib/aarch64-unknown-linux-gnu/mplapack/benchmark/*.pdf .
rsync -arv maho@172.27.109.65:/data/mnt/MPLAPACK/lib/aarch64-unknown-linux-gnu/mplapack/benchmark/*.plt logs/
rsync -arv maho@172.27.109.65:/data/mnt/MPLAPACK/lib/aarch64-unknown-linux-gnu/mplapack/benchmark/log.* logs/
