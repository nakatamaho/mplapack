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
mkdir -p 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X
mkdir -p 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X/logs
cd 2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X
rsync -arv ~/MPLAPACK/lib/*/mplapack/benchmark/*.pdf .
rsync -arv ~/MPLAPACK/lib/*/mplapack/benchmark/*.plt logs/
rsync -arv ~/MPLAPACK/lib/*/mplapack/benchmark/log.* logs/

