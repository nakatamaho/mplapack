cd /home/docker/mplapack/benchmark/results
mkdir -p 2022/x86_64-apple-darwin20.6.0-Corei5-8500B
cd 2022/x86_64-apple-darwin20.6.0-Corei5-8500B
rsync -arv maho@172.27.109.87:~/MPLAPACK/lib/*/mplapack/benchmark/log.* .
rsync -arv maho@172.27.109.87:~/MPLAPACK/lib/*/mplapack/benchmark/*.plt .
rsync -arv maho@172.27.109.87:~/MPLAPACK/lib/*/mplapack/benchmark/*.pdf .
 
cd /home/docker/mplapack/benchmark/results
mkdir -p 2022/x86_64-pc-linux-gnu-Xeon-E5-2623-TeslaV100
cd 2022/x86_64-pc-linux-gnu-Xeon-E5-2623-TeslaV100
rsync -arv maho@dadeba:~/MPLAPACK/lib/*/mplapack/benchmark/log.* .
rsync -arv maho@dadeba:~/MPLAPACK/lib/*/mplapack/benchmark/*.plt .
rsync -arv maho@dadeba:~/MPLAPACK/lib/*/mplapack/benchmark/*.pdf .

