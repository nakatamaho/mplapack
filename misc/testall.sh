#!/bin/sh
#ssh -N -L2003:tesla:22 maho@sdpa.indsys.chuo-u.ac.jp & 
#ssh -N -L2002:fedora13:22 maho@sdpa.indsys.chuo-u.ac.jp &
sh misc/testlinux.sh
sh misc/testmingw.sh
sh misc/testlinux_sdpa.sh
sh misc/testlinux_tesla.sh
sh misc/testlinux_fedora13.sh
sh misc/testfbsd.sh
sh misc/testmaci.sh
sh misc/testmacp.sh
