# 2026 measurements: AMD Ryzen Threadripper 3970X + NVIDIA A100 80GB PCIe

Host: AMD Ryzen Threadripper 3970X (32 cores, 64 threads), NVIDIA A100 80GB
PCIe (sm_80), CUDA 13.1, GCC 15.2, Ubuntu 26.04.  The host was shared with
other jobs; each directory's README states the load and how it was measured.
Same machine as ../../2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100.

| Directory | Contents |
|---|---|
| `rgemm/` | October 2026: Rgemm/Cgemm GPU against CPU (dd/td/qd/mpfr) and the MPFR 512-bit size sweep up to n = 9997, through the public Rgemm (transfers included) |
