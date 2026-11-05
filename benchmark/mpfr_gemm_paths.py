#!/usr/bin/env python3
# Path histogram of the additions in the reference GEMM NN accumulation
# C(i,j) = beta*C(i,j); for l: C(i,j) += (alpha*B(l,j)) * A(i,l), alpha=-1, beta=1,
# inputs uniform in [-1,1) (as benchmark/Rkernel_bench.cpp), exponents/signs only.
import math, random, sys, collections
random.seed(1)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 1024
samples = int(sys.argv[2]) if len(sys.argv) > 2 else 400
for PB in (256, 512, 1024):
    h = collections.Counter(); tot = 0
    for s in range(samples):
        c = random.uniform(-1, 1)
        for l in range(n):
            p = -random.uniform(-1, 1) * random.uniform(-1, 1)
            if c == 0 or p == 0:
                h['zero'] += 1; tot += 1; c += p; continue
            ec, ep = math.frexp(c)[1], math.frexp(p)[1]
            g = abs(ec - ep)
            same = (c > 0) == (p > 0)
            cat = '0' if g == 0 else ('1_63' if g < 64 else ('64_pb' if g < PB else 'big'))
            h[('add_same_g' if same else 'add_opp_g') + cat.replace('0', '0') ] += 1
            r = c + p
            if not same and r != 0 and math.frexp(r)[1] < max(ec, ep) - 63:
                h['cancel64'] += 1
            tot += 1
            c = r
    print(PB, {k: round(v / tot, 4) for k, v in sorted(h.items())})
