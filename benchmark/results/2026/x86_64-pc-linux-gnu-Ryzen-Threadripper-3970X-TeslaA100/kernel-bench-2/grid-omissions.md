# Additional omissions from the reduced23-point matrix /57-point vector grids

All other series completed their baseline grid. Dense matrix grids would have712 points;
the small powers and interior multiples were deliberately reduced before the baseline.
Remaining MPFR omissions below follow reported>10s cutoffs; see coverage.csv and reductions.log.

| Type | Routine | Method | Points | Max n | Omitted baseline sizes |
|---|---|---|---:|---:|---|
| mpfr512 | RCnrm2 | cpu | 56 | 14680064 | 16777216 |
| mpfr512 | Rtrtri | gpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Rtrtri | cpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Rlauum | gpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Rlauum | cpu | 18 | 2044 | 2046;2048;4092;4095;4096 |
| mpfr512 | Rpotri | gpu | 13 | 511 | 512;1022;1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Rpotri | cpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Rgetri | gpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Rgetri | cpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Ctrtri | gpu | 12 | 496 | 511;512;1022;1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Ctrtri | cpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Clauum | gpu | 12 | 496 | 511;512;1022;1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Clauum | cpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Cpotri | gpu | 9 | 248 | 252;256;496;511;512;1022;1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Cpotri | cpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Cgetri | gpu | 12 | 496 | 511;512;1022;1023;1024;2044;2046;2048;4092;4095;4096 |
| mpfr512 | Cgetri | cpu | 15 | 1022 | 1023;1024;2044;2046;2048;4092;4095;4096 |
