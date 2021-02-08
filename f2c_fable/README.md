# This directory contains docker image for fable; Automatic Fortran to C++ converter

https://github.com/cctbx/cctbx_project/tree/master/fable
FABLE: Automatic Fortran to C++ conversion (https://doi.org/10.1186/1751-0473-7-5).

## How to build
```docker image build -t maho/f2c_fable:latest .```

## How to run
```
docker run -it maho/f2c_fable:latest
```

* How to convert a sample fortran.
```
$ cd ; fable.cout sample.f
```
