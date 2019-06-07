FC=gfortran
FCFLAGS=-O2
CPP=g++
CPPFLAGS=-O2
OBJS =dlas2.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o dlasq6.o test_dlasq1.o
MPOBJS=Rlas2.o Rlasq1.o Rlasq2.o Rlasq3.o Rlasq4.o Rlasq5.o Rlasq6.o test_Rlasq1.o
MPINCLUDEFLAGS=-I/home/maho/mplapack/include -I/home/maho/mplapack/mpfrc++/
MPLDFLAGS=-L/home/maho/mplapack/lib/ -Wl,-rpath /home/maho/mplapack/lib/

.SUFFIXES: .f .cxx

all: test_dlasq1 test_Rlasq1

test_dlasq1: $(OBJS)
	$(CPP) -o $@ $(FCFLAGS) $(OBJS) -llapack -lblas

test_Rlasq1: $(MPOBJS)
	$(CPP) -fopenmp -o $@ $(CPPFLAGS) $(MPOBJS) $(MPLDFLAGS) -lmlapack_mpfr -lmblas_mpfr -lmpfrcxx -lmpc -lmpfr -lgmp

.f.o:
	$(FC) -o $@ $(FCFLAGS) -c $<

.cpp.o:
	$(CPP) -D___MPACK_BUILD_WITH_MPFR___ $(MPINCLUDEFLAGS) -o $@ $(CPPFLAGS) -c $<

clean:
	-rm $(OBJS) $(MPOBJS) *~ $(all)





