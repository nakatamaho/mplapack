FC=gfortran
FCFLAGS=-O2
CPP=g++
CPPFLAGS=-O2

TESTROUTINE=Rlasq1

OBJS =  dlas2.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o dlasq6.o
MPOBJS= Rlas2.o Rlasq1.o Rlasq2.o Rlasq3.o Rlasq4.o Rlasq5.o Rlasq6.o

MPLAPACKTOP=/home/maho/mplapack/
MPTESTLIB=libmplapack_mpfr_test.so
MPINCLUDEFLAGS=-I${MPLAPACKTOP}/include -I${MPLAPACKTOP}/mpfrc++/
MPLDFLAGS=-L${MPLAPACKTOP}/lib/ -Wl,-rpath=${MPLAPACKTOP}/lib/

.SUFFIXES: .f .cxx

all: test_$(TESTROUTINE)

test_$(TESTROUTINE): $(MPTESTLIB) $(OBJS) test_$(TESTROUTINE).cpp ${MPLAPACKTOP}/mpblas/debug/common/mplapack.debug.cpp
	$(CPP) -D___MPLAPACK_BUILD_WITH_MPFR___ $(MPINCLUDEFLAGS) -o mplapack.debug.o $(CPPFLAGS) -c ${MPLAPACKTOP}/mpblas/debug/common/mplapack.debug.cpp
	$(CPP) -std=c++11 -D___MPLAPACK_BUILD_WITH_MPFR___ $(MPINCLUDEFLAGS) -o test_$(TESTROUTINE).o $(CPPFLAGS) -c test_$(TESTROUTINE).cpp
	$(CPP) -fopenmp -o $@ test_$(TESTROUTINE).o $(OBJS) -L. -lmplapack_mpfr_test -Wl,-rpath=. $(MPLDFLAGS) mplapack.debug.o -lmplapack_mpfr -lmpblas_mpfr -lmpfrcxx -lmpc -lmpfr -lgmp -llapack -lblas -lgfortran

.f.o:
	$(FC) -o $@ $(FCFLAGS) -c $<

.cpp.o:
	$(CPP) -fPIC -DPIC -shared -D___MPLAPACK_BUILD_WITH_MPFR___ $(MPINCLUDEFLAGS) -o $@ $(CPPFLAGS) -c $<


$(MPTESTLIB): $(MPOBJS)
	rm -f $(MPTESTLIB).* $(MPTESTLIB)
	$(CPP) -shared -Wl,-soname=$(MPTESTLIB).1 -o $(MPTESTLIB).1.0.0 $(MPOBJS)
	ln -s $(MPTESTLIB).1.0.0 $(MPTESTLIB)
	ln -s $(MPTESTLIB).1.0.0 $(MPTESTLIB).1

clean:
	rm -f *.o *~ $(all) $(MPTESTLIB).* $(MPTESTLIB) test_$(TESTROUTINE)

