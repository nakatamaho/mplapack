FC=gfortran
FCFLAGS=-O2
CPP=g++
CPPFLAGS=-O2

TESTROUTINE=Rlasq1
TESTROUTINED=dlasq1

OBJS =  test_$(TESTROUTINED).o dlas2.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o dlasq6.o
MPOBJS=               Rlas2.o Rlasq1.o Rlasq2.o Rlasq3.o Rlasq4.o Rlasq5.o Rlasq6.o
MPTESTLIB=libmplapack_mpfr_test.so

MPINCLUDEFLAGS=-I/home/maho/mplapack/include -I/home/maho/mplapack/mpfrc++/
MPLDFLAGS=-L/home/maho/mplapack/lib/ -Wl,-rpath=/home/maho/mplapack/lib/

.SUFFIXES: .f .cxx

all: test_$(TESTROUTINED) test_$(TESTROUTINE)

test_dlasq1: $(OBJS)
	$(CPP) -o $@ $(FCFLAGS) $(OBJS) -llapack -lblas

test_Rlasq1: $(MPTESTLIB) test_$(TESTROUTINE).cpp
	$(CPP) -D___MPACK_BUILD_WITH_MPFR___ $(MPINCLUDEFLAGS) -o test_$(TESTROUTINE).o $(CPPFLAGS) -c test_$(TESTROUTINE).cpp
	$(CPP) -fopenmp -o $@ test_$(TESTROUTINE).o -L. -lmplapack_mpfr_test -Wl,-rpath=. $(MPLDFLAGS) -lmlapack_mpfr -lmblas_mpfr -lmpfrcxx -lmpc -lmpfr -lgmp

.f.o:
	$(FC) -o $@ $(FCFLAGS) -c $<

.cpp.o:
	$(CPP) -fPIC -DPIC -shared -D___MPACK_BUILD_WITH_MPFR___ $(MPINCLUDEFLAGS) -o $@ $(CPPFLAGS) -c $<

$(MPTESTLIB): $(MPOBJS)
	rm -f $(MPTESTLIB).* $(MPTESTLIB)
	$(CPP) -shared -Wl,-soname=$(MPTESTLIB).1 -o $(MPTESTLIB).1.0.0 $(MPOBJS)
	ln -s $(MPTESTLIB).1.0.0 $(MPTESTLIB)
	ln -s $(MPTESTLIB).1.0.0 $(MPTESTLIB).1

clean:
	rm -f $(OBJS) $(MPOBJS) *~ $(all) $(MPTESTLIB).* $(MPTESTLIB) test_$(TESTROUTINED) test_$(TESTROUTINE)

