# Sample Makefile
# by Nakata Maho
CXX=/usr/local/bin/g++

#MPACKDIR=/work/research/mpack/
#MPACKEXTDIR=/work/research/mpack/external/

MPACKDIR=/home/maho/work/mpack/trunk_daily_work/
MPACKEXTDIR=/home/maho/work/mpack/trunk_daily_work/external/i/
CXXFLAGS=-I$(MPACKDIR)/include -I$(MPACKEXTDIR)/GMP/include -I$(MPACKEXTDIR)/MPFR/include -I$(MPACKEXTDIR)/MPC/include -I$(MPACKDIR)/mpfrc++ -I. -fopenmp
LIBFLAGS=-L$(MPACKDIR)/lib -L$(MPACKEXTDIR)/GMP/lib -L$(MPACKEXTDIR)/MPFR/lib -L$(MPACKEXTDIR)/MPC/lib -fopenmp
OPTFLAGS= -O2
GMPLIBS=-lmblas_gmp -lgmpxx -lgmp
MP=$(MPACKDIR)/mpfrc++/mpreal.cpp $(MPACKDIR)/mpfrc++/mpcomplex.cpp
MPOBJ=mpreal.o mpcomplex.o
BLAS=-L/usr/local/lib -lblas

programs=test_complex___float128 test_complex_mpfr #test_complex_mpfr1 #test_complex_gmp 

all: $(programs)

clean: 
	rm $(programs) *.o
.cpp.o:
	$(CXX) -c $(OPTFLAGS) $(CXXFLAGS) $<

test_complex_gmp: test_complex_gmp.o $(MPOBJ)
	$(CXX) -static -o test_complex_gmp test_complex_gmp.o $(LIBFLAGS) -lgmpxx -lgmp

test_complex_mpfr: test_complex_mpfr.o
	$(CXX) -o test_complex_mpfr test_complex_mpfr.o $(CXXFLAGS) $(LIBFLAGS) $(MP) -lmpc -lmpfr -lgmpxx -lgmp

test_complex___float128: test_complex___float128.o
	$(CXX) -o test_complex___float128 test_complex___float128.o $(CXXFLAGS) $(LIBFLAGS) $(MP) -lmpc -lmpfr -lgmpxx -lgmp -lquadmath


