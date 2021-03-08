# Sample Makefile
# by Nakata Maho
CXX=/usr/local/bin/g++

#MPLAPACKDIR=/work/research/mplapack/
#MPLAPACKEXTDIR=/work/research/mplapack/external/

MPLAPACKDIR=/home/maho/work/mplapack/trunk_daily_work/
MPLAPACKEXTDIR=/home/maho/work/mplapack/trunk_daily_work/external/i/
CXXFLAGS=-I$(MPLAPACKDIR)/include -I$(MPLAPACKEXTDIR)/GMP/include -I$(MPLAPACKEXTDIR)/MPFR/include -I$(MPLAPACKEXTDIR)/MPC/include -I$(MPLAPACKDIR)/mpfrc++ -I. -fopenmp
LIBFLAGS=-L$(MPLAPACKDIR)/lib -L$(MPLAPACKEXTDIR)/GMP/lib -L$(MPLAPACKEXTDIR)/MPFR/lib -L$(MPLAPACKEXTDIR)/MPC/lib -fopenmp
OPTFLAGS= -O2
GMPLIBS=-lmpblas_gmp -lgmpxx -lgmp
MP=$(MPLAPACKDIR)/mpfrc++/mpreal.cpp $(MPLAPACKDIR)/mpfrc++/mpcomplex.cpp
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


