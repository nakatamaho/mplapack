# Sample Makefile
# by Nakata Maho
CXX=g++
CXXFLAGS= -std=c++17 -I/home/docker/mplapack/external/i/MPFR/include -I/home/docker/mplapack/external/i/MPC/include -I/home/docker/mplapack/external/i/GMP/include  -I/home/docker/mplapack/mpfrc++/
OPTFLAGS= -O2
LIBFLAGS= -L/home/docker/mplapack/external/i/MPFR/lib -L/home/docker/mplapack/external/i/MPC/lib -L/home/docker/mplapack/external/i/GMP/lib
MPLIBS=-lmpfr -lmpc -lgmp -lgmpxx -lgmp

programs=example example2 example3

all: $(programs)

$(programs): %: %.o
	$(CXX) -o $@ $< $(LIBFLAGS) $(MPLIBS) 

.cpp.o: ../mpreal.h ../mpcomplex.h
	$(CXX) -c $(OPTFLAGS) $(CXXFLAGS) $<

clean:
	rm -rf *.o *~ $(programs) *bak
