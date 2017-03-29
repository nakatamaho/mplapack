LDIR=/usr/local/lib/

example: example.o mpreal.o
	g++ example.o mpreal.o $(LDIR)libmpfr.a $(LDIR)libgmp.a -o example

mpreal.o: mpreal.cpp mpreal.h
	g++ -c mpreal.cpp

example.o: example.cpp mpreal.h
	g++ -c example.cpp

clean:
	rm -f *.o 

