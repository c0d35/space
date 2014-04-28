
CPP	=	g++ -O3  -std=c++11 -fexceptions -pedantic
INCPATH  =
LIBPATH = 
SYSLIBS = 	-lstdc++ -lpthread
PRGS	=	test

all: le_build test


test: le_build/test.o
	$(CPP)  $^ -I$(INCPATH) -L$(LIBPATH) $(SYSLIBS) -o $@

clean:
	rm -rf le_build test;\

le_build/test.o: test.cpp le_space.h
	$(CPP) -c $< -I$(INCPATH) -L$(LIBPATH) -o $@

le_build:
	mkdir ./le_build












