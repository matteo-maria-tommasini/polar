.SUFFIXES:
.SUFFIXES:	.cpp	.o

ifeq ($(shell uname),Linux)
   CC=g++
   OPT=-Wall -std=c++11 -O3 -march=native -I/usr/include/eigen3/
   LIB=
endif
ifeq ($(shell uname),Darwin)
#   CC=clang++-mp-11
   CC=clang++
   OPT=-Wall -std=c++11 -O3 -march=native -I/opt/local/include/ -I/opt/local/include/eigen3/
   LIB=
endif

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)

.cpp.o:
	$(CC) -c $(OPT) $<

polar.x: $(OBJ)
	$(CC) $(OPT) -o $@ $^ $(LIB)

clean:
	rm -rf $(OBJ) polar.ps polar.pdf 

distclean: clean
	rm -rf polar.x

remake:
	make distclean; make -j 20
