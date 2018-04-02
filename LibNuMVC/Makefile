RDFILES=../vendor/SFMT/SFMT.c
FILES=numvc.hpp indexed_heap.hpp main.cpp $(RDFILES)
RDFLAGS=-DUSE_SFMT -DSFMT_MEXP=607 -DHAVE_SSE2=1 -msse2 -fno-strict-aliasing
FLAGS=-lpthread -std=c++14 -Wall $(RDFLAGS)

# OpenMP support is experimental and might be slower (enable using -fopenmp)
numvc: $(FILES)
	g++ main.cpp $(RDFILES) -O3 -Ofast -Wall $(FLAGS) -o numvc

prof: $(FILES)
	g++ main.cpp $(RDFILES) -g -O2 -pg -no-pie $(FLAGS) -o numvc


debug: $(FILES)
	g++ main.cpp $(RDFILES) -g $(FLAGS) -o numvc

clean:
	rm numvc

all: numvc
