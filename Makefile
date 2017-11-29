FILES=numvc.hpp indexed_heap.hpp main.cpp

numvc: $(FILES)
	g++ main.cpp -O3 -Ofast -Wall -fopenmp -o numvc

debug: $(FILES)
	g++ main.cpp -g -Wall -o numvc

all: numvc
