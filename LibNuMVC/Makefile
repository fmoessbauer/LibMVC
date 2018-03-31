FILES=numvc.hpp indexed_heap.hpp main.cpp

# OpenMP support is experimental and might be slower (enable using -fopenmp)
numvc: $(FILES)
	g++ main.cpp -O3 -Ofast -Wall -std=c++11 -o numvc

prof: $(FILES)
	g++ main.cpp -g -O2 -pg -no-pie -Wall -o numvc


debug: $(FILES)
	g++ main.cpp -g -Wall -o numvc

clean:
	rm numvc

all: numvc
