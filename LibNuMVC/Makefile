FILES=numvc.hpp indexed_heap.hpp main.cpp

# OpenMP support is experimental and might be slower (enable using -fopenmp)
numvc: $(FILES)
	g++ main.cpp -g -Wall -std=c++14 -o numvc

numvc-parallel: $(FILES)
	g++ main_parallel.cpp -fopenmp -fopenmp-simd -O3 -Ofast -Wall -std=c++14 -pthread -o numvc

prof: $(FILES)
	g++ main.cpp -g -O2 -pg -no-pie -Wall -std=c++14 -o numvc


debug: $(FILES)
	g++ main.cpp -g -Wall -std=c++14 -o numvc

clean:
	rm numvc

all: numvc
