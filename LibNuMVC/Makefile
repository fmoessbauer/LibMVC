FILES=numvc.hpp indexed_heap.hpp main.cpp
FLAGS=-std=c++14 -Wall

# OpenMP support is experimental and might be slower (enable using -fopenmp)
numvc: $(FILES)
	g++ main.cpp -O3 -Ofast -Wall $(FLAGS) -o numvc

prof: $(FILES)
	g++ main.cpp -g -O2 -pg -no-pie $(FLAGS) -o numvc


debug: $(FILES)
	g++ main.cpp -g $(FLAGS) -o numvc

clean:
	rm numvc

all: numvc
