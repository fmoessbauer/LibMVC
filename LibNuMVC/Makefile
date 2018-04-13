CXX=g++
FILES=numvc.hpp internal/indexed_heap.hpp main.cpp
FLAGS=-std=c++14 -Wall

# OpenMP support is experimental and might be slower (enable using -fopenmp)
numvc: $(FILES)
	$(CXX) main.cpp -O3 -Ofast $(FLAGS) -o numvc

prof: $(FILES)
	$(CXX) main.cpp -g -O2 -pg -no-pie -fno-omit-frame-pointer $(FLAGS) -o numvc


debug: $(FILES)
	$(CXX) main.cpp -g $(FLAGS) -Wpedantic -o numvc

clean:
	rm numvc

all: numvc
