CXX=g++
FILES=fastvc.hpp main.cpp
FLAGS=-std=c++14 -Wall

# OpenMP support is experimental and might be slower (enable using -fopenmp)
fastvc: $(FILES)
	$(CXX) main.cpp -O3 -Ofast $(FLAGS) -o fastvc

prof: $(FILES)
	$(CXX) main.cpp -g -O2 -pg -no-pie $(FLAGS) -o fastvc


debug: $(FILES)
	$(CXX) main.cpp -g $(FLAGS) -Wpedantic -o fastvc

clean:
	rm fastvc

all: fastvc
