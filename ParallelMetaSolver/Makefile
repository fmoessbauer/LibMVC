CXX=g++
FILES=main.cpp
FLAGS=-std=c++14 -lpthread -Wall

solver: $(FILES)
	$(CXX) main.cpp -O3 -Ofast $(FLAGS) -o solver

prof: $(FILES)
	$(CXX) main.cpp -g -O2 -pg -no-pie $(FLAGS) -o solver

debug: $(FILES)
	$(CXX) main.cpp -g $(FLAGS) -Wpedantic -o solver

clean:
	rm solver

all: solver
