CXX=g++
FILES=main.cpp
ALGO=NUMVC
FLAGS=-std=c++14 -lpthread -Wall -Werror -DALGO_$(ALGO)

solver: $(FILES)
	$(CXX) main.cpp -O3 -Ofast $(FLAGS) -o solver-$(ALGO)

prof: $(FILES)
	$(CXX) main.cpp -g -O2 -pg -no-pie $(FLAGS) -o solver-$(ALGO)

debug: $(FILES)
	$(CXX) main.cpp -g $(FLAGS) -Wpedantic -o solver-$(ALGO)

clean:
	rm solver-$(ALGO)

all: solver
