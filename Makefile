fastvc: fastvc.hpp main.cpp
	g++ main.cpp -O3 -Ofast -Wall -o fastvc

all: fastvc
