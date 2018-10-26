CXX=g++
FILES=main.cpp
FLAGS=-std=c++14 -Wall -I ../LibMVC

release: $(FILES)
	$(CXX) main.cpp -O3 -Ofast $(FLAGS) -DSOLVER=NuMVC  -o numvc
	$(CXX) main.cpp -O3 -Ofast $(FLAGS) -DSOLVER=FastVC -o fastvc
	$(CXX) main.cpp -O3 -Ofast $(FLAGS) -DSOLVER=NuMVC  -DPARALLEL -lpthread -o parallel_numvc
	$(CXX) main.cpp -O3 -Ofast $(FLAGS) -DSOLVER=FastVC -DPARALLEL -lpthread -o parallel_fastvc

prof: $(FILES)
	$(CXX) main.cpp -g -O2 -pg -no-pie -fno-omit-frame-pointer $(FLAGS) -DSOLVER=NuMVC -o numvc

debug: $(FILES)
	$(CXX) main.cpp -g $(FLAGS) -Wpedantic -DSOLVER=NuMVC -o numvc

clean:
	rm numvc
	rm fastvc
	rm parallel_numvc
	rm parallel_fastvc

all: release
