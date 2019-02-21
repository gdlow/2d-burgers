default: compile

main.o: main.cpp
	g++ -std=c++11 -Wall -O2 -o main.o -c main.cpp

Burgers.o: Burgers.cpp Burgers.h
	g++ -std=c++11 -Wall -O2 -o Burgers.o -c Burgers.cpp

Model.o: Model.cpp Model.h
	g++ -std=c++11 -Wall -O2 -o Model.o -c Model.cpp

Helpers.o: Helpers.cpp Helpers.h
	g++ -std=c++11 -Wall -O2 -o Helpers.o -c Helpers.cpp	

compile: main.o Burgers.o Model.o Helpers.o
	g++ -o compile main.o Burgers.o Model.o Helpers.o -lblas

#invalid argument exception should be thrown
invalidArg: compile
	./compile 0 0 0 1 

diff: compile
	./compile 0 0 0 1 10 10 1

advx: compile
	./compile 1 0 0 0 10 10 1

advy: compile
	./compile 0 1 0 0 10 10 1

burg: compile
	./compile 1.0 0.5 1.0 0.02 10 10 1

all: compile

.PHONY: clean
clean:
	rm -f *.o compile