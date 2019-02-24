default: compile

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

Burgers.o: Burgers.cpp Burgers.h
	mpicxx -std=c++11 -Wall -O2 -o Burgers.o -c Burgers.cpp

Burgers2P.o: Burgers2P.cpp Burgers2P.h
	mpicxx -std=c++11 -Wall -O2 -o Burgers2P.o -c Burgers2P.cpp

Model.o: Model.cpp Model.h
	mpicxx -std=c++11 -Wall -O2 -o Model.o -c Model.cpp

Helpers.o: Helpers.cpp Helpers.h
	mpicxx -std=c++11 -Wall -O2 -o Helpers.o -c Helpers.cpp	

compile: main.o Burgers.o Burgers2P.o Model.o Helpers.o
	mpicxx -o compile main.o Burgers.o Burgers2P.o Model.o Helpers.o -lblas

#invalid argument exception should be thrown
invalidArg: compile
	./compile 0 0 0 1 

diff: compile
	./compile 0 0 0 1 10 10 1 0

advx: compile
	./compile 1 0 0 0 10 10 1 0

advy: compile
	./compile 0 1 0 0 10 10 1 0

burg: compile
	./compile 1.0 0.5 1.0 0.02 10 10 1 0

all: compile

## For running MPI processes
diffp: compile
	mpiexec -np 2 ./compile 0 0 0 1 10 10 1 1

advxp: compile
	mpiexec -np 2 ./compile 1 0 0 0 10 10 1 1

advyp: compile
	mpiexec -np 2 ./compile 0 1 0 0 10 10 1 1

burgp: compile
	mpiexec -np 2 ./compile 1.0 0.5 1.0 0.02 10 10 1 1

.PHONY: clean
clean:
	rm -f *.o compile compilep