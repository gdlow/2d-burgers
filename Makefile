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

## For running MPI processes (replaces g++ with mpicxx)
## TODO: build variable to replace g++ and mpicxx
mainp.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o mainp.o -c main.cpp

Burgersp.o: Burgers.cpp Burgers.h
	mpicxx -std=c++11 -Wall -O2 -o Burgersp.o -c Burgers.cpp

Modelp.o: Model.cpp Model.h
	mpicxx -std=c++11 -Wall -O2 -o Modelp.o -c Model.cpp

Helpersp.o: Helpers.cpp Helpers.h
	mpicxx -std=c++11 -Wall -O2 -o Helpersp.o -c Helpers.cpp	

compilep: mainp.o Burgersp.o Modelp.o Helpersp.o
	mpicxx -o compilep mainp.o Burgersp.o Model_mpi.o Helpersp.o -lblas

diffp: compilep
	mpiexec -np 2 ./compilep 0 0 0 1 10 10 1

advxp: compilep
	mpiexec -np 2 ./compilep 1 0 0 0 10 10 1

advyp: compilep
	mpiexec -np 2 ./compilep 0 1 0 0 10 10 1

burgp: compilep
	mpiexec -np 2 ./compilep 1.0 0.5 1.0 0.02 10 10 1

.PHONY: clean
clean:
	rm -f *.o compile compilep