CXX = g++
MPIXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O2
LDLIBS = -lblas

HDRS = Burgers2P.h Burgers.h Model.h Helpers.h
SRC = main.cpp Burgers2P.cpp Burgers.cpp Model.cpp Helpers.cpp
OBJS = $(SRC:.cpp=.o)

default: compile

%.o: %.cpp $(HDRS)
	$(MPIXX) $(CXXFLAGS) -o $@ -c $<

compile: $(OBJS)
	$(MPIXX) -o $@ $^ $(LDLIBS)

all: compile

diff: compile
	./compile 0 0 0 1 10 10 1 0

advx: compile
	./compile 1 0 0 0 10 10 1 0

advy: compile
	./compile 0 1 0 0 10 10 1 0

burg: compile
	./compile 1.0 0.5 1.0 0.02 10 10 1 0

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
	rm -f *.o compile