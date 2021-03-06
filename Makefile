# Compilers and flags
CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O3
LDLIBS = -lblas

# Serial variables
DIR_SER = serSrc
HDRS_SER = Burgers.h Model.h
SRC_SER = serialEntryPoint.cpp Burgers.cpp Model.cpp
OBJS_SER = $(addprefix $(DIR_SER)/,$(SRC_SER:.cpp=.o))

# Parallel variables
DIR_PAR = parSrc
HDRS_PAR = Burgers2P.h Model2P.h
SRC_PAR = parallelEntryPoint.cpp Burgers2P.cpp Model2P.cpp
OBJS_PAR = $(addprefix $(DIR_PAR)/,$(SRC_PAR:.cpp=.o))

# Build serial code
$(DIR_SER)/%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

compile: $(OBJS_SER)
	$(CXX) -o $@ $^ $(LDLIBS)

# Build parallel code
$(DIR_PAR)/%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

compilep: $(OBJS_PAR)
	$(CXX) -o $@ $^ $(LDLIBS)

# Serial targets
diff: compile
	./compile 0 0 0 1 10 10 1

advx: compile
	./compile 1 0 0 0 10 10 1

advy: compile
	./compile 0 1 0 0 10 10 1

burg: compile
	./compile 1.0 0.5 1.0 0.02 10 10 1

# Parallel targets (-np = 2)
diffp: compilep
	mpiexec -np 2 ./compilep 0 0 0 1 10 10 1 2 1

advxp: compilep
	mpiexec -np 2 ./compilep 1 0 0 0 10 10 1 2 1

advyp: compilep
	mpiexec -np 2 ./compilep 0 1 0 0 10 10 1 2 1

burgp: compilep
	mpiexec -np 2 ./compilep 1.0 0.5 1.0 0.02 10 10 1 2 1

report: compilep
	mpiexec -np 12 ./compilep 1.0 0.5 1.0 0.02 10 10 1 3 4

# Misc
default: compile

all: compile compilep

.PHONY: clean
clean:
	rm -f $(DIR_SER)/*.o $(DIR_PAR)/*.o compile compilep