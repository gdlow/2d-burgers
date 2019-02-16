default: burgers_exe

main.o: main.cpp
	g++ -std=c++11 -Wall -O2 -o main.o -c main.cpp

Burgers.o: Burgers.cpp Burgers.h
	g++ -std=c++11 -Wall -O2 -o Burgers.o -c Burgers.cpp

Model.o: Model.cpp Model.h
	g++ -std=c++11 -Wall -O2 -o Model.o -c Model.cpp

ParseException.o: ParseExcepton.h
	g++ -std=c++11 -Wall -O2 -o ParseException.o -c ParseException.h	

burgers_exe: main.o Burgers.o Model.o
	g++ -o burgers_exe main.o Burgers.o Model.o

#invalid argument exception should be thrown
Test_0: burgers_exe
	./burgers_exe 0 0 0 1 

Test_1: burgers_exe
	./burgers_exe 0 0 0 1 10 10 1

Test_2: burgers_exe
	./burgers_exe 1 0 0 0 10 10 1

Test_3: burgers_exe
	./burgers_exe 0 1 0 0 10 10 1

Test_4: burgers_exe
	./burgers_exe 1.0 0.5 1.0 0.02 10 10 1

all: burgers_exe

.PHONY: clean
clean:
	rm -f *.o burgers_exe