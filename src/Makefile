INCLUDEPATH=/home/mjb865/include
LIBPATH=/home/mjb865/lib
CC=g++
CFLAGS=-c -g -Wall -I$(INCLUDEPATH)
LDFLAGS=-L$(LIBPATH) -lgsl -lgslcblas
SOURCES=mainSerial.cpp mainMPI.cpp PEG.cpp PESolver.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) pegSerial pegMPI

pegSerial: mainSerial.o PEG.o PESolver.o
	$(CC) $(LDFLAGS) mainSerial.o PEG.o PESolver.o -o $@

pegMPI: mainMPI.o PEG.o PESolver.o
	$(CC) $(LDFLAGS) mainMPI.o PEG.o PESolver.o -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
