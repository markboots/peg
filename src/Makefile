INCLUDEPATH=/home/mjb865/include
LIBPATH=/home/mjb865/lib
CC=mpic++
CFLAGS=-c -g -Wall -I$(INCLUDEPATH)
LDFLAGS=-L$(LIBPATH) -lgsl -lgslcblas
SOURCES=mainSerial.cpp mainMPI.cpp PEG.cpp PESolver.cpp PEMainSupport.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) pegSerial pegMPI

pegSerial: mainSerial.o PEMainSupport.o PEG.o PESolver.o
	$(CC) $(LDFLAGS) mainSerial.o PEMainSupport.o PEG.o PESolver.o -o $@

pegMPI: mainMPI.o PEMainSupport.o PEG.o PESolver.o
	$(CC) $(LDFLAGS) mainMPI.o PEMainSupport.o PEG.o PESolver.o -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
