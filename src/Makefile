INCLUDEPATH=/home/mjb865/include
LIBPATH=/home/mjb865/lib
CC=g++
CFLAGS=-c -Wall -I$(INCLUDEPATH)
LDFLAGS=-L$(LIBPATH) -lgsl -lgslcblas
SOURCES=main.cpp PEG.cpp PESolver.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) ../build/pegSerial

../build/pegSerial: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
