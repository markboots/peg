INCLUDEPATH=/home/mjb865/include
CC=g++
CFLAGS=-c -Wall -I$(INCLUDEPATH)
LDFLAGS=-lgsl
SOURCES=main.cpp PEG.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) pegSerial

pegSerial: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
