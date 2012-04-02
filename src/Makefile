CC=g++
CFLAGS=-c -Wall
LDFLAGS=-lgsl
SOURCES=main.cpp PEG.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) pegSerial

pegSerial: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
