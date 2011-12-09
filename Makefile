CC=g++
CFLAGS=-Wall -O3 -funroll-loops
LDFLAGS=-O2 -lm
SOURCES=$(wildcard *.cc)
HEADERS=$(wildcard *.h)
OBJECTS=$(addsuffix .o, $(basename ${SOURCES}))
EXECUTABLE=MyBot

#Uncomment the following to enable debugging
#CFLAGS+=-g #-DDEBUG

all: $(OBJECTS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cc $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f ${EXECUTABLE} ${OBJECTS} *.d
	-rm -f debug.txt

.PHONY: all clean

