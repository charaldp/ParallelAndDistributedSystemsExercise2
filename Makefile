TARGET=vptree
CC=mpicc
CFLAGS=-Wall -O3 -std=gnu99 -I.
OBJ=vptree.o mpiFindMedian.o
DEPS=mpiFindMedian.h

.PHONY: default all clean

default: $(TARGET)
all: default

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -lm

.PRECIOUS: $(TARGET) $(OBJ)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -lm

clean:
	$(RM) *.o *~
