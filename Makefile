CC = clang
#CFLAGS = -O0 -march=native -Wall -ggdb3
CFLAGS = -O3 -march=native -Wall
TARGET = main
SRC = kmer_variance.c

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)
