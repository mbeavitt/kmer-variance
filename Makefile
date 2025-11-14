CC = clang
CFLAGS = -ggdb3 -O3 -march=native -fPIC
LDFLAGS = -shared

SRC = kmer_variance/kmer_variance_lib.c
TARGET = kmer_variance/libkmer_variance.so

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -f $(TARGET)
