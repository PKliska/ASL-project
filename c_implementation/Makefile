CC=gcc
CFLAGS=-O3 -Wall -Wextra -Iinclude

.PHONY: all clean benchmark

all: bin/baseline

bin/baseline: src/main.c src/baseline_simulation.c
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $^

benchmark:
	$(CC) $(CFLAGS) -o bin/runtime benchmarking/runtime.c src/baseline_simulation.c
	./bin/runtime

clean:
	rm -rf bin
