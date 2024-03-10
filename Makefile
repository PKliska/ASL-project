CC=gcc

bin/baseline: src/main.c src/baseline_simulation.c
	mkdir bin
	$(CC) -O3 -Wall -Wextra -o bin/baseline src/main.c src/baseline_simulation.c -Iinclude
