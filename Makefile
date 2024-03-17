CC=gcc

bin/baseline: src/main.c src/baseline_simulation.c
	# "mkdir bin" fails when the bin directory exists. We don't want that so swallow the
	# error with "|| true".
	mkdir bin || true

	$(CC) -O3 -Wall -Wextra -o bin/baseline src/main.c src/baseline_simulation.c -Iinclude
