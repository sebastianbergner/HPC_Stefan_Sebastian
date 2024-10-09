CFLAGS = -std=gnu11 -O2 -Werror -Wall -Wpedantic -Wextra

.PHONY: all
all: monte_carlo_serial monte_carlo_parallel

.PHONY: clean
clean:
	$(RM) monte_carlo_serial monte_carlo_parallel

monte_carlo_serial: monte_carlo_serial.c

monte_carlo_parallel: monte_carlo_parallel.c
	gcc $(CFLAGS) monte_carlo_parallel.c -o monte_carlo_parallel -pthread
