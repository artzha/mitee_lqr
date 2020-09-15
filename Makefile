CC=gcc
CFLAGS=-lm -lgsl -lgslcblas -std=c99 -Wconversion -Wall -Werror -Wextra -pedantic
DEBUG=-g3

DEPS=controller.o
OBJ=main.o controller.o

main: main.o
	$(CC) $< $(DEBUG) -o main $(CFLAGS)

main.o: main.c controller.o
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $<
	
controller.o: controller.c controller.h
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $<

.PHONY: clean

clean:
	rm main $(OBJ)
