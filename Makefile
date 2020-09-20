CC=gcc
CFLAGS=-lm -lgsl -lgslcblas -std=c99 -Wconversion -Wall -Werror -Wextra -pedantic
DEBUG=-g3

DEPS=controller.o
OBJ=main.o controller.o external.o

main: main.o controller.o external.o
	$(CC) $^ $(DEBUG) -o main $(CFLAGS)
	
main.o: main.c controller.h external.h
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $<	

controller.o: controller.c controller.h external.h
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $<
	
external.o: external.c external.h
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $<

.PHONY: clean

clean:
	rm main $(OBJ)
