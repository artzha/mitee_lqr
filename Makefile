CC=gcc
CFLAGS=-lm -lgsl -lgslcblas -std=c99
DEBUG=-g

DEPS=controller.o
OBJ=main.o controller.o

%.o: %.c $(DEPS) 
	$(CC) $(CFLAGS) -c -o  $@ $<

lqr.o: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean

clean:
	rm -rf $(ODIR)/*.o *~ core *.dSYM *.o
