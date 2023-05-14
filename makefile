CC = gcc
OBJS = Main.o allvars.o funzioni.o MainFunctions.o
CFLAGS = -Wall -Wextra -O3 -ggdb3 -fno-omit-frame-pointer
all: Main
Main: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o Main -lfftw3 -lm
Main.o: Main.c
	$(CC) $(CFLAGS) -c Main.c
allvars.o: allvars.c allvars.h
	$(CC) $(CFLAGS) -c allvars.c
funzioni.o: funzioni.c funzioni.h
	$(CC) $(CFLAGS) -c funzioni.c
MainFunctions.o: MainFunctions.c MainFunctions.h
	$(CC) $(CFLAGS) -c MainFunctions.c
clean:
	rm -f $(OBJS) $(FOBJS) Main
