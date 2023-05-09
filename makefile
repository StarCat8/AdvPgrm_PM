CC = gcc
OBJS = Main.o allvars.o funzioni.o MainFunctions.o
CFLAGS = -Wall -Wextra
all: Main
Main: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o Main -lm
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
