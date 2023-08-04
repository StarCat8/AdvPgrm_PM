CC = gcc
OBJS1 = InitCondi.o allvars.o funzioni.o
OBJS2 = Main.o allvars.o funzioni.o MainFunctions.o
CFLAGS = -Wall -Wextra -O3 -ggdb3 -fno-omit-frame-pointer -fopenmp
all: Init Main
Init: $(OBJS1)
	$(CC) $(CFLAGS) $(OBJS1) -o Init -lfftw3 -lm
Main: $(OBJS2)
	$(CC) $(CFLAGS) $(OBJS2) -o Main -lfftw3 -lm
InitCondi.o: InitCondi.c
	$(CC) $(CFLAGS) -c InitCondi.c
Main.o: Main.c
	$(CC) $(CFLAGS) -c Main.c
allvars.o: allvars.c allvars.h
	$(CC) $(CFLAGS) -c allvars.c
funzioni.o: funzioni.c funzioni.h
	$(CC) $(CFLAGS) -c funzioni.c
MainFunctions.o: MainFunctions.c MainFunctions.h
	$(CC) $(CFLAGS) -c MainFunctions.c
clean:
	rm -f $(OBJS1) $(OBJS2)
