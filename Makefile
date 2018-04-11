CC		= mpicc

CFLAGS		= -std=c99

CLIBS		= -lm 

all: gs

gs: 
	$(CC) $(CFLAGS) -o gs source/gs.c  
clean:
	rm -f gs
