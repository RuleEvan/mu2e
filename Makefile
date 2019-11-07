CC=gcc -m64 
CFLAGS=-c -Wall -lm -ldl

all: Mu2E

Mu2E: main.o wfn.o charge.o harmonic.o file_io.o angular.o romberg.o
	$(CC) main.o wfn.o charge.o harmonic.o file_io.o angular.o romberg.o -o Mu2E -lm -ldl -lgsl -lgslcblas

main.o: main.c
	$(CC) $(CFLAGS) main.c

wfn.o: wfn.c
	$(CC) $(CFLAGS) wfn.c

harmonic.o: harmonic.c
	$(CC) $(CFLAGS) harmonic.c

charge.o: charge.c
	$(CC) $(CFLAGS) charge.c

file_io.o: file_io.c
	$(CC) $(CFLAGS) file_io.c

angular.o: angular.c
	$(CC) $(CFLAGS) angular.c

romberg.o: romberg.c
	$(CC) $(CFLAGS) romberg.c

clean:
	rm -rf *.o Mu2E
