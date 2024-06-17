CFLAGS=-I$PWD -lm  -Wall  -O3  -march=native

a.out: utils.o main.o
	mpicc -o a.out main.o utils.o $(CFLAGS)

utils.o: utils.c utils.h
	mpicc -c utils.c $(CFLAGS)

main.o: main.c utils.h
	mpicc -c main.c $(CFLAGS)
	
delete:
	rm -f out_*

clean:
	rm -f *~ *.o a.out
