# make file
all:
	g++ -g -pthread -c coarse_alignment.c -o coarse_alignment.o
	g++ -g coarse_alignment.o -o test
clean:
	rm -rf *.o
	rm test

