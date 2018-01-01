all:
	gcc -o simplex simplex.c -Wall -pedantic -fopenmp -llapack -lblas

clean:
	rm -rf simplex;
