CC = mpicc
CFLAGS = -Wall -g

run:
	make clean
	make multiply
	mpirun -np 4 ./multiply
clean:
	rm -r multiply
