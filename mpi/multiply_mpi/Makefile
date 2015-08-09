CC = mpicc
CFLAGS = -Wall -g

run:
	make multiply
	mpirun -np 4 ./multiply
clean:
	rm -r multiply
