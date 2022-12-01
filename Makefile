np?=1

all: compile

compile: main.cpp
	mpicxx main.cpp

run: 
	mpirun -np $(np) ./a.out

clean:
	rm -rf a.out
