all: mpi mpi_omp

mpi: mainMPI.cpp
	mpicxx -O3 -std=c++11 mainMPI.cpp -o mainMPI

mpi_omp: mainMPI_OMP.cpp
	mpicxx -O3 -fopenmp -std=c++11 mainMPI_OMP.cpp -o mainMPI_OMP

clean:
	rm -rf mainMPI mainMPI_OMP
