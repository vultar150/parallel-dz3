#!/bin/bash

mkdir -p outputMPI 
mkdir -p errMPI

PI=3.14159265358979323846

for i in 1 4 8 16 32
do
for N in 128 256 512
do
mpisubmit.pl -p $i -w 00:15 --stdout outputMPI/mainMPIoutP$i-N$N-L1.txt --stderr errMPI/mainMPIerrP$i-N$N.txt mainMPI -- 1.0 1.0 1.0 $N $N $N
mpisubmit.pl -p $i -w 00:15 --stdout outputMPI/mainMPIoutP$i-N$N-LPI.txt --stderr errMPI/mainMPIerrP$i-N$N.txt mainMPI -- $PI $PI $PI $N $N $N
done
done
