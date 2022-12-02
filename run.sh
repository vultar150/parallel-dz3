#!/bin/bash

mkdir -p output 
mkdir -p err

PI=3.14159265358979323846

for i in 10 20 40
do
for N in 128 256 512
do
mpisubmit.pl -p $i -w 00:15 --stdout output/mainoutP$i-N$N-L1.out --stderr err/mainerrP$i-N$N.err a.out -- 1.0 1.0 1.0 $N $N $N
mpisubmit.pl -p $i -w 00:15 --stdout output/mainoutP$i-N$N-LPI.out --stderr err/mainerrP$i-N$N.err a.out -- $PI $PI $PI $N $N $N
done
done
