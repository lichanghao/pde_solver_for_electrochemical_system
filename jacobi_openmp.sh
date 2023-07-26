#! /bin/bash
#
g++ -O3 -std=c++11 -mcmodel=medium -fopenmp jacobi_openmp.cpp -o jacobi_openmp
# g++ -c -Wall -fopenmp jacobi_openmp.cpp
# if [ $? -ne 0 ]; then
#   echo "Compile error."
#   exit
# fi
# #
# g++ -fopenmp -o jacobi_openmp jacobi_openmp.o -lm
# if [ $? -ne 0 ]; then
#   echo "Load error."
#   exit
# fi
# #
# rm jacobi_openmp.o
# #
# mv jacobi_openmp $HOME/bincpp
# #
# echo "Normal end of execution."
