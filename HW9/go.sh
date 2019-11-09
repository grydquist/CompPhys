#!/bin/bash
rm -f a.out
rm -f *.gch
#rm -f *txt
#DBGFLAG="-ffree-line-length-none -g -fbacktrace -fcheck=all -pedantic -Wall -Wextra -W -Wno-unused-function -fopenmp"
#OPTFLAG="-ffree-line-length-none -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe"
g++ -std=c++11  HW9.cpp
rm -f *.gch *.o
./a.out