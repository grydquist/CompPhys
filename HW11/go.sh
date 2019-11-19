#!/bin/bash
# Remove executable
rm -f hw11.exe

# Compile the source
FLAGS="-o"
g++ $FLAGS hw11.exe HW11.cpp


# Run the code
./hw11.exe
