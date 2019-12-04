#!/bin/bash
# Remove executable
rm -f hw11.exe

# Compile the source
FLAGS="-o"
g++ $FLAGS hw11.exe hw11.cpp


# Run the code
./hw11.exe
