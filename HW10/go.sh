#!/bin/bash
# Remove executable
rm -f hw10.exe

# Compile the source
FLAGS="-o"
g++ $FLAGS hw10.exe hw10.cpp


# Run the code
./hw10.exe
