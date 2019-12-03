#!/bin/bash
# Remove executable
rm -f hw12.exe

# Compile the source
FLAGS="-o"
g++ $FLAGS hw12.exe hw12.cpp


# Run the code
./hw12.exe
