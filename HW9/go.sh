#!/bin/bash
# Remove executable
rm -f hw9.exe

# Compile the source
FLAGS="-o"
g++ $FLAGS hw9.exe HW9.cpp

# Move executable back out
mv hw9.exe ../

# Move back out
cd ../

# Run the code
./hw9.exe
