#!/bin/bash
# Remove executable
rm -f hw10.exe

# Compile the source
FLAGS="-o"
g++ $FLAGS hw10.exe hw10.cpp

# Move executable back out
mv hw10.exe ../

# Move back out
cd ../

# Run the code
./hw10.exe
