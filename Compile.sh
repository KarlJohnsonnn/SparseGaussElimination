#!/bin/bash
#===============================================================================
# Build script for SparseGaussElimination
#===============================================================================

set -e  # Exit on error

# Compiler settings
FC=gfortran
TARGET=SparseGaussElimination.exe
SOURCES="Kind_Mod.f90 mo_unirnk.f90 Sparse_Mod.f90 Main_SparseGaussElimination.f90"

# Compiler flags
FFLAGS_DEBUG="-g -O0 -Wall -Wextra -fbounds-check -fbacktrace"
FFLAGS_RELEASE="-O3 -march=native"

# Default to release build
FFLAGS=$FFLAGS_RELEASE

# Parse command line arguments
if [ "$1" == "debug" ]; then
    FFLAGS=$FFLAGS_DEBUG
    echo "Building in DEBUG mode..."
elif [ "$1" == "clean" ]; then
    echo "Cleaning build artifacts..."
    rm -f *.o *.mod $TARGET
    echo "Clean complete."
    exit 0
else
    echo "Building in RELEASE mode..."
fi

# Compile
echo "Compiling $TARGET..."
$FC $FFLAGS -o $TARGET $SOURCES

if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo "Executable: $TARGET"
else
    echo "Build failed!"
    exit 1
fi
