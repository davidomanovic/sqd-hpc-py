# sqd-hpc-py

Python bindings for the Qiskit **SQD-HPC** addon and the RIKEN **SBD** (sample-based diagonalization) eigensolver.

This lets you do, from Python:

1. run a circuit on IBM / Qiskit,
2. grab the measured bitstrings,
3. feed them to a single Python function,
4. under the hood it runs the C++ SQD-HPC postprocessing **and** the SBD diagonalizer,

## Requirements

- C++17 compiler
- CMake ≥ 3.20
- MPI (e.g. OpenMPI)
- BLAS + LAPACK (e.g. OpenBLAS)
- OpenMP
- Python ≥ 3.9

## Easy installation
This repo uses submodules

```bash
git clone --recursive https://github.com/you/sqd-hpc-py.git
cd sqd-hpc-py
# if you forgot --recursive:
git submodule update --init --recursive
```

Change directory into the repo and install simply:
```
pip install .
```

## Installing on a HPC cluster
The submodules work on having 32-bit LAPACK integers. Some clusters have by default BLAS/MKL that are ILP64 (64-bit ints). To ensure you have the right environment do the following:

1. Load toolchain
```
module purge
module load GCC/11.2.0
module load OpenMPI/4.1.1-GCC-11.2.0
module load FlexiBLAS/3.0.4-GCC-11.2.0
```

2. Build and install
```
CMAKE_PREFIX_PATH=$EBROOTFLEXIBLAS \
CMAKE_ARGS="-DBLA_VENDOR=FlexiBLAS -DBLA_SIZEOF_INTEGER=4" \
pip install .
```

## Usage

See `tutorials/test.ipynb`
