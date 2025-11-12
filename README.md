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

## Clone

This repo uses submodules:

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

## Usage

See `tutorials/test.ipynb`
