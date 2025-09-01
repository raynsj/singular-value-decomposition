# C++ SVD with Python Wrapper

This project demonstrates a high-precision Singular Value Decomposition (SVD) calculation implemented in C++ using the Boost Multiprecision library. It can be run as a standalone command-line tool or as a high-performance backend for a Python script.

## Prerequisites

- A C++ compiler (e.g., `g++`)
- The Boost C++ libraries (specifically `multiprecision` and `system`)
- Python 3
- The NumPy Python library (`pip install numpy`)

## How to Use

The project can be used in two ways:

### 1. Standalone C++ Executable

This mode generates a random matrix and performs the SVD calculation directly in the terminal.

**Build the executable:**
```bash
make
```
This creates an executable file named `svd`.

**Run the program:**
```bash
./svd
```
The program will prompt you to enter a matrix size (e.g., `3`). It will then generate a random 3x3 matrix and print the resulting U, S, and V matrices.

You can also pipe the input directly:
```bash
echo "3" | ./svd
```

### 2. Python Wrapper (via Shared Library)

This mode allows you to call the high-performance C++ SVD function from a Python script, passing a NumPy array as input.

**Build the shared library:**
```bash
make shared
```
This creates a shared library file named `libsvd.so`.

**Run the Python script:**
```bash
python3 run_svd.py
```
The script will:
1. Define a sample matrix in NumPy.
2. Call the C++ function in `libsvd.so` to perform the SVD.
3. Print the resulting U, S, and V matrices.
4. Verify the result by reconstructing the original matrix (`U * S * V^T`).

### Cleaning Up

To remove all compiled files (`svd` and `libsvd.so`), run:
```bash
make clean
```
