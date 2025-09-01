import ctypes
import numpy as np
import os

# 1. Load Shared Library
# Use the full path to the shared library
lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libsvd.so')
svd_lib = ctypes.CDLL(lib_path)

# 2. Define Argument and Return Types
# The function signature in C++ is:
# extern "C" void calculate_svd_c(const double* input_matrix_flat, int rows, int cols,
#                                   double* output_s_flat, double* output_u_flat, double* output_v_flat);
calculate_svd_c = svd_lib.calculate_svd_c
calculate_svd_c.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
]
calculate_svd_c.restype = None

# 3. Prepare Data
# Create a sample Python matrix
input_matrix = np.array([
    [4.0, 1.0, 1.0],
    [1.0, 2.0, 3.0],
    [1.0, 3.0, 6.0]
], dtype=np.float64)

rows, cols = input_matrix.shape

# Create empty NumPy arrays to receive the output
output_s = np.empty_like(input_matrix)
output_u = np.empty_like(input_matrix)
output_v = np.empty_like(input_matrix)

# 4. Call the C++ Function
calculate_svd_c(input_matrix, rows, cols, output_s, output_u, output_v)

# 5. Print Results
print("Original Matrix (A):")
print(input_matrix)
print("\nSingular Values (S):")
print(np.diag(output_s)) # Print S as a 1D array of singular values
print("\nMatrix S:")
print(output_s)
print("\nMatrix U:")
print(output_u)
print("\nMatrix V:")
print(output_v)

# Verification (optional)
# Reconstruct the original matrix from U, S, V
# Note: In SVD, A = U * S * V^T
reconstructed_matrix = output_u @ output_s @ output_v.T
print("\nReconstructed Matrix (U * S * V^T):")
print(reconstructed_matrix)

print("\nVerification (A is close to U * S * V^T):")
print(np.allclose(input_matrix, reconstructed_matrix))
