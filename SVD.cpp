// SVD.cpp : Singular Value Decomposition demo (refactored to be standalone and portable)
// Original math/algorithm preserved; only build-simplification edits applied.

#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <cstdlib>

#include <boost/multiprecision/cpp_bin_float.hpp>

using namespace boost::multiprecision;
using namespace std;

typedef number<backends::cpp_bin_float<4096, backends::digit_base_2, void, boost::int16_t, -16382, 16383>, et_off>  cpp_bin_float_4096;

// Forward declarations (order kept; math unchanged)
template<typename Arg = cpp_bin_float_4096> void matrix_transpose(std::vector<std::vector<Arg>> matrix1, std::vector<std::vector<Arg>>& matrix2);
template<typename Arg = cpp_bin_float_4096> void matrix_by_matrix(std::vector<std::vector<Arg>> matrix1, std::vector<std::vector<Arg>>& matrix2, std::vector<std::vector<Arg>>& matrix3);
template<typename Arg = cpp_bin_float_4096> void get_hermitian_matrix(std::vector<Arg> eigenvector, std::vector<std::vector<Arg>>& h_matrix);
template<typename Arg = cpp_bin_float_4096> void get_hermitian_matrix_inverse(std::vector<Arg> eigenvector, std::vector<std::vector<Arg>>& ih_matrix);
template<typename Arg = cpp_bin_float_4096> void jordan_gaussian_transform(std::vector<std::vector<Arg>> matrix, std::vector<Arg>& eigenvector);
template<typename Arg = cpp_bin_float_4096> void get_inverse_diagonal_matrix(std::vector<std::vector<Arg>> matrix, std::vector<std::vector<Arg>>& inv_matrix);
template<typename Arg = cpp_bin_float_4096> void get_reduced_matrix(std::vector<std::vector<Arg>> matrix, std::vector<std::vector<Arg>>& r_matrix, std::size_t new_size);
template<typename Arg = cpp_bin_float_4096> void generate_matrix(std::vector<std::vector<cpp_bin_float_4096>>& matrix, std::size_t rows, std::size_t cols);
template<typename Arg = cpp_bin_float_4096> void print_matrix(std::vector<std::vector<cpp_bin_float_4096>> matrix);
template<typename Arg = cpp_bin_float_4096> void compute_evd(std::vector<std::vector<Arg>> matrix, std::vector<Arg>& eigenvalues, std::vector<std::vector<Arg>>& eigenvectors, std::size_t eig_count);

// SVD driver (unchanged math; removed unused matrix_product1 & u_1)
template<typename Arg = cpp_bin_float_4096>
void svd(std::vector<std::vector<Arg>> matrix, std::vector<std::vector<Arg>>& s,
	std::vector<std::vector<Arg>>& u, std::vector<std::vector<Arg>>& v)
{
	std::vector<std::vector<Arg>> matrix_t;
	matrix_transpose(matrix, matrix_t);

	// Removed unused matrix_product1
	std::vector<std::vector<Arg>> matrix_product2;
	matrix_by_matrix(matrix_t, matrix, matrix_product2);

	std::vector<std::vector<Arg>> v_1; // u_1 removed (unused)
	std::vector<Arg> eigenvalues;
	compute_evd(matrix_product2, eigenvalues, v_1, 0);

	matrix_transpose(v_1, v);

	s.resize(matrix.size());
	for (std::size_t index = 0; index < eigenvalues.size(); index++)
	{
		s[index].resize(eigenvalues.size());
		s[index][index] = eigenvalues[index];
	}

	std::vector<std::vector<Arg>> s_inverse;
	get_inverse_diagonal_matrix(s, s_inverse);

	std::vector<std::vector<Arg>> av_matrix;
	matrix_by_matrix(matrix, v, av_matrix);
	matrix_by_matrix(av_matrix, s_inverse, u);
}

std::vector<std::vector<cpp_bin_float_4096>> matrix_i;

template<typename Arg>
void compute_evd(std::vector<std::vector<Arg>> matrix,
	std::vector<Arg>& eigenvalues, std::vector<std::vector<Arg>>& eigenvectors, std::size_t eig_count)
{
	std::size_t m_size = matrix.size();
	std::vector<Arg> vec; vec.resize(m_size);
	std::generate(vec.begin(), vec.end(), []() {
		return rand() % 100;
	});

	if (eigenvalues.size() == 0 && eigenvectors.size() == 0)
	{
		eigenvalues.resize(m_size);
		eigenvectors.resize(eigenvalues.size());
		matrix_i = matrix;
	}

	std::vector<std::vector<Arg>> m; m.resize(m_size);
	for (std::size_t row = 0; row < m_size; row++)
		m[row].resize(100);

	Arg lambda_old = 0;

	std::size_t index = 0; bool is_eval = false;
	while (is_eval == false)
	{
		for (std::size_t row = 0; row < m_size && (index % 100) == 0; row++)
			m[row].resize(m[row].size() + 100);

		for (std::size_t row = 0; row < m_size; row++)
		{
			m[row][index] = 0;
			for (std::size_t col = 0; col < m_size; col++)
				m[row][index] += matrix[row][col] * vec[col];
		}

		for (std::size_t col = 0; col < m_size; col++)
			vec[col] = m[col][index];

		if (index > 0)
		{
			Arg lambda = ((index - 1) > 0) ? \
				(m[0][index] / m[0][index - 1]) : m[0][index];
			is_eval = (boost::multiprecision::fabs(lambda - lambda_old) < /*10e-15*/10e-10) ? true : false;

			eigenvalues[eig_count] = lambda; lambda_old = lambda;
		}

		index++;
	}

	std::vector<std::vector<Arg>> matrix_new;

	if (m_size > 1)
	{
		std::vector<std::vector<Arg>> matrix_target;
		matrix_target.resize(m_size);

		for (std::size_t row = 0; row < m_size; row++)
		{
			matrix_target[row].resize(m_size);
			for (std::size_t col = 0; col < m_size; col++)
				matrix_target[row][col] = (row == col) ? \
				(matrix[row][col] - eigenvalues[eig_count]) : matrix[row][col];
		}

		std::vector<Arg> eigenvector;
		jordan_gaussian_transform(matrix_target, eigenvector);

		std::vector<std::vector<Arg>> hermitian_matrix;
		get_hermitian_matrix(eigenvector, hermitian_matrix);

		std::vector<std::vector<Arg>> ha_matrix_product;
		matrix_by_matrix(hermitian_matrix, matrix, ha_matrix_product);

		std::vector<std::vector<Arg>> inverse_hermitian_matrix;
		get_hermitian_matrix_inverse(eigenvector, inverse_hermitian_matrix);

		std::vector<std::vector<Arg>> iha_matrix_product;
		matrix_by_matrix(ha_matrix_product, inverse_hermitian_matrix, iha_matrix_product);

		get_reduced_matrix(iha_matrix_product, matrix_new, m_size - 1);
	}

	if (m_size <= 1)
	{
		for (std::size_t index = 0; index < eigenvalues.size(); index++)
		{
			Arg lambda = eigenvalues[index];
			std::vector<std::vector<Arg>> matrix_target;
			matrix_target.resize(matrix_i.size());

			for (std::size_t row = 0; row < matrix_i.size(); row++)
			{
				matrix_target[row].resize(matrix_i.size());
				for (std::size_t col = 0; col < matrix_i.size(); col++)
					matrix_target[row][col] = (row == col) ? \
					(matrix_i[row][col] - lambda) : matrix_i[row][col];
			}

			eigenvectors.resize(matrix_i.size());
			jordan_gaussian_transform(matrix_target, eigenvectors[index]);

			Arg eigsum_sq = 0;
			for (std::size_t v = 0; v < eigenvectors[index].size(); v++)
				eigsum_sq += boost::multiprecision::pow(eigenvectors[index][v], 2.0);

			for (std::size_t v = 0; v < eigenvectors[index].size(); v++)
				eigenvectors[index][v] /= boost::multiprecision::sqrt(eigsum_sq);

			eigenvalues[index] = boost::multiprecision::sqrt(eigenvalues[index]);
		}

		return;
	}

	compute_evd(matrix_new, eigenvalues, eigenvectors, eig_count + 1);

	return;
}

template<typename Arg>
void matrix_transpose(std::vector<std::vector<Arg>> matrix1,
	std::vector<std::vector<Arg>>& matrix2)
{
	matrix2.resize(matrix1.size());
	for (std::size_t row = 0; row < matrix1.size(); row++)
	{
		matrix2[row].resize(matrix1[row].size());
		for (std::size_t col = 0; col < matrix1[row].size(); col++)
			matrix2[row][col] = matrix1[col][row];
	}
}

template<typename Arg>
void matrix_by_matrix(std::vector<std::vector<Arg>> matrix1,
	std::vector<std::vector<Arg>>& matrix2, std::vector<std::vector<Arg>>& matrix3)
{
	matrix3.resize(matrix1.size());
	for (std::size_t row = 0; row < matrix1.size(); row++)
	{
		matrix3[row].resize(matrix1[row].size());
		for (std::size_t col = 0; col < matrix1[row].size(); col++)
		{
			matrix3[row][col] = 0.00;
			for (std::size_t k = 0; k < matrix1[row].size(); k++)
				matrix3[row][col] += matrix1[row][k] * matrix2[k][col];
		}
	}
}

template<typename Arg>
void get_hermitian_matrix(std::vector<Arg> eigenvector,
	std::vector<std::vector<Arg>>& h_matrix)
{
	h_matrix.resize(eigenvector.size());
	for (std::size_t row = 0; row < eigenvector.size(); row++)
		h_matrix[row].resize(eigenvector.size());

	h_matrix[0][0] = 1.0 / eigenvector[0];
	for (std::size_t row = 1; row < eigenvector.size(); row++)
		h_matrix[row][0] = -eigenvector[row] / eigenvector[0];

	for (std::size_t row = 1; row < eigenvector.size(); row++)
		h_matrix[row][row] = 1;
}

template<typename Arg>
void get_hermitian_matrix_inverse(std::vector<Arg> eigenvector,
	std::vector<std::vector<Arg>>& ih_matrix)
{
	ih_matrix.resize(eigenvector.size());
	for (std::size_t row = 0; row < eigenvector.size(); row++)
		ih_matrix[row].resize(eigenvector.size());

	ih_matrix[0][0] = eigenvector[0];
	for (std::size_t row = 1; row < eigenvector.size(); row++)
		ih_matrix[row][0] = -eigenvector[row];

	for (std::size_t row = 1; row < eigenvector.size(); row++)
		ih_matrix[row][row] = 1;
}

template<typename Arg>
void jordan_gaussian_transform(
	std::vector<std::vector<Arg>> matrix, std::vector<Arg>& eigenvector)
{
	const Arg eps = 10e-6; bool eigenv_found = false;
	for (std::size_t s = 0; s < matrix.size() - 1 && !eigenv_found; s++)
	{
		std::size_t col = s; Arg alpha = matrix[s][s];
		while (col < matrix[s].size() && alpha != 0 && alpha != 1)
			matrix[s][col++] /= alpha;

		for (std::size_t col2 = s; col2 < matrix[s].size() && !alpha; col2++)
			std::swap(matrix[s][col2], matrix[s + 1][col2]);

		for (std::size_t row = 0; row < matrix.size(); row++)
		{
			Arg gamma = matrix[row][s];
			for (std::size_t col3 = s; col3 < matrix[row].size() && row != s; col3++)
				matrix[row][col3] = matrix[row][col3] - matrix[s][col3] * gamma;
		}
	}

	std::size_t row = 0;
	while (row < matrix.size())
		eigenvector.push_back(-matrix[row++][matrix.size() - 1]);

	eigenvector[eigenvector.size() - 1] = 1;
}

template<typename Arg>
void get_inverse_diagonal_matrix(std::vector<std::vector<Arg>> matrix,
	std::vector<std::vector<Arg>>& inv_matrix)
{
	inv_matrix.resize(matrix.size());
	for (std::size_t index = 0; index < matrix.size(); index++)
	{
		inv_matrix[index].resize(matrix[index].size());
		inv_matrix[index][index] = 1.0 / matrix[index][index];
	}
}

template<typename Arg>
void get_reduced_matrix(std::vector<std::vector<Arg>> matrix,
	std::vector<std::vector<Arg>>& r_matrix, std::size_t new_size)
{
	if (new_size > 1)
	{
		r_matrix.resize(new_size);
		std::size_t index_d = matrix.size() - new_size;
		std::size_t row = index_d, row_n = 0;
		while (row < matrix.size())
		{
			r_matrix[row_n].resize(new_size);
			std::size_t col = index_d, col_n = 0;
			while (col < matrix.size())
				r_matrix[row_n][col_n++] = matrix[row][col++];

			row++; row_n++;
		}
	}

	else if (new_size == 1)
	{
		r_matrix.resize(new_size);
		r_matrix[0].resize(new_size);
		r_matrix[0][0] = matrix[1][1];
	}
}

template<typename Arg>
void generate_matrix(std::vector<std::vector<cpp_bin_float_4096>>& \
	matrix, std::size_t rows, std::size_t cols)
{
	matrix.resize(rows);
	for (std::size_t row = 0; row < matrix.size(); row++)
	{
		matrix[row].resize(cols);
		for (std::size_t col = 0; col < matrix[row].size(); col++)
			matrix[row][col] = std::rand() % 10;
	}
}

template<typename Arg>
void print_matrix(std::vector<std::vector<cpp_bin_float_4096>>	matrix)
{
	for (std::size_t row = 0; row < matrix.size(); row++)
	{
		for (std::size_t col = 0; col < matrix[row].size(); col++)
			std::cout << std::setprecision(5) << std::fixed << matrix[row][col] << " ";

		std::cout << "\n";
	}

	std::cout << "\n";
}

#ifndef PYTHON_WRAPPER
int main(int argc, char* argv[])
{
	std::srand((unsigned int)std::time(nullptr)); // Seed once at start
	
	std::size_t matrix_size = 0;
	std::vector<std::vector<cpp_bin_float_4096>> u, v;
	std::vector<std::vector<cpp_bin_float_4096>> matrix, s;
	std::cout << "Singular Value Decomposition (SVD):\n\n";
	
	std::cout << "Enter size of matrix N = (50x50 max): "; std::cin >> matrix_size;

	if (matrix_size <= 50)
	{
		generate_matrix(matrix, matrix_size, matrix_size);

		std::cout << "\nA = \n"; print_matrix(matrix);

		svd(matrix, s, u, v);

		std::cout << "\nS = \n"; print_matrix(s);
		std::cout << "\nU = \n"; print_matrix(u);
		std::cout << "\nV = \n"; print_matrix(v);
	}

	else std::cout << "Wrong matrix size... (matrix decomposition recommended)";

	std::cin.get();

	return 0;
}
#endif

extern "C" void calculate_svd_c(
    const double* input_matrix_flat, int rows, int cols,
    double* output_s_flat, double* output_u_flat, double* output_v_flat
) {
    // 1. Input Conversion
    std::vector<std::vector<cpp_bin_float_4096>> matrix(rows, std::vector<cpp_bin_float_4096>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = input_matrix_flat[i * cols + j];
        }
    }

    // 2. Call C++ SVD
    std::vector<std::vector<cpp_bin_float_4096>> s, u, v;
    svd(matrix, s, u, v);

    // 3. Output Conversion
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            output_s_flat[i * cols + j] = s[i][j].convert_to<double>();
            output_u_flat[i * cols + j] = u[i][j].convert_to<double>();
            output_v_flat[i * cols + j] = v[i][j].convert_to<double>();
        }
    }
}


