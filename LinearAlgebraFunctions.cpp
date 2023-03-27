#include "LinearAlgebraFunctions.hpp"
#include <iostream>
#include <fstream>

namespace gfd {

	//print the matrix to the console
	void printMatrix(const Buffer<Buffer<double>>& mat) {
		for (uint i = 0; i < mat.size(); ++i) {
			for (uint j = 0; j < mat[i].size(); ++j) {
				std::cout << mat[i][j] << ' ';
			}
			std::cout << '\n';
		}
	}

	//print the matrix to the console
	void printMatrix(const Buffer<VectorN>& mat) {
		for (uint i = 0; i < mat.size(); ++i) {
			for (uint j = 0; j < mat[i].size(); ++j) {
				std::cout << mat[i][j] << ' ';
			}
			std::cout << '\n';
		}
	}

	//compute the LU decomposition and the permutation vector (see https://en.wikipedia.org/wiki/LU_decomposition)
	bool decomposeLUP(Buffer<Buffer<double>>& A, Buffer<uint>& P, const double tol) {
		uint i, j, k, imax;
		double maxA, absA;
		Buffer<double> row;
		uint n = A.size();
		P.resize(n);

		for (i = 0; i < n; i++)
			P[i] = i;

		for (i = 0; i < n; i++) {
			maxA = 0.0;
			imax = i;

			for (k = i; k < n; k++)
				if ((absA = std::abs(A[k][i])) > maxA) {
					maxA = absA;
					imax = k;
				}

			if (maxA < tol)
				return false; //matrix is degenerate

			if (imax != i) {
				//pivoting P
				j = P[i];
				P[i] = P[imax];
				P[imax] = j;

				//pivoting rows of A
				row = A[i];
				A[i] = A[imax];
				A[imax] = row;
			}

			for (j = i + 1; j < n; j++) {
				A[j][i] /= A[i][i];

				for (k = i + 1; k < n; k++)
					A[j][k] -= A[j][i] * A[i][k];
			}
		}

		return true;  //decomposition succeeded 
	}

	//solve the linear system using the LU decomposition and the permutation vector (see https://en.wikipedia.org/wiki/LU_decomposition)
	void solveLUP(const Buffer<Buffer<double>>& A, const Buffer<uint>& P, const Buffer<double>& b, Buffer<double>& x) {
		const int n = b.size();
		x.resize(n);
		for (int i = 0; i < n; i++) {
			x[i] = b[P[i]];

			for (int k = 0; k < i; k++)
				x[i] -= A[i][k] * x[k];
		}

		for (int i = n - 1; i >= 0; i--) {
			for (int k = i + 1; k < n; k++)
				x[i] -= A[i][k] * x[k];

			x[i] /= A[i][i];
		}
	}

	//solve the linear system using the LU decomposition and the permutation vector (see https://en.wikipedia.org/wiki/LU_decomposition)
	void solveLUP(const MatrixN& A, const Buffer<uint>& P, const VectorN& b, VectorN& x) {
		int n = b.size();
		x.toVectorN(n);
		for (int i = 0; i < n; i++) {
			x[i] = b[P[i]];

			for (int k = 0; k < i; k++)
				x[i] -= A[i][k] * x[k];
		}

		for (int i = n - 1; i >= 0; i--) {
			for (int k = i + 1; k < n; k++)
				x[i] -= A[i][k] * x[k];

			x[i] /= A[i][i];
		}
	}

	//convert buffer of buffers to MatrixN
	void convertMatrix(const Buffer<Buffer<double>>& source, MatrixN& target) {
		uint n = source.size();
		target.toMatrixN(n);
		uint i, j;
		for (i = 0; i < n; ++i) {
			target[i].toVectorN(n);
			for (j = 0; j < n; ++j) {
				target[i][j] = source[i][j];
			}
		}
	}

	//save matrix to file
	bool saveMatrix(const Buffer<Buffer<double>>& mat, const std::string& path) {
		std::ofstream fs(path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
		if (fs.fail()) return false;
		uint rows = mat.size();
		uint cols = mat[0].size();
		fs.write((char*)&rows, sizeof(uint));
		fs.write((char*)&cols, sizeof(uint));
		for (uint i = 0; i < rows; ++i)
			fs.write((char*)&mat[i][0], cols * sizeof(double));
		fs.close();
		return true;
	}

	//save matrix to file
	bool saveMatrix(const MatrixN& mat, const std::string& path) {
		std::ofstream fs(path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
		if (fs.fail()) return false;
		uint size = mat.size();
		fs.write((char*)&size, sizeof(uint));
		for (uint i = 0; i < size; ++i)
			fs.write((char*)&mat[i].val[0], size * sizeof(double));
		fs.close();
		return true;
	}

	//save matrix to file
	bool saveMatrix(const Buffer<VectorN>& mat, const std::string& path) {
		std::ofstream fs(path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
		if (fs.fail()) return false;
		uint rows = mat.size();
		uint cols = mat[0].size();
		fs.write((char*)&rows, sizeof(uint));
		fs.write((char*)&cols, sizeof(uint));
		for (uint i = 0; i < rows; ++i)
			fs.write((char*)&mat[i][0], cols * sizeof(double));
		fs.close();
		return true;
	}

	//save matrix of vectors to file
	bool saveMatrix(const Buffer<Buffer<VectorN>>& mat, const std::string& path) {
		std::ofstream fs(path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
		if (fs.fail()) return false;
		uint size_i = mat.size();
		uint size_j = mat[0].size();
		uint size_k = mat[0][0].size();
		fs.write((char*)&size_i, sizeof(uint));
		fs.write((char*)&size_j, sizeof(uint));
		fs.write((char*)&size_k, sizeof(uint));
		for (uint i = 0; i < size_i; ++i)
			for (uint j = 0; j < size_j; ++j)
				fs.write((char*)&mat[i][j].val[0], size_k * sizeof(double));
		fs.close();
		return true;
	}

	//load matrix from file
	bool loadMatrix(Buffer<Buffer<double>>& mat, const std::string& path) {
		std::ifstream fs(path.c_str(), std::ios::binary | std::ios::in);
		if (fs.fail()) return false;
		uint rows, cols;
		fs.read((char*)&rows, sizeof(uint));
		fs.read((char*)&cols, sizeof(uint));
		mat.resize(rows);
		for (uint i = 0; i < rows; ++i) {
			mat[i].resize(cols);
			fs.read((char*)&mat[i][0], cols * sizeof(double));
		}
		fs.close();
		return true;
	}

	//load matrix from file
	bool loadMatrix(MatrixN& mat, const std::string& path) {
		std::ifstream fs(path.c_str(), std::ios::binary | std::ios::in);
		if (fs.fail()) return false;
		uint size;
		fs.read((char*)&size, sizeof(uint));
		mat.val.resize(size);
		for (uint i = 0; i < size; ++i) {
			mat[i].val.resize(size);
			fs.read((char*)&mat[i][0], size * sizeof(double));
		}
		fs.close();
		return true;
	}

	//load matrix from file
	bool loadMatrix(Buffer<VectorN>& mat, const std::string& path) {
		std::ifstream fs(path.c_str(), std::ios::binary | std::ios::in);
		if (fs.fail()) return false;
		uint rows, cols;
		fs.read((char*)&rows, sizeof(uint));
		fs.read((char*)&cols, sizeof(uint));
		mat.resize(rows);
		for (uint i = 0; i < rows; ++i) {
			mat[i].toVectorN(cols);
			fs.read((char*)&mat[i][0], cols * sizeof(double));
		}
		fs.close();
		return true;
	}

	//load matrix of vectors from file
	bool loadMatrix(Buffer<Buffer<VectorN>>& mat, const std::string& path) {
		std::ifstream fs(path.c_str(), std::ios::binary | std::ios::in);
		if (fs.fail()) return false;
		uint size_i, size_j, size_k;
		fs.read((char*)&size_i, sizeof(uint));
		fs.read((char*)&size_j, sizeof(uint));
		fs.read((char*)&size_k, sizeof(uint));
		mat.resize(size_i);
		for (uint i = 0; i < size_i; ++i) {
			mat[i].resize(size_j);
			for (uint j = 0; j < size_j; ++j) {
				mat[i][j].val.resize(size_k);
				fs.read((char*)&mat[i][j][0], size_k * sizeof(double));
			}		
		}
		fs.close();
		return true;
	}

	//save buffer to file
	template <class T> bool saveBuffer(const Buffer<T>& buf, const std::string& path) {
		std::ofstream fs(path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
		if (fs.fail()) return false;
		uint size = buf.size();
		fs.write((char*)&size, sizeof(uint));
		fs.write((char*)&buf[0], size * sizeof(T));
		fs.close();
		return true;
	}
	template bool saveBuffer<double>(const Buffer<double>& buf, const std::string& path);
	template bool saveBuffer<uint>(const Buffer<uint>& buf, const std::string& path);

	//load buffer from file
	template <class T> bool loadBuffer(Buffer<T>& buf, const std::string& path) {
		std::ifstream fs(path.c_str(), std::ios::binary | std::ios::in);
		if (fs.fail()) return false;
		uint size;
		fs.read((char*)&size, sizeof(uint));
		buf.resize(size);
		fs.read((char*)&buf[0], size * sizeof(T));
		fs.close();
		return true;
	}
	template bool loadBuffer<double>(Buffer<double>& buf, const std::string& path);
	template bool loadBuffer<uint>(Buffer<uint>& buf, const std::string& path);

	//save a single element to file
	template <class T> bool saveElement(const T& element, const std::string& path) {
		std::ofstream fs(path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
		if (fs.fail()) return false;
		fs.write((char*)&element, sizeof(T));
		fs.close();
		return true;
	}
	template bool saveElement<double>(const double& element, const std::string& path);

	//load a single element from file
	template <class T> bool loadElement(T& element, const std::string& path) {
		std::ifstream fs(path.c_str(), std::ios::binary | std::ios::in);
		if (fs.fail()) return false;
		fs.read((char*)&element, sizeof(T));
		fs.close();
		return true;
	}
	template bool loadElement<double>(double& element, const std::string& path);

	//compute the infinity norm of the matrix
	double computeInfinityNorm(const Buffer<Buffer<double>>& mat) {
		double maxRowSum = 0.0;
		for (uint i = 0; i < mat.size(); ++i) {
			double rowSum = 0.0;
			for (uint j = 0; j < mat[i].size(); ++j)
				rowSum += std::abs(mat[i][j]);
			if (rowSum > maxRowSum)
				maxRowSum = rowSum;
		}
		return maxRowSum;
	}

	//compute the infinity norm of the matrix
	double computeInfinityNorm(const MatrixN& mat) {
		double maxRowSum = 0.0;
		for (uint i = 0; i < mat.size(); ++i) {
			double rowSum = 0.0;
			for (uint j = 0; j < mat[i].size(); ++j)
				rowSum += std::abs(mat[i][j]);
			if (rowSum > maxRowSum)
				maxRowSum = rowSum;
		}
		return maxRowSum;
	}

	//compute a weighted infinity norm of the matrix
	double computeWeightedInfinityNorm(const MatrixN& mat, const Buffer<double>& weights) {
		double maxRowSum = 0.0;
		for (uint i = 0; i < mat.size(); ++i) {
			double rowSum = 0.0;
			for (uint j = 0; j < mat[i].size(); ++j)
				rowSum += std::abs(mat[i][j]) * weights[j];
			if (rowSum > maxRowSum)
				maxRowSum = rowSum;
		}
		return maxRowSum;
	}

	//count the number of nonzero elements in the matrix
	uint numberOfNonzeros(const Buffer<Buffer<double>>& mat) {
		uint nz = 0;
		for (uint i = 0; i < mat.size(); ++i)
			for (uint j = 0; j < mat[i].size(); ++j)
				if (mat[i][j] != 0)
					++nz;
		return nz;
	}
}