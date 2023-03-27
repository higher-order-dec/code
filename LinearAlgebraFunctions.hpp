/*
LinearAlgebraFunctions contains functions related to vectors, matrices, and linear algebra.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _LINEARALGEBRAFUNCTIONS_HPP_INCLUDED_
#define _LINEARALGEBRAFUNCTIONS_HPP_INCLUDED_

#include "GFD/Types/Types.hpp"
#include "GFD/Types/Buffer.hpp"
#include "GFD/Types/Matrix.hpp"
#include <string>
#include <unordered_map>

namespace gfd
{
	void printMatrix(const Buffer<Buffer<double>>& mat);
	void printMatrix(const Buffer<VectorN>& mat);
	bool decomposeLUP(Buffer<Buffer<double>>& A, Buffer<uint>& P, const double tol);
	void solveLUP(const Buffer<Buffer<double>>& A, const Buffer<uint>& P, const Buffer<double>& b, Buffer<double>& x);
	void solveLUP(const MatrixN& A, const Buffer<uint>& P, const VectorN& b, VectorN& x);
	void convertMatrix(const Buffer<Buffer<double>>& source, MatrixN& target);
	bool saveMatrix(const Buffer<Buffer<double>>& mat, const std::string& path);
	bool saveMatrix(const MatrixN& mat, const std::string& path);
	bool saveMatrix(const Buffer<VectorN>& mat, const std::string& path);
	bool saveMatrix(const Buffer<Buffer<VectorN>>& mat, const std::string& path);
	bool loadMatrix(Buffer<Buffer<double>>& mat, const std::string& path);
	bool loadMatrix(MatrixN& mat, const std::string& path);
	bool loadMatrix(Buffer<VectorN>& mat, const std::string& path);
	bool loadMatrix(Buffer<Buffer<VectorN>>& mat, const std::string& path);
	template <class T> bool saveBuffer(const Buffer<T>& buf, const std::string& path);
	template <class T> bool loadBuffer(Buffer<T>& buf, const std::string& path);
	template <class T> bool saveElement(const T& element, const std::string& path);
	template <class T> bool loadElement(T& element, const std::string& path);
	double computeInfinityNorm(const Buffer<Buffer<double>>& mat);
	double computeInfinityNorm(const MatrixN& mat);
	double computeWeightedInfinityNorm(const MatrixN& mat, const Buffer<double>& weights);
	uint numberOfNonzeros(const Buffer<Buffer<double>>& mat);
}

#endif //_LINEARALGEBRAFUNCTIONS_HPP_INCLUDED_