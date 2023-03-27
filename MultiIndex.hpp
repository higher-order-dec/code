/*
MultiIndex.hpp provides functions for handling multi-indices.
An n-dimensional multi-index is an n-tuple of non-negative integers.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _MULTIINDEX_HPP_INCLUDED_
#define _MULTIINDEX_HPP_INCLUDED_

#include "GFD/Types/Types.hpp"
#include "GFD/Types/Vector.hpp"
#include "GFD/Types/Buffer.hpp"
#include <vector>
#include <array>

namespace gfd {
	typedef std::vector<uint> mi_t;
	void printMultiIndex(const mi_t& multiIndex);
	std::ostream& operator<< (std::ostream& out, const mi_t& multiIndex);
	uint binCoef(uint n, uint k);
	uint binCoef(const mi_t& n, const mi_t& k);
	unsigned long long factorial(uint n);
	unsigned long long factorial(const mi_t& multiIndex);
	uint multiIndexCardinality(uint n, uint k);
	void getMultiIndices(uint n, uint k, Buffer<mi_t>& multiIndexList);
	void multiIndexTest(uint n, uint k);
	void nonzeroIndices(const mi_t& multiIndex, Buffer<uint>& nzIndices);
	bool contains(const mi_t& multiIndex, uint element);
	bool equals(const mi_t& mi1, const mi_t& mi2);
	uint sum(const mi_t& multiIndex);
	uint indexOf(const mi_t& multiIndex, uint val);
	mi_t subtract(const mi_t& mi1, const mi_t& mi2);
	mi_t insert(const mi_t& multiIndex, uint index, uint value);
	mi_t addZeros(const mi_t& multiIndex, uint index1, uint index2);
	template <class T> double pow(const T& base, const mi_t& exponent);
	template <class T> T permute(const T& list, const mi_t& permutation);
	mi_t permuteOtherWay(const mi_t& multiIndex, const mi_t& permutation);
	uint getPermutation(const Buffer<uint>& buf, const std::vector<uint>& vec);
	double getEdgePermutationSign(uint permutationIndex);
	double getFacePermutationSign(uint permutationIndex, const mi_t& edge);
	double getFacePermutationSign(uint permutationIndex);
	uint getEdgeIndex(const mi_t& multiIndex);
	uint getFaceIndex(const mi_t& multiIndex);
	void getSmallerMultiIndices(const mi_t& multiIndex, std::vector<mi_t>& multiIndices);

	//some special multi-indices

	const std::array<mi_t, 2> edgeNodes{ { {1,0}, {0,1} } };
	const mi_t edgeEdge{ {1,1} };
	const std::array<mi_t, 3> faceNodes{ { {1,0,0}, {0,1,0}, {0,0,1} } };
	const std::array<mi_t, 3> faceEdges{ { {1,1,0}, {1,0,1}, {0,1,1} } };
	const mi_t faceFace{ {1,1,1} };
	const std::array<mi_t, 4> bodyNodes{ { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1} } };
	const std::array<mi_t, 6> bodyEdges{ { {1,1,0,0}, {1,0,1,0}, {1,0,0,1}, {0,1,1,0}, {0,1,0,1}, {0,0,1,1} } };
	const std::array<mi_t, 4> bodyFaces{ { {1,1,1,0}, {1,1,0,1}, {1,0,1,1}, {0,1,1,1} } };
	const mi_t bodyBody{ {1,1,1,1} };
	const std::array<mi_t, 2> edgePermutations{ { {0,1}, {1,0} } };
	const std::array<mi_t, 6> facePermutations{ { {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0} } };
}

#endif //_MULTIINDEX_HPP_INCLUDED_