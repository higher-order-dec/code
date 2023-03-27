#include "MultiIndex.hpp"
#include <iostream>
#include <vector>
#include <numeric>

namespace gfd {

	//print the multi-index to the console
	void printMultiIndex(const mi_t& multiIndex) {
		std::cout << '(';
		for (uint i = 0; i < multiIndex.size() - 1; ++i) {
			std::cout << multiIndex[i] << ',';
		}
		std::cout << multiIndex[multiIndex.size() - 1] << ")";
	}

	//define how multi-indices are printed to the console (std::cout << multiIndex)
	std::ostream& operator<< (std::ostream& out, const mi_t& multiIndex) {
		out << '(';
		for (uint i = 0; i < multiIndex.size() - 1; ++i) {
			out << multiIndex[i] << ',';
		}
		out << multiIndex[multiIndex.size() - 1] << ")";
		return out;
	}

	//returns the binomial coefficient
	uint binCoef(uint n, uint k) {
		if (k > n - k)
			k = n - k;
		uint res = 1;
		for (uint i = n - k + 1; i <= n; ++i)
			res *= i;
		for (uint i = 2; i <= k; ++i)
			res /= i;
		return res;
	}

	//Returns the product of binCoef(n[i], k[i]). Assumes n is greater than or equal to k at all slots!
	uint binCoef(const mi_t& n, const mi_t& k) {
		uint result = 1;
		for (uint i = 0; i < n.size(); ++i)
			result *= binCoef(n[i], k[i]);
		return result;
	}

	//returns n!
	unsigned long long factorial(uint n) {
		unsigned long long result = 1;
		for (uint i = 2; i <= n; ++i) {
			result *= i;
		}
		return result;
	}

	//returns the product of multiIndex[i]!
	unsigned long long factorial(const mi_t& multiIndex) {
		unsigned long long result = 1;
		for (uint i = 0; i < multiIndex.size(); ++i) {
			result *= factorial(multiIndex[i]);
		}
		return result;
	}

	//returns the number of multi-indices of length n that sum to k
	uint multiIndexCardinality(uint n, uint k) {
		if (n == 0) return 0;
		if (n == 1) return 1;
		return binCoef(n + k - 1, k);
	}

	/*
	Returns an array (of length binCoef(n+k-1,n)) of arrays (of length n) that represent all multi-indices of length n that sum to k (this set is denoted I(n,k)).
	Multi-indices are contructed by iterating on length.
	n should be a positive integer and k should be a nonnegative integer; otherwise nullptr is returned
	*/
	static uint** formMultiIndices(uint n, uint k) {
		
		//cases n=0 and n=1 are handled separately

		if (n == 0)
			return nullptr;
		if (n == 1) {
			uint** res = new uint * [1];
			res[0] = new uint[1];
			res[0][0] = k;
			return res;
		}

		//the general case n>1

		uint i, j, l;
		i = 2; //i is the length of multi-indices that are currently under construction
		uint*** list_of_previous_indices = new uint * *[k + 1]; //jth element is the list of elements in I(i-1,j)
		for (j = 0; j <= k; ++j) {
			list_of_previous_indices[j] = new uint * [1];
			list_of_previous_indices[j][0] = new uint[1];
			list_of_previous_indices[j][0][0] = j;
		}
		for (i = 2; i <= n; ++i) {
			uint*** list_of_new_indices = new uint * *[k + 1]; //jth element of this will be the list of elements in I(i,j)
			for (j = 0; j <= k; ++j) {
				list_of_new_indices[j] = new uint * [multiIndexCardinality(i, j)];
				uint jj = 0;
				for (l = 0; l <= j; ++l) {
					for (uint m = 0; m < multiIndexCardinality(i - 1, j - l); ++m) {
						uint* buf = new uint[i];
						for (uint mm = 0; mm < i - 1; ++mm)
							buf[mm] = list_of_previous_indices[j - l][m][mm];
						buf[i - 1] = l;
						list_of_new_indices[j][jj++] = buf;
					}
				}
			}
			//replace previous indices with new ones and free memory
			for (j = 0; j <= k; ++j) {
				for (l = 0; l < multiIndexCardinality(i - 1, j); ++l) {
					delete[] list_of_previous_indices[j][l];
				}
				delete[] list_of_previous_indices[j];
			}
			delete[] list_of_previous_indices;
			list_of_previous_indices = list_of_new_indices;
		}

		//free memory and return result
		uint** res{ list_of_previous_indices[k] };
		for (j = 0; j < k; ++j) {
			for (l = 0; l < multiIndexCardinality(n, j); ++l) {
				delete[] list_of_previous_indices[j][l];
			}
			delete[] list_of_previous_indices[j];
		}
		delete[] list_of_previous_indices;
		return res;
	}

	//make multiIndexList contain all multi-indices of length n that sum to k
	void getMultiIndices(uint n, uint k, Buffer<mi_t>& multiIndexList) {
		uint** mi = formMultiIndices(n, k);
		uint mi_num = multiIndexCardinality(n, k);
		multiIndexList.resize(mi_num);
		for (uint i = 0; i < mi_num; ++i) {
			multiIndexList[i].resize(n);
			for (uint j = 0; j < n; ++j) {
				multiIndexList[i][j] = mi[i][j];
			}
			delete[] mi[i];
		}
		delete[] mi;
	}

	//print all multi-indices from the set I(n,k) for testing purposes
	void multiIndexTest(uint n, uint k) {
		std::cout << "Printing all multi-indices from the set I(" << n << ',' << k << ")." << std::endl;
		uint** multiIndices = formMultiIndices(n, k);
		uint size = multiIndexCardinality(n, k);
		for (uint i = 0; i < size; i++) {
			std::cout << '(';
			for (uint j = 0; j < n; j++) {
				std::cout << multiIndices[i][j];
				if (j < n - 1)
					std::cout << ',';
			}
			std::cout << ')' << std::endl;
			delete[] multiIndices[i];
		}
		delete[] multiIndices;
	}

	//put the indices of nonzero elements in multiIndex to nzIndices
	void nonzeroIndices(const mi_t& multiIndex, Buffer<uint>& nzIndices) {
		for (uint i = 0; i < multiIndex.size(); ++i)
			if (multiIndex[i] != 0)
				nzIndices.push_back(i);
	}

	//returns true if multiIndex has element at some slot, false otherwise
	bool contains(const mi_t& multiIndex, uint element) {
		for (uint i = 0; i < multiIndex.size(); ++i)
			if (multiIndex[i] == element)
				return true;
		return false;
	}

	//returns true if the multi-indices are the same
	bool equals(const mi_t& mi1, const mi_t& mi2) {
		if (mi1.size() != mi2.size())
			return false;
		for (uint i = 0; i < mi1.size(); ++i)
			if (mi1[i] != mi2[i])
				return false;
		return true;
	}

	//returns the sum of the elements of multiIndex
	uint sum(const mi_t& multiIndex) {
		return std::accumulate(multiIndex.begin(), multiIndex.end(), 0);
	}

	//returns the index of val in multiIndex
	uint indexOf(const mi_t& multiIndex, uint val) {
		for (uint i = 0; i < multiIndex.size(); ++i)
			if (multiIndex[i] == val)
				return i;
		return NONE;
	}

	//Returns the difference mi1 - mi2. Assumes mi1 is greater than or equal to mi2 at all slots!
	mi_t subtract(const mi_t& mi1, const mi_t& mi2) {
		mi_t result(mi1);
		for (uint i = 0; i < mi1.size(); ++i)
			result[i] -= mi2[i];
		return result;
	}

	//returns new multi-index with value inserted at slot index in multiIndex
	mi_t insert(const mi_t& multiIndex, uint index, uint value) {
		if (index > multiIndex.size()) {
			return multiIndex; //element not inserted
		}		
		mi_t result(multiIndex.size() + 1);
		uint i;
		for (i = 0; i < index; ++i) {
			result[i] = multiIndex[i];
		}
		result[index] = value;
		for (i = index + 1; i < result.size(); ++i) {
			result[i] = multiIndex[i-1];
		}
		return result;
	}

	//create new multi-index with additional zeros (see the implementation for details)
	mi_t addZeros(const mi_t& multiIndex, uint index1, uint index2) {
		mi_t result(multiIndex.size() + 2, 0);
		result[index1] = multiIndex[0];
		result[index2] = multiIndex[1];
		return result;
	}

	//returns the product of vec[i] ^ multiIndex[i]
	template <class T> double pow(const T& base, const mi_t& exponent) {
		double result = 1.0;
		for (uint i = 0; i < base.size(); ++i)
			result *= std::pow(base[i], exponent[i]);
		return result;
	}
	template double pow<std::vector<double>>(const std::vector<double>& base, const mi_t& exponent);
	template double pow<mi_t>(const mi_t& base, const mi_t& exponent);

	//returns new list whose element at slot i is the element of the parameter list at permutation[i]
	template <class T> T permute(const T& list, const mi_t& permutation) {
		T result(list.size());
		for (uint i = 0; i < result.size(); ++i) {
			result[i] = list[permutation[i]];
		}
		return result;
	}
	template std::vector<double> permute<std::vector<double>>(const std::vector<double>& list, const mi_t& permutation);
	template mi_t permute<mi_t>(const mi_t& list, const mi_t& permutation);

	//returns new multi-index (result) such that the parameter (multiIndex) is its permutation (multiIndex = permute(result, permutation))
	mi_t permuteOtherWay(const mi_t& multiIndex, const mi_t& permutation) {
		mi_t result(multiIndex.size());
		for (uint i = 0; i < result.size(); ++i) {
			result[permutation[i]] = multiIndex[i];
		}
		return result;
	}

	//returns the index of the permutation in edgePermutations or facePermutations such that buf = permute(vec, permutation) 
	uint getPermutation(const Buffer<uint>& buf, const std::vector<uint>& vec) {
		uint i, j;
		switch (vec.size()) {
		case 2:
			for (i = 0; i < edgePermutations.size(); ++i) {
				for (j = 0; j < vec.size(); ++j) {
					if (buf[j] != vec[edgePermutations[i][j]])
						break;
				}
				if (j == vec.size())
					break;
			}
			return i;
			break;
		case 3:
			for (i = 0; i < facePermutations.size(); ++i) {
				for (j = 0; j < vec.size(); ++j) {
					if (buf[j] != vec[facePermutations[i][j]])
						break;
				}
				if (j == vec.size())
					break;
			}
			return i;
			break;
		default:
			return 0;
			break;
		}
		
	}	

	//returns the sign of the permutation edgePermutations[permutationIndex]
	double getEdgePermutationSign(uint permutationIndex) {
		if (permutationIndex == 1)
			return -1.0;
		else return 1.0;
	}

	//returns 1 if the permutation facePermutations[permutationIndex] preserves the orientation of edge and -1 if not
	double getFacePermutationSign(uint permutationIndex, const mi_t& edge) {
		const mi_t& permutation = facePermutations[permutationIndex];
		Buffer<uint> nodes;
		nonzeroIndices(edge, nodes);
		if (permutation[nodes[0]] > permutation[nodes[1]])
			return -1.0;
		else return 1.0;
	}

	//returns the sign of the permutation facePermutations[permutationIndex]
	double getFacePermutationSign(uint permutationIndex) {
		if (permutationIndex == 0 || permutationIndex == 3 || permutationIndex == 4)
			return 1.0;
		else return -1.0;
	}

	//returns the index of multiIndex in faceEdges or bodyEdges (assuming it is in one of them)
	uint getEdgeIndex(const mi_t& multiIndex) {
		if (multiIndex.size() == 3) {
			if (multiIndex[0] == 1) {
				if (multiIndex[1] == 1)
					return 0;
				else return 1;
			}
			else return 2;
		}
		else {
			if (multiIndex[0] == 1) {
				if (multiIndex[1] == 1)
					return 0;
				else if (multiIndex[2] == 1)
					return 1;
				else return 2;
			}
			else if (multiIndex[1] == 1) {
				if (multiIndex[2] == 1)
					return 3;
				else return 4;
			}
			else return 5;
		}
	}

	//returns the index of multiIndex in bodyFaces (assuming it is there)
	uint getFaceIndex(const mi_t& multiIndex) {
		if (multiIndex[3] == 0)
			return 0;
		else if (multiIndex[2] == 0)
			return 1;
		else if (multiIndex[1] == 0)
			return 2;
		else return 3;
	}

	/*
	Fills the list multiIndices with all multi-indices whose each component is less than or equal to the corresponding component in multiIndex.
	Assumes multiIndex has length 2, 3, or 4.
	*/
	void getSmallerMultiIndices(const mi_t& multiIndex, std::vector<mi_t>& multiIndices) {
		uint i, j, k, l, m;
		if (multiIndex.size() == 2) {
			multiIndices.resize((multiIndex[0] + 1) * (multiIndex[1] + 1));
			i = 0;
			j = 0;
			k = 0;
			while (i < multiIndices.size()) {
				multiIndices[i].resize(2);
				multiIndices[i][0] = j;
				multiIndices[i][1] = k;
				++i;
				if (++j > multiIndex[0]) {
					j = 0;
					++k;
				}
			}
		}
		else if (multiIndex.size() == 3) {
			multiIndices.resize((multiIndex[0] + 1) * (multiIndex[1] + 1) * (multiIndex[2] + 1));
			i = 0;
			j = 0;
			k = 0;
			l = 0;
			while (i < multiIndices.size()) {
				multiIndices[i].resize(3);
				multiIndices[i][0] = j;
				multiIndices[i][1] = k;
				multiIndices[i][2] = l;
				++i;
				if (++j > multiIndex[0]) {
					j = 0;
					if (++k > multiIndex[1]) {
						k = 0;
						++l;
					}
				}
			}
		}
		else {
			multiIndices.resize((multiIndex[0] + 1) * (multiIndex[1] + 1) * (multiIndex[2] + 1) * (multiIndex[3] + 1));
			i = 0;
			j = 0;
			k = 0;
			l = 0;
			m = 0;
			while (i < multiIndices.size()) {
				multiIndices[i].resize(4);
				multiIndices[i][0] = j;
				multiIndices[i][1] = k;
				multiIndices[i][2] = l;
				multiIndices[i][3] = m;
				++i;
				if (++j > multiIndex[0]) {
					j = 0;
					if (++k > multiIndex[1]) {
						k = 0;
						if (++l > multiIndex[2]) {
							l = 0;
							++m;
						}
					}
				}
			}
		}
	}
}