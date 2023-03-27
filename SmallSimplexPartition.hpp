/*
SmallSimplexPartition is a base class for the partition of a simplicial mesh into small simplices.
The base class contains features that are required in both SmallSimplexPartition2D and SmallSimplexPartition3D.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _SMALLSIMPLEXPARTITION_HPP_INCLUDED_
#define _SMALLSIMPLEXPARTITION_HPP_INCLUDED_

#include "MultiIndex.hpp"
#include "GFD/Types/Types.hpp"
#include "GFD/Types/Vector.hpp"
#include "GFD/Types/Buffer.hpp"
#include "GFD/Mesh/BuilderMesh.hpp"
#include <functional>

namespace gfd {

	class SmallSimplexPartition {

	public:
		static std::function<double(Vector2)> getBarycentricProductRefSimplex2D(const mi_t& multiIndex);
		static std::function<double(Vector3)> getBarycentricProductRefSimplex3D(const mi_t& multiIndex);
		static double integrateWhitneyForm(const mi_t& multiIndex, const mi_t& nodeIndices, const mi_t& whitneyForm);
		static double integrateWhitneyForm(const mi_t& miSimplex, const mi_t& niSimplex, const mi_t& miIntegrand, const mi_t& niIntegrand);
		static double integrateWhitneyFormNumerically(const mi_t& miSimplex, const mi_t& niSimplex, const mi_t& miIntegrand, const mi_t& niIntegrand);
		static double integralAverage(const mi_t& smallSimplex, const mi_t& barycentricProduct);
		static double integralAverage(const mi_t& multiIndex, const mi_t& nodeIndices, const mi_t& barycentricProduct);
		static double integralAverageNumerically(const mi_t& multiIndex, const mi_t& nodeIndices, const mi_t& barycentricProduct);
		static bool mapsToInterior(const mi_t& multiIndex, const mi_t& subsimplex);
		static Vector2 getNodeImageRefSimplex2D(const mi_t& multiIndex, uint node, uint sum = 0);
		static Vector3 getNodeImageRefSimplex3D(const mi_t& multiIndex, uint node, uint sum = 0);

	protected:
		/*
		Small simplices are labelled with respect to a big simplex by giving a multi-index and the indices of subsimplex whose image the small simplex is
		multiIndex is given with respect to the big simplex
		nodeIndices has 1 in the node slots of the subsimplex whose image this small simplex is, 0 in other slots
		Examples:
		(multiIndex; nodeIndices) = (1,2,0; 1,0,1) in triangle ijk is the small edge that is the image of the edge ik through the multi-index (1,2,0)
		(multiIndex; nodeIndices) = (2,0,0,2; 0,0,1,0) in tetrahedron ijkl is the small node that is the image of the node k through the multi-index (2,0,0,2)
		(multiIndex; nodeIndices) = (2,0,0,2; 0,1,1,1) in tetrahedron ijkl is small face that is the image of the face jkl through the multi-index (2,0,0,2)
		*/
		struct SmallSimplex {
			const mi_t* multiIndex;
			const mi_t* nodeIndices;
		};

		SmallSimplexPartition(uint ord);
		Vector4 getNodeImage(const mi_t& multiIndex, const Mesh& mesh, const Buffer<uint>& nodes, const mi_t& node) const;
		static uint indexOf(const Buffer<mi_t>& multiIndexList, const mi_t& multiIndex);
		static uint indexOf(const Buffer<SmallSimplex>& smallSimplexList, const mi_t& multiIndex);
		static uint indexOf(const Buffer<SmallSimplex>& smallSimplexList, const SmallSimplex& smallSimplex);
		static void printSmallSimplexList(const Buffer<SmallSimplex>& ssList);

		const uint order;
		const BuilderMesh* mesh_old_ptr;
		const BuilderMesh* mesh_ptr;

	};
}

#endif //_SMALLSIMPLEXPARTITION_HPP_INCLUDED_