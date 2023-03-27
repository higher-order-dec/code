/*
SmallCubePartition2D enables the partition of a cubical mesh into small cubes in 2D.
Similar to SmallSimplexPartition2D but simpler and more incomplete.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _SMALLCUBEPARTITION2D_HPP_INCLUDED_
#define _SMALLCUBEPARTITION2D_HPP_INCLUDED_

#include "MultiIndex.hpp"
#include "GFD/Types/Types.hpp"
#include "GFD/Types/Matrix.hpp"
#include "GFD/Types/Vector.hpp"
#include "GFD/Types/Buffer.hpp"
#include "GFD/Mesh/BuilderMesh.hpp"
#include <unordered_map>

namespace gfd {

	class SmallCubePartition2D {

	public:
		SmallCubePartition2D(uint order = 1);
		void refineMesh(const BuilderMesh& mesh_old, BuilderMesh& mesh);
		void solve1FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients) const;
		Vector2 evaluate1FormWithCoefficients(const Vector2& evaluationPoint, const Buffer<VectorN>& edgeCoefficients,
			const Buffer<VectorN>& faceCoefficients, uint element) const;
		void formHodgeMatrix1Forms(Buffer<std::unordered_map<uint, double>>& star, bool minkowskiMetric = false) const;
		uint smallNodesEdgeList(uint i, uint j) const;
		uint smallEdgesEdgeList(uint i, uint j) const;
		uint smallNodesFaceList(uint i, uint j) const;
		uint smallEdgesFaceList(uint i, uint j) const;
		uint smallFacesFaceList(uint i, uint j) const;

	protected:
		const uint order;
		const BuilderMesh* mesh_old_ptr;
		const BuilderMesh* mesh_ptr;

		struct SmallEdge {
			mi_t mi;
			uint edge;
		};

		uint getSmallNodeIndex(const Vector2& pos, uint element, uint dim);
		void initialiseSmallCells();
		Vector2 evaluate1FormWithEdgeCoefficients(const Vector2& evaluationPoint, const Buffer<Vector2>& pos, const VectorN& edgeCoefficients, uint edge) const;
		Vector2 evaluate1FormWithFaceCoefficients(const Vector2& evaluationPoint, const Buffer<Vector2>& pos, const VectorN& faceCoefficients) const;
		void formMatrices();

		Buffer<Buffer<double>> integrals1Forms;
		MatrixN matrix1FormsEdges;
		Buffer<uint> matrix1FormsEdges_p;
		MatrixN matrix1FormsFaces;
		Buffer<uint> matrix1FormsFaces_p;
		std::vector<SmallEdge> smallEdgesInEdges;
		std::vector<SmallEdge> smallEdgesInFaces;
		Buffer<Buffer<VectorN>> faceEdgeValues1Forms;
	};
}

#endif //_SMALLCUBEPARTITION2D_HPP_INCLUDED_
