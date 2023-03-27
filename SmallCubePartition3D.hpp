/*
SmallCubePartition3D enables the partition of a cubical mesh into small cubes in 3D.
Similar to SmallSimplexPartition3D but simpler and more incomplete.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _SMALLCUBEPARTITION3D_HPP_INCLUDED_
#define _SMALLCUBEPARTITION3D_HPP_INCLUDED_

#include "MultiIndex.hpp"
#include "GFD/Types/Types.hpp"
#include "GFD/Types/Matrix.hpp"
#include "GFD/Types/Vector.hpp"
#include "GFD/Types/Buffer.hpp"
#include "GFD/Mesh/BuilderMesh.hpp"
#include <unordered_map>

namespace gfd {

	class SmallCubePartition3D {

	public:
		SmallCubePartition3D(uint order = 1);
		void refineMesh(const BuilderMesh& mesh_old, BuilderMesh& mesh);
		void solve1FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients, Buffer<VectorN>& bodyCoefficients) const;
		Vector3 evaluate1FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<VectorN>& edgeCoefficients,
			const Buffer<VectorN>& faceCoefficients, const Buffer<VectorN>& bodyCoefficients, uint element) const;
		void formHodgeMatrix1Forms(Buffer<std::unordered_map<uint, double>>& star, bool minkowskiMetric = false) const;
		uint smallNodesEdgeList(uint i, uint j) const;
		uint smallEdgesEdgeList(uint i, uint j) const;
		uint smallNodesFaceList(uint i, uint j) const;
		uint smallEdgesFaceList(uint i, uint j) const;
		uint smallFacesFaceList(uint i, uint j) const;
		uint smallNodesBodyList(uint i, uint j) const;
		uint smallEdgesBodyList(uint i, uint j) const;
		uint smallFacesBodyList(uint i, uint j) const;
		uint smallBodiesBodyList(uint i, uint j) const;

	protected:
		const uint order;
		const BuilderMesh* mesh_old_ptr;
		const BuilderMesh* mesh_ptr;

		struct SmallEdge {
			mi_t mi;
			uint face;
		};

		void initialiseSmallCells();
		uint getSmallNodeIndex(const Vector3& pos, uint element, uint dim);
		Vector3 evaluate1FormWithEdgeCoefficients(const Vector3& evaluationPoint, const Buffer<Vector3>& pos, const VectorN& edgeCoefficients, uint edge) const;
		Vector3 evaluate1FormWithFaceCoefficients(const Vector3& evaluationPoint, const Buffer<Vector3>& pos, const VectorN& faceCoefficients, uint face) const;
		Vector3 evaluate1FormWithBodyCoefficients(const Vector3& evaluationPoint, const Buffer<Vector3>& pos, const VectorN& bodyCoefficients) const;
		void formMatrices();

		Buffer<Buffer<double>> integrals1Forms;
		MatrixN matrix1FormsEdges;
		Buffer<uint> matrix1FormsEdges_p;
		MatrixN matrix1FormsFaces;
		Buffer<uint> matrix1FormsFaces_p;
		MatrixN matrix1FormsBodies;
		Buffer<uint> matrix1FormsBodies_p;
		std::vector<SmallEdge> smallEdgesInEdges;
		std::vector<SmallEdge> smallEdgesInFaces;
		std::vector<SmallEdge> smallEdgesInBodies;
		Buffer<Buffer<VectorN>> faceEdgeValues1Forms;
		Buffer<Buffer<VectorN>> bodyEdgeValues1Forms;
		Buffer<Buffer<VectorN>> bodyFaceValues1Forms;
	};
}

#endif //_SMALLCUBEPARTITION3D_HPP_INCLUDED_