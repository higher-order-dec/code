/*
SmallSimplexPartition2D extends the base class SmallSimplexPartition and enables the partition of a simplicial mesh into small simplices in 2D.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _SMALLSIMPLEXPARTITION2D_HPP_INCLUDED_
#define _SMALLSIMPLEXPARTITION2D_HPP_INCLUDED_

#include "SmallSimplexPartition.hpp"
#include "MultiIndex.hpp"
#include "GFD/Types/Types.hpp"
#include "GFD/Types/Matrix.hpp"
#include "GFD/Types/Vector.hpp"
#include "GFD/Types/Buffer.hpp"
#include "GFD/Mesh/BuilderMesh.hpp"
#include <unordered_map>
#include <vector>
#include <iostream>
#include <functional>

namespace gfd {

	class SmallSimplexPartition2D : public SmallSimplexPartition {

	public:
		SmallSimplexPartition2D(uint order = 1, bool loadMatricesFromFile = false);
		void refineMesh(const BuilderMesh& mesh_old, BuilderMesh& mesh);
		double evaluate0Form(const Buffer<double>& discreteForm, const Vector2& evaluationPoint) const;
		Vector2 evaluate1Form(const Buffer<double>& discreteForm, const Vector2& evaluationPoint) const;
		double evaluate2Form(const Buffer<double>& discreteForm, const Vector2& evaluationPoint) const;
		void solve0FormCoefficients(const Buffer<double>& discreteForm, Buffer<double>& nodeCoefficients, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients) const;
		void solve1FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients) const;
		void solve2FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& faceCoefficients) const;
		double evaluate0FormWithCoefficients(const Vector2& evaluationPoint, const Buffer<double>& nodeCoefficients, const Buffer<VectorN>& edgeCoefficients,
			const Buffer<VectorN>& faceCoefficients, uint element = 0) const;
		Vector2 evaluate1FormWithCoefficients(const Vector2& evaluationPoint, const Buffer<VectorN>& edgeCoefficients, const Buffer<VectorN>& faceCoefficients, uint element = 0) const;
		double evaluate2FormWithCoefficients(const Vector2& evaluationPoint, const Buffer<VectorN>& faceCoefficients, uint element = 0) const;
		void formHodgeMatrix1Forms(Buffer<std::unordered_map<uint, double>>& star, bool circumcentric = false, bool minkowskiMetric = false) const;
		void formHodgeMatrix1FormsCustom(Buffer<std::unordered_map<uint, double>>& star, std::string foldername, bool circumcentric = false, bool loadIntegralsFromFile = false) const;
		void formHodgeMatrix1FormsCustom(Buffer<std::unordered_map<uint, double>>& star, const Buffer<uint>& refElementIndices, const Buffer<Buffer<uint>>& refElements,
			std::string foldername, bool circumcentric = false, bool minkowskiMetric = false, bool loadIntegralsFromFile = false) const;
		uint smallNodesEdgeList(uint i, uint j) const;
		uint smallEdgesEdgeList(uint i, uint j) const;
		uint smallNodesFaceList(uint i, uint j) const;
		uint smallEdgesFaceList(uint i, uint j) const;
		uint smallFacesFaceList(uint i, uint j) const;

	private:
		void initialiseMultiIndices();
		void formActiveSmallSimplices();
		Buffer<uint> getFaceNodesCustom(uint f) const;
		uint firstEdgeOfFace(uint i) const;
		mi_t* getHoleEdgeMultiIndex(const mi_t& holeIndex, uint edgeIndex) const;
		uint getSmallNodeIndex(mi_t multiIndex, uint nodeIndex, Buffer<uint> bigSimplexNodes) const;
		Vector2 evaluate1FormWithEdgeCoefficients(const Vector2& evaluationPoint, const VectorN& edgeCoefficients, uint edgeIndex, uint edgePermutation,
			const Buffer<uint>& nodes, const Buffer<uint>& edges, const Vector2& gradLambda0, const Vector2& gradLambda1, const Vector2& gradLambda2) const;
		Vector2 evaluate1FormWithCoefficients(const Vector2& evaluationPoint, const VectorN& edgeCoefficients, uint edgeIndex, uint edgePermutation, const VectorN& faceCoefficients,
			uint element, const Buffer<uint>& nodes, const Buffer<uint>& edges, const Vector2& gradLambda0, const Vector2& gradLambda1, const Vector2& gradLambda2) const;
		Vector2 evaluate1FormWithFaceCoefficients(const Vector2& evaluationPoint, const VectorN& faceCoefficients,
			const Buffer<uint>& nodes, const Buffer<uint>& edges, const Vector2& gradLambda0, const Vector2& gradLambda1, const Vector2& gradLambda2) const;
		void formMatrices();
		void saveMatrices() const;
		void loadMatrices();
		void saveHodgeIntegrals1Forms(std::string foldername, bool circumcentric, const Buffer<VectorN>& hodgeIntegralsFaceEdges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge0Edges,
			const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge1Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge2Edges, int ref = -1) const;
		void loadHodgeIntegrals1Forms(std::string foldername, bool circumcentric, Buffer<VectorN>& hodgeIntegralsFaceEdges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge0Edges,
			Buffer<Buffer<VectorN>>& hodgeIntegralsEdge1Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge2Edges, int ref = -1) const;

		Buffer<mi_t> edgeMultiIndices;
		Buffer<mi_t> faceMultiIndices;
		Buffer<mi_t> faceEdgeHoles;
		Buffer<SmallSimplex> activeSmallNodesInEdges;
		Buffer<SmallSimplex> activeSmallEdgesInEdges;
		Buffer<SmallSimplex> activeSmallNodesInFaces;
		Buffer<SmallSimplex> activeSmallEdgesInFaces;
		Buffer<SmallSimplex> activeSmallFacesInFaces;
		Buffer<uint> smallEdgesFaceOffsets;
		MatrixN matrix0FormsEdges;
		Buffer<uint> matrix0FormsEdges_p;
		MatrixN matrix0FormsFaces;
		Buffer<uint> matrix0FormsFaces_p;
		MatrixN matrix1FormsEdges;
		Buffer<uint> matrix1FormsEdges_p;
		MatrixN matrix1FormsFaces;
		Buffer<uint> matrix1FormsFaces_p;
		MatrixN matrix2FormsFaces;
		Buffer<uint> matrix2FormsFaces_p;
		VectorN edgeValuesNode0;
		VectorN edgeValuesNode1;
		VectorN faceValuesNode0;
		VectorN faceValuesNode1;
		VectorN faceValuesNode2;
		Buffer<Buffer<VectorN>> faceNodeValuesEdge0;
		Buffer<Buffer<VectorN>> faceNodeValuesEdge1;
		Buffer<Buffer<VectorN>> faceNodeValuesEdge2;
		Buffer<Buffer<VectorN>> faceEdgeValuesEdge0;
		Buffer<Buffer<VectorN>> faceEdgeValuesEdge1;
		Buffer<Buffer<VectorN>> faceEdgeValuesEdge2;
	};
}

#endif //_SMALLSIMPLEXPARTITION2D_HPP_INCLUDED_
