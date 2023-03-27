/*
SmallSimplexPartition3D extends the base class SmallSimplexPartition and enables the partition of a simplicial mesh into small simplices in 3D.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _SMALLSIMPLEXPARTITION3D_HPP_INCLUDED_
#define _SMALLSIMPLEXPARTITION3D_HPP_INCLUDED_

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

	class SmallSimplexPartition3D : public SmallSimplexPartition {

	public:
		SmallSimplexPartition3D(uint order = 1, bool loadMatricesFromFile = false);
		void refineMesh(const BuilderMesh& mesh_old, BuilderMesh& mesh);
		double evaluate0Form(const Buffer<double>& discreteForm, const Vector3& evaluationPoint) const;
		Vector3 evaluate1Form(const Buffer<double>& discreteForm, const Vector3& evaluationPoint) const;
		Vector3 evaluate2Form(const Buffer<double>& discreteForm, const Vector3& evaluationPoint) const;
		double evaluate3Form(const Buffer<double>& discreteForm, const Vector3& evaluationPoint) const;
		void solve0FormCoefficients(const Buffer<double>& discreteForm, Buffer<double>& nodeCoefficients, Buffer<VectorN>& edgeCoefficients,
			Buffer<VectorN>& faceCoefficients, Buffer<VectorN>& bodyCoefficients) const;
		void solve1FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients, Buffer<VectorN>& bodyCoefficients) const;
		void solve2FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& faceCoefficients, Buffer<VectorN>& bodyCoefficients) const;
		void solve3FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& bodyCoefficients) const;
		double evaluate0FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<double>& nodeCoefficients, const Buffer<VectorN>& edgeCoefficients,
			const Buffer<VectorN>& faceCoefficients, const Buffer<VectorN>& bodyCoefficients, uint element = 0) const;
		Vector3 evaluate1FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<VectorN>& edgeCoefficients,
			const Buffer<VectorN>& faceCoefficients, const Buffer<VectorN>& bodyCoefficients, uint element = 0) const;
		Vector3 evaluate2FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<VectorN>& faceCoefficients, const Buffer<VectorN>& bodyCoefficients, uint element = 0) const;
		double evaluate3FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<VectorN>& bodyCoefficients, uint element = 0) const;
		uint estimateNonzeros() const;
		void formHodgeMatrix1Forms(Buffer<std::unordered_map<uint, double>>& star, bool circumcentric = false, bool minkowskiMetric = false) const;
		void formHodgeMatrix1FormsCustom(Buffer<std::unordered_map<uint, double>>& star, std::string foldername, bool circumcentric = false, bool loadIntegralsFromFile = false) const;
		void formHodgeMatrix1FormsCustom(Buffer<std::unordered_map<uint, double>>& star, const Buffer<uint>& refElementIndices, const Buffer<Buffer<uint>>& refElements,
			std::string foldername, bool circumcentric = false, bool minkowskiMetric = false, bool loadIntegralsFromFile = false) const;
		uint smallNodesEdgeList(uint i, uint j) const;
		uint smallEdgesEdgeList(uint i, uint j) const;
		uint smallNodesFaceList(uint i, uint j) const;
		uint smallEdgesFaceList(uint i, uint j) const;
		uint smallFacesFaceList(uint i, uint j) const;
		uint smallNodesBodyList(uint i, uint j) const;
		uint smallEdgesBodyList(uint i, uint j) const;
		uint smallFacesBodyList(uint i, uint j) const;
		uint smallBodiesBodyList(uint i, uint j) const;

	private:
		void initialiseMultiIndices();
		void formActiveSmallSimplices();
		Buffer<uint> getBodyNodesCustom(uint b) const;
		uint firstEdgeOfFace(uint i) const;
		uint firstEdgeOfBody(uint i) const;
		uint firstFaceOfBody(uint i) const;
		mi_t* getHoleEdgeMultiIndex(const mi_t& holeIndex, uint edgeIndex) const;
		mi_t* getHoleFaceMultiIndex(const mi_t& holeIndex, uint faceIndex) const;
		uint getOppositeEdge(uint edge) const;
		uint getSmallNodeIndex(mi_t multiIndex, uint nodeIndex, Buffer<uint> bigSimplexNodes) const;
		Vector3 evaluate1FormWithEdgeCoefficients(const Vector3& evaluationPoint, const VectorN& edgeCoefficients, uint edgeIndex, uint edgePermutation,
			const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
			const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const;
		Vector3 evaluate1FormWithCoefficients(const Vector3& evaluationPoint, const VectorN& edgeCoefficients, uint edgeIndex, uint edgePermutation, const VectorN& face0Coefficients,
			uint face0Index, uint face0Permutation, const VectorN& face1Coefficients, uint face1Index, uint face1Permutation, const VectorN& bodyCoefficients,
			const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
			const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const;
		Vector3 evaluate1FormWithFaceCoefficients(const Vector3& evaluationPoint, const VectorN& faceCoefficients, uint faceIndex, uint facePermutation,
			const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
			const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const;
		Vector3 evaluate1FormWithCoefficients(const Vector3& evaluationPoint, const VectorN& faceCoefficients, uint faceIndex, uint facePermutation, const VectorN& bodyCoefficients,
			const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
			const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const;
		Vector3 evaluate1FormWithBodyCoefficients(const Vector3& evaluationPoint, const VectorN& bodyCoefficients,
			const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
			const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const;
		void formMatrices();
		void saveMatrices() const;
		void loadMatrices();
		void saveHodgeIntegrals1Forms(std::string foldername, bool circumcentric, double scalingFactor, const Buffer<VectorN>& hodgeIntegralsBodyEdges,
			const Buffer<Buffer<VectorN>>& hodgeIntegralsFace0Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsFace1Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsFace2Edges,
			const Buffer<Buffer<VectorN>>& hodgeIntegralsFace3Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge0Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge1Edges,
			const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge2Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge3Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge4Edges, 
			const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge5Edges, int ref = -1) const;
		double loadHodgeIntegrals1Forms(std::string foldername, bool circumcentric, Buffer<VectorN>& hodgeIntegralsBodyEdges, Buffer<Buffer<VectorN>>& hodgeIntegralsFace0Edges,
			Buffer<Buffer<VectorN>>& hodgeIntegralsFace1Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsFace2Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsFace3Edges,
			Buffer<Buffer<VectorN>>& hodgeIntegralsEdge0Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge1Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge2Edges,
			Buffer<Buffer<VectorN>>& hodgeIntegralsEdge3Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge4Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge5Edges, int ref = -1) const;

		Buffer<mi_t> edgeMultiIndices;
		Buffer<mi_t> faceMultiIndices;
		Buffer<mi_t> bodyMultiIndices;
		Buffer<mi_t> faceEdgeHoles;
		Buffer<mi_t> bodyEdgeHoles;
		Buffer<mi_t> bodyFaceHoles;
		Buffer<SmallSimplex> activeSmallNodesInEdges;
		Buffer<SmallSimplex> activeSmallEdgesInEdges;
		Buffer<SmallSimplex> activeSmallNodesInFaces;
		Buffer<SmallSimplex> activeSmallEdgesInFaces;
		Buffer<SmallSimplex> activeSmallFacesInFaces;
		Buffer<SmallSimplex> activeSmallNodesInBodies;
		Buffer<SmallSimplex> activeSmallEdgesInBodies;
		Buffer<SmallSimplex> activeSmallFacesInBodies;
		Buffer<SmallSimplex> activeSmallBodiesInBodies;
		Buffer<uint> smallEdgesFaceOffsets;
		Buffer<uint> smallEdgesBodyOffsets;
		Buffer<uint> smallFacesBodyOffsets;
		MatrixN matrix0FormsEdges;
		Buffer<uint> matrix0FormsEdges_p;
		MatrixN matrix0FormsFaces;
		Buffer<uint> matrix0FormsFaces_p;
		MatrixN matrix0FormsBodies;
		Buffer<uint> matrix0FormsBodies_p;
		MatrixN matrix1FormsEdges;
		Buffer<uint> matrix1FormsEdges_p;
		MatrixN matrix1FormsFaces;
		Buffer<uint> matrix1FormsFaces_p;
		MatrixN matrix1FormsBodies;
		Buffer<uint> matrix1FormsBodies_p;
		MatrixN matrix2FormsFaces;
		Buffer<uint> matrix2FormsFaces_p;
		MatrixN matrix2FormsBodies;
		Buffer<uint> matrix2FormsBodies_p;
		MatrixN matrix3Forms;
		Buffer<uint> matrix3Forms_p;
		VectorN edgeValuesNode0;
		VectorN edgeValuesNode1;
		VectorN faceValuesNode0;
		VectorN faceValuesNode1;
		VectorN faceValuesNode2;
		VectorN bodyValuesNode0;
		VectorN bodyValuesNode1;
		VectorN bodyValuesNode2;
		VectorN bodyValuesNode3;
		Buffer<Buffer<VectorN>> faceNodeValuesEdge0;
		Buffer<Buffer<VectorN>> faceNodeValuesEdge1;
		Buffer<Buffer<VectorN>> faceNodeValuesEdge2;
		Buffer<Buffer<VectorN>> bodyNodeValuesEdge0;
		Buffer<Buffer<VectorN>> bodyNodeValuesEdge1;
		Buffer<Buffer<VectorN>> bodyNodeValuesEdge2;
		Buffer<Buffer<VectorN>> bodyNodeValuesEdge3;
		Buffer<Buffer<VectorN>> bodyNodeValuesEdge4;
		Buffer<Buffer<VectorN>> bodyNodeValuesEdge5;
		Buffer<Buffer<VectorN>> bodyNodeValuesFace0;
		Buffer<Buffer<VectorN>> bodyNodeValuesFace1;
		Buffer<Buffer<VectorN>> bodyNodeValuesFace2;
		Buffer<Buffer<VectorN>> bodyNodeValuesFace3;
		Buffer<Buffer<VectorN>> faceEdgeValuesEdge0;
		Buffer<Buffer<VectorN>> faceEdgeValuesEdge1;
		Buffer<Buffer<VectorN>> faceEdgeValuesEdge2;
		Buffer<Buffer<VectorN>> bodyEdgeValuesEdge0;
		Buffer<Buffer<VectorN>> bodyEdgeValuesEdge1;
		Buffer<Buffer<VectorN>> bodyEdgeValuesEdge2;
		Buffer<Buffer<VectorN>> bodyEdgeValuesEdge3;
		Buffer<Buffer<VectorN>> bodyEdgeValuesEdge4;
		Buffer<Buffer<VectorN>> bodyEdgeValuesEdge5;
		Buffer<Buffer<VectorN>> bodyEdgeValuesFace0;
		Buffer<Buffer<VectorN>> bodyEdgeValuesFace1;
		Buffer<Buffer<VectorN>> bodyEdgeValuesFace2;
		Buffer<Buffer<VectorN>> bodyEdgeValuesFace3;
		Buffer<Buffer<VectorN>> bodyFaceValuesFace0;
		Buffer<Buffer<VectorN>> bodyFaceValuesFace1;
		Buffer<Buffer<VectorN>> bodyFaceValuesFace2;
		Buffer<Buffer<VectorN>> bodyFaceValuesFace3;
	};
}

#endif //_SMALLSIMPLEXPARTITION3D_HPP_INCLUDED_
