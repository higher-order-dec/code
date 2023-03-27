#include "SmallSimplexPartition2D.hpp"
#include "WhitneyForm.hpp"
#include "NumericalIntegration.hpp"
#include "LinearAlgebraFunctions.hpp"

using namespace gfd;

SmallSimplexPartition2D::SmallSimplexPartition2D(uint order, bool loadMatricesFromFile) : SmallSimplexPartition(order) {
	initialiseMultiIndices();
	formActiveSmallSimplices();
	if (loadMatricesFromFile) {
		loadMatrices();
	}
	else {
		formMatrices();
		saveMatrices();
	}
}

//initialise relevant multi-indices
void SmallSimplexPartition2D::initialiseMultiIndices() {
	//initialise the multi-indices in edges and faces
	getMultiIndices(2, order - 1, edgeMultiIndices);
	getMultiIndices(3, order - 1, faceMultiIndices);
	//initialise the hole indices in faces
	if (order >= 2) {
		getMultiIndices(3, order - 2, faceEdgeHoles);
	}
}

//index the small simplices that are in the interior of a big simplex, choosing and omitting redundant ones
void SmallSimplexPartition2D::formActiveSmallSimplices() {

	uint i, j;

	//edge as big simplex
	activeSmallEdgesInEdges.resize(edgeMultiIndices.size());
	for (i = 0; i < edgeMultiIndices.size(); ++i) { //small nodes and small edges
		if (i > 0) {
			activeSmallNodesInEdges.push_back({ &edgeMultiIndices[i], &edgeNodes[0] });
		}
		activeSmallEdgesInEdges[i].multiIndex = &edgeMultiIndices[i];
		activeSmallEdgesInEdges[i].nodeIndices = &edgeEdge;
	}

	//face as big simplex
	activeSmallFacesInFaces.resize(faceMultiIndices.size());
	for (i = 0; i < faceMultiIndices.size(); ++i) { //small nodes and small faces
		if (mapsToInterior(faceMultiIndices[i], faceNodes[0])) {
			activeSmallNodesInFaces.push_back({ &faceMultiIndices[i], &faceNodes[0] });
		}
		activeSmallFacesInFaces[i].multiIndex = &faceMultiIndices[i];
		activeSmallFacesInFaces[i].nodeIndices = &faceFace;
	}
	if (order >= 2) { //small edges
		for (i = 0; i < faceEdgeHoles.size(); ++i) {
			for (j = 0; j < faceEdges.size(); ++j) {
				if (j == 0) { // omit small edges in the direction of edge 0
					continue;
				}
				//if (i % faceEdges.size() == j) { // alternative: omit a similar number of edges in each direction
				//	continue;
				//}
				activeSmallEdgesInFaces.push_back({ getHoleEdgeMultiIndex(faceEdgeHoles[i], j), &faceEdges[j] });
			}
		}
	}
}

//refine mesh and save the required information that implies indices of small simplices
void SmallSimplexPartition2D::refineMesh(const BuilderMesh& mesh_old, BuilderMesh& mesh) {

	uint i, j, k;
	//SmallSimplexPartition2D is now associated with the old mesh (mesh_old) and its refinement (mesh)
	mesh_old_ptr = &mesh_old;
	mesh_ptr = &mesh;

	//in the lowest order case, the mesh remains the same
	if (order == 1) {
		mesh.createCopy(mesh_old);
		return;
	}

	//first add the existing nodes
	for (i = 0; i < mesh_old.getNodeSize(); ++i) {
		mesh.addNode(mesh_old.getNodePosition(i));
	}

	//refine edges
	for (i = 0; i < mesh_old.getEdgeSize(); ++i) {
		const Buffer<uint>& nodes = mesh_old.getEdgeNodes(i);
		//small nodes in edge i
		for (j = 0; j < activeSmallNodesInEdges.size(); ++j) {
			mesh.addNode(getNodeImage(*activeSmallNodesInEdges[j].multiIndex, mesh_old, nodes, *activeSmallNodesInEdges[j].nodeIndices));
		}
		//small edges in edge i
		mesh.addEdge(nodes[0], smallNodesEdgeList(i, 0));
		for (j = 1; j < activeSmallEdgesInEdges.size() - 1; ++j) {
			mesh.addEdge(smallNodesEdgeList(i, j - 1), smallNodesEdgeList(i, j));
		}
		mesh.addEdge(smallNodesEdgeList(i, j - 1), nodes[1]);
	}

	//refine faces
	smallEdgesFaceOffsets.resize(activeSmallEdgesInFaces.size());
	for (i = 0; i < mesh_old.getFaceSize(); ++i) {
		const Buffer<uint> nodes = getFaceNodesCustom(i);
		//small nodes in face i
		for (j = 0; j < activeSmallNodesInFaces.size(); ++j) {
			mesh.addNode(getNodeImage(*activeSmallNodesInFaces[j].multiIndex, mesh_old, nodes, *activeSmallNodesInFaces[j].nodeIndices));
		}
		//small edges and small faces in face i
		for (j = 0; j < faceMultiIndices.size(); ++j) {
			uint node0 = getSmallNodeIndex(faceMultiIndices[j], 0, nodes);
			uint node1 = getSmallNodeIndex(faceMultiIndices[j], 1, nodes);
			uint node2 = getSmallNodeIndex(faceMultiIndices[j], 2, nodes);
			Buffer<uint> edges(3);
			edges[0] = mesh.addEdge(node0, node1);
			edges[1] = mesh.addEdge(node0, node2);
			uint edge2 = edges[2] = mesh.addEdge(node1, node2);
			if (i == 0) { //the information implying the indices of small simplices is collected in face 0
				for (k = 0; k < edges.size(); ++k) {
					if (mapsToInterior(faceMultiIndices[j], faceEdges[k])) {
						SmallSimplex smallEdge{ &faceMultiIndices[j], &faceEdges[k] };
						uint smallEdgeIndex = indexOf(activeSmallEdgesInFaces, smallEdge);
						if (smallEdgeIndex != NONE)
							smallEdgesFaceOffsets[smallEdgeIndex] = edges[k];
					}
				}
			}
			edges[2] = edges[1]; //swap these to orient the face as expected
			edges[1] = edge2;
			mesh.addFace(edges);
		}
		//fill the holes in face i
		for (j = 0; j < faceEdgeHoles.size(); ++j) {
			Buffer<uint> holeNodes(3);
			for (k = 0; k < holeNodes.size(); ++k) {
				mi_t multiIndex{ faceEdgeHoles[j] };
				multiIndex[k] += 1;
				holeNodes[k] = getSmallNodeIndex(multiIndex, (k + 1) % holeNodes.size(), nodes);
			}
			mesh.addFace(findEdges(holeNodes[0], holeNodes[1], holeNodes[2], mesh));
		}
	}
}

/*
Computes the value of the (higher order) Whitney 0-form (interpolant of discreteForm) at evaluationPoint.
When evaluating many times, it is more efficient to first solve the coefficients with solve0FormCoefficients and then evaluate with evaluate0FormWithCoefficients.
*/
double SmallSimplexPartition2D::evaluate0Form(const Buffer<double>& discreteForm, const Vector2& evaluationPoint) const {

	if (!mesh_old_ptr)
		return 0.0; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4

	//find the element containing evaluationPoint
	static uint element = 0;
	if (!mesh_old.findElement(evaluationPoint4, element))
		return 0.0; //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getFaceNodesCustom(element);
	Buffer<double> nodeCoefficients(nodes.size()); //coefficients for small nodes that are also big nodes are obtained at once
	uint i, j;
	for (i = 0; i < nodeCoefficients.size(); ++i) {
		nodeCoefficients[i] = discreteForm[nodes[i]];
	}

	Buffer<uint> edges(3);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[1], nodes[2]);
	Buffer<VectorN> edgeCoefficients(edges.size()); //coefficients for small nodes that are on big edges
	for (i = 0; i < edgeCoefficients.size(); ++i) {
		Buffer<uint> edgeNodes = mesh_old.getEdgeNodes(edges[i]);
		VectorN edgeValues = edgeValuesNode0 * discreteForm[edgeNodes[0]] + edgeValuesNode1 * discreteForm[edgeNodes[1]];
		for (j = 0; j < activeSmallNodesInEdges.size(); ++j) {
			edgeValues[j] += discreteForm[smallNodesEdgeList(edges[i], j)];
		}
		solveLUP(matrix0FormsEdges, matrix0FormsEdges_p, edgeValues, edgeCoefficients[i]);
	}

	VectorN faceCoefficients; //coefficients for small nodes that are on big faces
	VectorN faceValues = faceValuesNode0 * discreteForm[nodes[0]] + faceValuesNode1 * discreteForm[nodes[1]] + faceValuesNode2 * discreteForm[nodes[2]];
	uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
	uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
	uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[1], nodes[2] });
	for (j = 0; j < activeSmallNodesInFaces.size(); ++j) {
		faceValues[j] += faceNodeValuesEdge0[permutationEdge0][j].dot(edgeCoefficients[edges.findFirst(edges[0])]);
		faceValues[j] += faceNodeValuesEdge1[permutationEdge1][j].dot(edgeCoefficients[edges.findFirst(edges[1])]);
		faceValues[j] += faceNodeValuesEdge2[permutationEdge2][j].dot(edgeCoefficients[edges.findFirst(edges[2])]);
		faceValues[j] += discreteForm[smallNodesFaceList(element, j)];
	}
	solveLUP(matrix0FormsFaces, matrix0FormsFaces_p, faceValues, faceCoefficients);

	//find barycentric coordinates of evaluationPoint
	double h0 = mesh_old.getEdgeDeviation(edges[2], mesh_old.getNodePosition(nodes[0])).len();
	double lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() / h0;
	double h1 = mesh_old.getEdgeDeviation(edges[1], mesh_old.getNodePosition(nodes[1])).len();
	double lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() / h1;
	double h2 = mesh_old.getEdgeDeviation(edges[0], mesh_old.getNodePosition(nodes[2])).len();
	double lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() / h2;

	//compute the result using the coefficients solved above
	std::vector<double> barycentricCoords;

	//contribution from small nodes that are also big nodes
	double result = std::pow(lambda0, order) * nodeCoefficients[0] + std::pow(lambda1, order) * nodeCoefficients[1]
		+ std::pow(lambda2, order) * nodeCoefficients[2];

	//small nodes on edges
	//edge0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1 }, edgePermutations[permutationEdge0]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[0][i];
	}
	//edge1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2 }, edgePermutations[permutationEdge1]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[1][i];
	}
	//edge2
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[permutationEdge2]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[2][i];
	}

	//small nodes in the interior
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(lambda0, mi[0] + 1) * std::pow(lambda1, mi[1]) * std::pow(lambda2, mi[2]) * faceCoefficients[i];
	}

	return result;
}

/*
Computes the value of the (higher order) Whitney 1-form (interpolant of discreteForm) at evaluationPoint.
When evaluating many times, it is more efficient to first solve the coefficients with solve1FormCoefficients and then evaluate with evaluate1FormWithCoefficients.
*/
Vector2 SmallSimplexPartition2D::evaluate1Form(const Buffer<double>& discreteForm, const Vector2& evaluationPoint) const {

	if (!mesh_old_ptr)
		return Vector2(0.0, 0.0); //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4

	//find the element containing evaluationPoint
	static uint element = 0;
	if (!mesh_old.findElement(evaluationPoint4, element))
		return Vector2(0.0, 0.0); //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getFaceNodesCustom(element);
	uint i, j;

	Buffer<uint> edges(3);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[1], nodes[2]);
	Buffer<VectorN> edgeCoefficients(edges.size()); //coefficients for small edges that are on big edges
	for (i = 0; i < edgeCoefficients.size(); ++i) {
		VectorN edgeValues(activeSmallEdgesInEdges.size(), 0.0);
		for (j = 0; j < activeSmallEdgesInEdges.size(); ++j) {
			edgeValues[j] = discreteForm[smallEdgesEdgeList(edges[i], j)];
		}
		solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeValues, edgeCoefficients[i]);
	}

	VectorN faceCoefficients; //coefficients for small edges that are on big faces
	VectorN faceValues(activeSmallEdgesInFaces.size(), 0.0);
	uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
	uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
	uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[1], nodes[2] });
	for (j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
		faceValues[j] += faceEdgeValuesEdge0[permutationEdge0][j].dot(edgeCoefficients[edges.findFirst(edges[0])]);
		faceValues[j] += faceEdgeValuesEdge1[permutationEdge1][j].dot(edgeCoefficients[edges.findFirst(edges[1])]);
		faceValues[j] += faceEdgeValuesEdge2[permutationEdge2][j].dot(edgeCoefficients[edges.findFirst(edges[2])]);
		faceValues[j] += discreteForm[smallEdgesFaceList(element, j)];
	}
	solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues, faceCoefficients);

	//find barycentric coordinates and their gradients at evaluationPoint
	double lambda0, lambda1, lambda2;
	Vector2 gradLambda0 = mesh_old.getEdgeDeviation(edges[2], mesh_old.getNodePosition(nodes[0])).toVector2();
	double lensq0 = gradLambda0.lensq();
	lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() / std::sqrt(lensq0);
	gradLambda0 /= lensq0;
	Vector2 gradLambda1 = mesh_old.getEdgeDeviation(edges[1], mesh_old.getNodePosition(nodes[1])).toVector2();
	double lensq1 = gradLambda1.lensq();
	lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() / std::sqrt(lensq1);
	gradLambda1 /= lensq1;
	Vector2 gradLambda2 = mesh_old.getEdgeDeviation(edges[0], mesh_old.getNodePosition(nodes[2])).toVector2();
	double lensq2 = gradLambda2.lensq();
	lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() / std::sqrt(lensq2);
	gradLambda2 /= lensq2;

	//compute the values of the lowest order Whitney forms
	Vector2 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector2 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector2 w2 = lambda1 * gradLambda2 - lambda2 * gradLambda1;

	//compute the result using the coefficients solved above
	Buffer<double> d(3, 0.0);
	std::vector<double> barycentricCoords;

	//small simplices of edges
	//edge0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1 }, edgePermutations[permutationEdge0]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[0] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[0][i];
	}
	//edge1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2 }, edgePermutations[permutationEdge1]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[1] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[1][i];
	}
	//edge2
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[permutationEdge2]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[2][i];
	}

	//take orientations into account
	if (getEdgePermutationSign(permutationEdge0) < 0)
		d[0] = -d[0];
	if (getEdgePermutationSign(permutationEdge1) < 0)
		d[1] = -d[1];
	if (getEdgePermutationSign(permutationEdge2) < 0)
		d[2] = -d[2];

	//small simplices in the interior
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * faceCoefficients[i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2;
}

/*
Computes the value of the (higher order) Whitney 2-form (interpolant of discreteForm) at evaluationPoint.
When evaluating many times, it is more efficient to first solve the coefficients with solve2FormCoefficients and then evaluate with evaluate2FormWithCoefficients.
*/
double SmallSimplexPartition2D::evaluate2Form(const Buffer<double>& discreteForm, const Vector2& evaluationPoint) const {

	if (!mesh_old_ptr)
		return 0.0; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4

	//find the element containing evaluationPoint
	static uint element = 0;
	if (!mesh_old.findElement(evaluationPoint4, element))
		return 0.0; //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getFaceNodesCustom(element);
	uint i;

	VectorN faceValues(activeSmallFacesInFaces.size(), 0.0);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		faceValues[i] += discreteForm[smallFacesFaceList(element, i)];
	}
	VectorN faceCoefficients; //coefficients for small faces in the interior of the tetrahedron
	solveLUP(matrix2FormsFaces, matrix2FormsFaces_p, faceValues, faceCoefficients);

	Buffer<uint> edges(3);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[1], nodes[2]);

	//find barycentric coordinates of evaluationPoint
	double h0 = mesh_old.getEdgeDeviation(edges[2], mesh_old.getNodePosition(nodes[0])).len();
	double lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() / h0;
	double h1 = mesh_old.getEdgeDeviation(edges[1], mesh_old.getNodePosition(nodes[1])).len();
	double lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() / h1;
	double h2 = mesh_old.getEdgeDeviation(edges[0], mesh_old.getNodePosition(nodes[2])).len();
	double lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() / h2;

	double result = 0.0; //compute the result using the coefficients solved above

	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		result += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1]) * std::pow(lambda2, mi[2]) * faceCoefficients[i];
	}
	return result / std::abs(mesh_old.getFaceVector2(element).determinant());
}

//solves the coefficients of the (higher order) Whitney 0-form (interpolant of discreteForm)
void SmallSimplexPartition2D::solve0FormCoefficients(const Buffer<double>& discreteForm, Buffer<double>& nodeCoefficients, Buffer<VectorN>& edgeCoefficients,
	Buffer<VectorN>& faceCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	uint i, j;

	nodeCoefficients.resize(mesh_old.getNodeSize()); //coefficients for small nodes that are also big nodes are obtained at once
	for (i = 0; i < nodeCoefficients.size(); ++i) {
		nodeCoefficients[i] = discreteForm[i];
	}

	edgeCoefficients.resize(mesh_old.getEdgeSize()); //coefficients for small nodes that are on big edges
	for (i = 0; i < edgeCoefficients.size(); ++i) {
		Buffer<uint> edgeNodes = mesh_old.getEdgeNodes(i);
		VectorN edgeValues = edgeValuesNode0 * discreteForm[edgeNodes[0]] + edgeValuesNode1 * discreteForm[edgeNodes[1]];
		for (j = 0; j < activeSmallNodesInEdges.size(); ++j) {
			edgeValues[j] += discreteForm[smallNodesEdgeList(i, j)];
		}
		solveLUP(matrix0FormsEdges, matrix0FormsEdges_p, edgeValues, edgeCoefficients[i]);
	}

	faceCoefficients.resize(mesh_old.getFaceSize()); //coefficients for small nodes that are on big faces
	for (i = 0; i < faceCoefficients.size(); ++i) {
		const Buffer<uint> faceNodes = getFaceNodesCustom(i);
		Buffer<uint> faceEdges(3);
		faceEdges[0] = mesh_old.findEdge(faceNodes[0], faceNodes[1]);
		faceEdges[1] = mesh_old.findEdge(faceNodes[0], faceNodes[2]);
		faceEdges[2] = mesh_old.findEdge(faceNodes[1], faceNodes[2]);
		VectorN faceValues = faceValuesNode0 * discreteForm[faceNodes[0]] + faceValuesNode1 * discreteForm[faceNodes[1]] + faceValuesNode2 * discreteForm[faceNodes[2]];
		uint permutationIndexEdge0 = getPermutation(mesh_old.getEdgeNodes(faceEdges[0]), { faceNodes[0], faceNodes[1] });
		uint permutationIndexEdge1 = getPermutation(mesh_old.getEdgeNodes(faceEdges[1]), { faceNodes[0], faceNodes[2] });
		uint permutationIndexEdge2 = getPermutation(mesh_old.getEdgeNodes(faceEdges[2]), { faceNodes[1], faceNodes[2] });
		for (j = 0; j < activeSmallNodesInFaces.size(); ++j) {
			faceValues[j] += faceNodeValuesEdge0[permutationIndexEdge0][j].dot(edgeCoefficients[faceEdges[0]]);
			faceValues[j] += faceNodeValuesEdge1[permutationIndexEdge1][j].dot(edgeCoefficients[faceEdges[1]]);
			faceValues[j] += faceNodeValuesEdge2[permutationIndexEdge2][j].dot(edgeCoefficients[faceEdges[2]]);
			faceValues[j] += discreteForm[smallNodesFaceList(i, j)];
		}
		solveLUP(matrix0FormsFaces, matrix0FormsFaces_p, faceValues, faceCoefficients[i]);
	}
}

//solves the coefficients of the (higher order) Whitney 1-form (interpolant of discreteForm)
void SmallSimplexPartition2D::solve1FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	uint i, j;

	edgeCoefficients.resize(mesh_old.getEdgeSize()); //coefficients for small edges that are on big edges
	for (i = 0; i < edgeCoefficients.size(); ++i) {
		VectorN edgeValues(activeSmallEdgesInEdges.size(), 0.0);
		for (j = 0; j < activeSmallEdgesInEdges.size(); ++j) {
			edgeValues[j] = discreteForm[smallEdgesEdgeList(i, j)];
		}
		solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeValues, edgeCoefficients[i]);
	}

	faceCoefficients.resize(mesh_old.getFaceSize()); //coefficients for small edges that are on big faces
	for (i = 0; i < faceCoefficients.size(); ++i) {
		const Buffer<uint> faceNodes = getFaceNodesCustom(i);
		Buffer<uint> faceEdges(3);
		faceEdges[0] = mesh_old.findEdge(faceNodes[0], faceNodes[1]);
		faceEdges[1] = mesh_old.findEdge(faceNodes[0], faceNodes[2]);
		faceEdges[2] = mesh_old.findEdge(faceNodes[1], faceNodes[2]);
		VectorN faceValues(activeSmallEdgesInFaces.size(), 0.0);
		uint permutationIndexEdge0 = getPermutation(mesh_old.getEdgeNodes(faceEdges[0]), { faceNodes[0], faceNodes[1] });
		uint permutationIndexEdge1 = getPermutation(mesh_old.getEdgeNodes(faceEdges[1]), { faceNodes[0], faceNodes[2] });
		uint permutationIndexEdge2 = getPermutation(mesh_old.getEdgeNodes(faceEdges[2]), { faceNodes[1], faceNodes[2] });
		for (j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
			faceValues[j] += faceEdgeValuesEdge0[permutationIndexEdge0][j].dot(edgeCoefficients[faceEdges[0]]);
			faceValues[j] += faceEdgeValuesEdge1[permutationIndexEdge1][j].dot(edgeCoefficients[faceEdges[1]]);
			faceValues[j] += faceEdgeValuesEdge2[permutationIndexEdge2][j].dot(edgeCoefficients[faceEdges[2]]);
			faceValues[j] += discreteForm[smallEdgesFaceList(i, j)];
		}
		solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues, faceCoefficients[i]);
	}
}

//solves the coefficients of the (higher order) Whitney 2-form (interpolant of discreteForm)
void SmallSimplexPartition2D::solve2FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& faceCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	uint i, j;

	faceCoefficients.resize(mesh_old.getFaceSize()); //coefficients for small faces that are on big faces
	for (i = 0; i < faceCoefficients.size(); ++i) {
		VectorN faceValues(activeSmallFacesInFaces.size(), 0.0);
		for (j = 0; j < activeSmallFacesInFaces.size(); ++j) {
			faceValues[j] = discreteForm[smallFacesFaceList(i, j)];
		}
		solveLUP(matrix2FormsFaces, matrix2FormsFaces_p, faceValues, faceCoefficients[i]);
	}
}

//computes the value of the (higher order) Whitney 0-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known
double SmallSimplexPartition2D::evaluate0FormWithCoefficients(const Vector2& evaluationPoint, const Buffer<double>& nodeCoefficients, const Buffer<VectorN>& edgeCoefficients,
	const Buffer<VectorN>& faceCoefficients, uint element) const {

	if (!mesh_old_ptr)
		return 0.0; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4
	uint i;

	//find the element containing evaluationPoint
	if (!mesh_old.findElement(evaluationPoint4, element))
		return 0.0; //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getFaceNodesCustom(element);
	Buffer<uint> edges(3);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[1], nodes[2]);

	uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
	uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
	uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[1], nodes[2] });

	//find barycentric coordinates of evaluationPoint
	double h0 = mesh_old.getEdgeDeviation(edges[2], mesh_old.getNodePosition(nodes[0])).len();
	double lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() / h0;
	double h1 = mesh_old.getEdgeDeviation(edges[1], mesh_old.getNodePosition(nodes[1])).len();
	double lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() / h1;
	double h2 = mesh_old.getEdgeDeviation(edges[0], mesh_old.getNodePosition(nodes[2])).len();
	double lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() / h2;

	//compute the result using the given coefficients
	std::vector<double> barycentricCoords;

	//contribution from small nodes that are also big nodes
	double result = std::pow(lambda0, order) * nodeCoefficients[nodes[0]] + std::pow(lambda1, order) * nodeCoefficients[nodes[1]]
		+ std::pow(lambda2, order) * nodeCoefficients[nodes[2]];

	//small nodes on edges
	//edge0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1 }, edgePermutations[permutationEdge0]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[0]][i];
	}
	//edge1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2 }, edgePermutations[permutationEdge1]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[1]][i];
	}
	//edge2
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[permutationEdge2]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[2]][i];
	}

	//small nodes in the interior
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(lambda0, mi[0] + 1) * std::pow(lambda1, mi[1]) * std::pow(lambda2, mi[2]) * faceCoefficients[element][i];
	}

	return result;
}

//computes the value of the (higher order) Whitney 1-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known
Vector2 SmallSimplexPartition2D::evaluate1FormWithCoefficients(const Vector2& evaluationPoint, const Buffer<VectorN>& edgeCoefficients,
	const Buffer<VectorN>& faceCoefficients, uint element) const {

	if (!mesh_old_ptr)
		return Vector2(0.0, 0.0); //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4
	uint i;

	//find the element containing evaluationPoint
	if (!mesh_old.findElement(evaluationPoint4, element))
		return Vector2(0.0, 0.0); //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getFaceNodesCustom(element);
	Buffer<uint> edges(3);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[1], nodes[2]);

	uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
	uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
	uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[1], nodes[2] });

	//find barycentric coordinates and their gradients at evaluationPoint
	double lambda0, lambda1, lambda2;
	Vector2 gradLambda0 = mesh_old.getEdgeDeviation(edges[2], mesh_old.getNodePosition(nodes[0])).toVector2();
	double lensq0 = gradLambda0.lensq();
	lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() / std::sqrt(lensq0);
	gradLambda0 /= lensq0;
	Vector2 gradLambda1 = mesh_old.getEdgeDeviation(edges[1], mesh_old.getNodePosition(nodes[1])).toVector2();
	double lensq1 = gradLambda1.lensq();
	lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() / std::sqrt(lensq1);
	gradLambda1 /= lensq1;
	Vector2 gradLambda2 = mesh_old.getEdgeDeviation(edges[0], mesh_old.getNodePosition(nodes[2])).toVector2();
	double lensq2 = gradLambda2.lensq();
	lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() / std::sqrt(lensq2);
	gradLambda2 /= lensq2;

	//compute the values of the lowest order Whitney forms
	Vector2 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector2 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector2 w2 = lambda1 * gradLambda2 - lambda2 * gradLambda1;

	//compute the result using the coefficients solved above
	Buffer<double> d(3, 0.0);
	std::vector<double> barycentricCoords;

	//small simplices of edges
	//edge0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1 }, edgePermutations[permutationEdge0]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[0] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[0]][i];
	}
	//edge1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2 }, edgePermutations[permutationEdge1]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[1] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[1]][i];
	}
	//edge2
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[permutationEdge2]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[2]][i];
	}

	//take orientations into account
	if (getEdgePermutationSign(permutationEdge0) < 0)
		d[0] = -d[0];
	if (getEdgePermutationSign(permutationEdge1) < 0)
		d[1] = -d[1];
	if (getEdgePermutationSign(permutationEdge2) < 0)
		d[2] = -d[2];

	//small simplices in the interior
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * faceCoefficients[element][i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2;
}

/*
Computes the value of the (higher order) Whitney 2-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known.
N.B. discreteForm[i] is interpreted with the positive sign as if getFaceVector(i).xy > 0.
*/
double SmallSimplexPartition2D::evaluate2FormWithCoefficients(const Vector2& evaluationPoint, const Buffer<VectorN>& faceCoefficients, uint element) const {

	if (!mesh_old_ptr)
		return 0.0; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4
	uint i;

	//find the element containing evaluationPoint
	if (!mesh_old.findElement(evaluationPoint4, element))
		return 0.0; //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getFaceNodesCustom(element);
	Buffer<uint> edges(3);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[1], nodes[2]);

	//find barycentric coordinates of evaluationPoint
	double h0 = mesh_old.getEdgeDeviation(edges[2], mesh_old.getNodePosition(nodes[0])).len();
	double lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() / h0;
	double h1 = mesh_old.getEdgeDeviation(edges[1], mesh_old.getNodePosition(nodes[1])).len();
	double lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() / h1;
	double h2 = mesh_old.getEdgeDeviation(edges[0], mesh_old.getNodePosition(nodes[2])).len();
	double lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() / h2;

	double result = 0.0; //compute the result using the given coefficients

	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		result += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1]) * std::pow(lambda2, mi[2]) * faceCoefficients[element][i];
	}
	return result / std::abs(mesh_old.getFaceVector2(element).determinant());
}

/*
Computes the discrete Hodge for higher order Whitney 1-forms to a sparse matrix (of type Buffer<std::unordered_map<uint, double>>).
Integrals over dual edges are computed in each element separately. If all elements have same shape, use formHodgeMatrix1FormsCustom for efficiency.
*/
void SmallSimplexPartition2D::formHodgeMatrix1Forms(Buffer<std::unordered_map<uint, double>>& star, bool circumcentric, bool minkowskiMetric) const {
	if (!mesh_old_ptr)
		return; //SmallSimplexPartition2D has not been associated with a mesh
	const BuilderMesh& mesh_old = *mesh_old_ptr;
	const BuilderMesh& mesh = *mesh_ptr;

	star.resize(mesh.getEdgeSize());

	struct DualEdgePart
	{
		uint index;
		Vector2 p0;
		Vector2 p1;
		Vector2 normalVector;
	};

	for (uint f = 0; f < mesh_old.getFaceSize(); ++f) {
		//find the subsimplices and their orientations
		const Buffer<uint> nodes = getFaceNodesCustom(f);
		Buffer<uint> edges(3);
		edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
		edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
		edges[2] = mesh_old.findEdge(nodes[1], nodes[2]);
		uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
		uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
		uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[1], nodes[2] });

		//save node positions, and outward normal vectors of the opposite edges, and gradients of barycentric functions
		Vector2 p0 = mesh_old.getNodePosition2(nodes[0]);
		Vector2 p1 = mesh_old.getNodePosition2(nodes[1]);
		Vector2 p2 = mesh_old.getNodePosition2(nodes[2]);
		Vector2 p0OppositeNormal = (p2 - p1).dual();
		if (p0OppositeNormal.dot(p1 - p0) < 0)
			p0OppositeNormal *= -1.0;
		Vector2 p1OppositeNormal = (p2 - p0).dual();
		if (p1OppositeNormal.dot(p0 - p1) < 0)
			p1OppositeNormal *= -1.0;
		Vector2 p2OppositeNormal = (p1 - p0).dual();
		if (p2OppositeNormal.dot(p0 - p2) < 0)
			p2OppositeNormal *= -1.0;
		Vector2 gradLambda0 = mesh_old.getEdgeDeviation(edges[2], mesh_old.getNodePosition(nodes[0])).toVector2();
		gradLambda0 /= gradLambda0.lensq();
		Vector2 gradLambda1 = mesh_old.getEdgeDeviation(edges[1], mesh_old.getNodePosition(nodes[1])).toVector2();
		gradLambda1 /= gradLambda1.lensq();
		Vector2 gradLambda2 = mesh_old.getEdgeDeviation(edges[0], mesh_old.getNodePosition(nodes[2])).toVector2();
		gradLambda2 /= gradLambda2.lensq();

		//form the list of dual edge parts contained in face f
		std::vector<DualEdgePart> edgeList;
		edgeList.reserve(3 * edgeMultiIndices.size() + 6 * faceEdgeHoles.size());
		for (uint e = 0; e < edges.size(); ++e) {
			for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
				uint edgeIndex = smallEdgesEdgeList(edges[e], i);
				const Buffer<uint>& smallFaces = mesh.getEdgeFaces(edgeIndex);
				uint sf = 0;
				while (!isInsideTriangle(mesh.getFaceAverage2(smallFaces[sf]), p0, p1, p2, p0OppositeNormal, p1OppositeNormal, p2OppositeNormal))
					++sf;
				Vector2 sf_bc = (circumcentric ? mesh.getFacePosition2(smallFaces[sf]) : mesh.getFaceAverage2(smallFaces[sf]));
				Vector2 se_bc = mesh.getEdgeAverage2(edgeIndex);
				Vector2 normalVector(sf_bc.y - se_bc.y, se_bc.x - sf_bc.x);
				if (minkowskiMetric)
					normalVector.y = -normalVector.y;
				if (mesh.getFaceIncidence(smallFaces[sf], edgeIndex) * mesh.getFaceDualVector2(smallFaces[sf]) < 0)
					normalVector *= -1.0;
				edgeList.push_back({ edgeIndex, sf_bc, se_bc, normalVector });
			}
		}
		for (uint i = firstEdgeOfFace(f); i < firstEdgeOfFace(f) + 3 * faceEdgeHoles.size(); ++i) {
			const Buffer<uint>& smallFaces = mesh.getEdgeFaces(i);
			Vector2 sf0_bc = (circumcentric ? mesh.getFacePosition2(smallFaces[0]) : mesh.getFaceAverage2(smallFaces[0]));
			Vector2 sf1_bc = (circumcentric ? mesh.getFacePosition2(smallFaces[1]) : mesh.getFaceAverage2(smallFaces[1]));
			Vector2 se_bc = mesh.getEdgeAverage2(i);
			Vector2 normalVector0(sf0_bc.y - se_bc.y, se_bc.x - sf0_bc.x);
			if (minkowskiMetric)
				normalVector0.y = -normalVector0.y;
			if (mesh.getFaceIncidence(smallFaces[0], i) * mesh.getFaceDualVector2(smallFaces[0]) < 0)
				normalVector0 *= -1.0;
			Vector2 normalVector1(sf1_bc.y - se_bc.y, se_bc.x - sf1_bc.x);
			if (minkowskiMetric)
				normalVector1.y = -normalVector1.y;
			if (mesh.getFaceIncidence(smallFaces[1], i) * mesh.getFaceDualVector2(smallFaces[1]) < 0)
				normalVector1 *= -1.0;
			edgeList.push_back({ i, sf0_bc, se_bc, normalVector0 });
			edgeList.push_back({ i, sf1_bc, se_bc, normalVector1 });
		}

		//integrate all basis functions over the edges in edgeList and add the values in Hodge matrix
		VectorN edgeCochain(activeSmallEdgesInEdges.size(), 0.0);
		VectorN faceCochain(activeSmallEdgesInFaces.size(), 0.0);
		VectorN edgeCoefficients(activeSmallEdgesInEdges.size(), 0.0);
		VectorN faceValues(activeSmallEdgesInFaces.size(), 0.0);
		VectorN faceCoefficients(activeSmallEdgesInFaces.size(), 0.0);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) { //edge basis functions
			edgeCochain[i] = 1.0;
			solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeCochain, edgeCoefficients);
			//edge 0
			for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
				faceValues[j] = faceEdgeValuesEdge0[permutationEdge0][j].dot(edgeCoefficients);
			}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues, faceCoefficients);
			for (uint j = 0; j < edgeList.size(); ++j) {
				std::function<double(Vector2)> interpolant{
						[&](Vector2 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 0, permutationEdge0, faceCoefficients, f, nodes, edges,
								gradLambda0, gradLambda1, gradLambda2).dot(edgeList[j].normalVector);
						}
				};
				star[edgeList[j].index][smallEdgesEdgeList(edges[0], i)] += gfd::integralAverage(interpolant, edgeList[j].p0, edgeList[j].p1);
			}
			faceValues.val.fill(0.0);
			faceCoefficients.val.fill(0.0);
			//edge1
			for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
				faceValues[j] = faceEdgeValuesEdge1[permutationEdge1][j].dot(edgeCoefficients);
			}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues, faceCoefficients);
			for (uint j = 0; j < edgeList.size(); ++j) {
				std::function<double(Vector2)> interpolant{
						[&](Vector2 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 1, permutationEdge1, faceCoefficients, f, nodes, edges,
								gradLambda0, gradLambda1, gradLambda2).dot(edgeList[j].normalVector);
						}
				};
				star[edgeList[j].index][smallEdgesEdgeList(edges[1], i)] += gfd::integralAverage(interpolant, edgeList[j].p0, edgeList[j].p1);
			}
			faceValues.val.fill(0.0);
			faceCoefficients.val.fill(0.0);
			//edge2
			for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
				faceValues[j] += faceEdgeValuesEdge2[permutationEdge2][j].dot(edgeCoefficients);
			}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues, faceCoefficients);
			for (uint j = 0; j < edgeList.size(); ++j) {
				std::function<double(Vector2)> interpolant{
						[&](Vector2 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 2, permutationEdge2, faceCoefficients, f, nodes, edges,
								gradLambda0, gradLambda1, gradLambda2).dot(edgeList[j].normalVector);
						}
				};
				star[edgeList[j].index][smallEdgesEdgeList(edges[2], i)] += gfd::integralAverage(interpolant, edgeList[j].p0, edgeList[j].p1);
			}
			faceValues.val.fill(0.0);
			faceCoefficients.val.fill(0.0);
			edgeCoefficients.val.fill(0.0);
			edgeCochain[i] = 0.0;
		}
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) { //face basis functions
			faceCochain[i] = 1.0;
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceCochain, faceCoefficients);
			for (uint j = 0; j < edgeList.size(); ++j) {
				std::function<double(Vector2)> interpolant{
						[&](Vector2 p) -> double {
							return evaluate1FormWithFaceCoefficients(p, faceCoefficients, nodes, edges,
								gradLambda0, gradLambda1, gradLambda2).dot(edgeList[j].normalVector);
						}
				};
				star[edgeList[j].index][smallEdgesFaceList(f, i)] += gfd::integralAverage(interpolant, edgeList[j].p0, edgeList[j].p1);
			}
			faceCoefficients.val.fill(0.0);
			faceCochain[i] = 0.0;
		}
	}
}

/*
Reimplementation of Mesh::getFaceNodes that defines a definite order of vertices for the specific element type currently in use.
When all elements are congruent, this enables one to compute the integrals for Hodge matrix only once in a single element.
*/
Buffer<uint> SmallSimplexPartition2D::getFaceNodesCustom(uint f) const {
	//implemention for right triangle
	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Buffer<uint> result(3);
	const Buffer<uint>& edges = mesh_old.getFaceEdges(f);
	double maxLensq = 0.0;
	uint maxInd = 0;
	for (uint i = 0; i < edges.size(); ++i) {
		const Buffer<uint>& nodes = mesh_old.getEdgeNodes(edges[i]);
		double lensq = (mesh_old.getNodePosition2(nodes[1]) - mesh_old.getNodePosition2(nodes[0])).lensq();
		if (lensq > maxLensq) {
			maxLensq = lensq;
			maxInd = i;
		}
	}
	const Buffer<uint>& baseEdgeNodes = mesh_old.getEdgeNodes(edges[maxInd]);
	const Buffer<uint>& sideEdgeNodes = mesh_old.getEdgeNodes(edges[(maxInd + 1) % 3]);
	if (sideEdgeNodes[0] != baseEdgeNodes[0] && sideEdgeNodes[0] != baseEdgeNodes[1])
		result[2] = sideEdgeNodes[0];
	else
		result[2] = sideEdgeNodes[1];
	const Vector2 p0 = mesh_old.getNodePosition2(baseEdgeNodes[0]);
	const Vector2 p1 = mesh_old.getNodePosition2(baseEdgeNodes[1]);
	const Vector2 p2 = mesh_old.getNodePosition2(result[2]);
	if ((p2 - p0).dot(Vector2(p0.y - p1.y, p1.x - p0.x)) > 0) {
		result[0] = baseEdgeNodes[0];
		result[1] = baseEdgeNodes[1];
	}
	else {
		result[0] = baseEdgeNodes[1];
		result[1] = baseEdgeNodes[0];
	}
	return result;

	//implemention for isosceles triangle with same base and height
	/*const BuilderMesh& mesh_old = *mesh_old_ptr;
	Buffer<uint> result(3);
	const Buffer<uint>& edges = mesh_old.getFaceEdges(f);
	double minLensq = 1.0e+300;
	uint minInd = 0;
	for (uint i = 0; i < edges.size(); ++i) {
		const Buffer<uint>& nodes = mesh_old.getEdgeNodes(edges[i]);
		double lensq = (mesh_old.getNodePosition2(nodes[1]) - mesh_old.getNodePosition2(nodes[0])).lensq();
		if (lensq < minLensq) {
			minLensq = lensq;
			minInd = i;
		}
	}
	const Buffer<uint>& baseEdgeNodes = mesh_old.getEdgeNodes(edges[minInd]);
	const Buffer<uint>& sideEdgeNodes = mesh_old.getEdgeNodes(edges[(minInd + 1) % 3]);
	if (sideEdgeNodes[0] != baseEdgeNodes[0] && sideEdgeNodes[0] != baseEdgeNodes[1])
		result[2] = sideEdgeNodes[0];
	else
		result[2] = sideEdgeNodes[1];
	const Vector2 p0 = mesh_old.getNodePosition2(baseEdgeNodes[0]);
	const Vector2 p1 = mesh_old.getNodePosition2(baseEdgeNodes[1]);
	const Vector2 p2 = mesh_old.getNodePosition2(result[2]);
	if ((p2 - p0).dot(Vector2(p0.y - p1.y, p1.x - p0.x)) > 0) {
		result[0] = baseEdgeNodes[0];
		result[1] = baseEdgeNodes[1];
	}
	else {
		result[0] = baseEdgeNodes[1];
		result[1] = baseEdgeNodes[0];
	}
	return result;*/
}

/*
Computes the discrete Hodge (with respect to Euclidean metric) for higher order Whitney 1-forms to a sparse matrix (of type Buffer<std::unordered_map<uint, double>>).
Integrals over dual edges are computed only once in single element.
Requires that all elements have same shape and getBodyNodesCustom has been implemented for that element shape.
Optionally the required integrals can be saved to or loaded from files in folder given by foldername.
*/
void SmallSimplexPartition2D::formHodgeMatrix1FormsCustom(Buffer<std::unordered_map<uint, double>>& star, std::string foldername, bool circumcentric, bool loadIntegralsFromFile) const {
	if (!mesh_old_ptr)
		return; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Buffer<uint> refElementIndices(1, 0);
	Buffer<Buffer<uint>> refElements(1);
	refElements[0].resize(mesh_old.getFaceSize());
	for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
		refElements[0][i] = i;
	}
	formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, foldername, circumcentric, false, loadIntegralsFromFile);
}

/*
Computes the discrete Hodge for higher order Whitney 1-forms to a sparse matrix (of type Buffer<std::unordered_map<uint, double>>).
Integrals over dual edges are computed once in each reference element whose index is given in refElementIndices.
Requires that every element should be in one refElements buffer; the computations done in the corresponding reference element are used for all elements in the buffer.
Optionally the required integrals can be saved to or loaded from files in folder given by foldername.
*/
void SmallSimplexPartition2D::formHodgeMatrix1FormsCustom(Buffer<std::unordered_map<uint, double>>& star, const Buffer<uint>& refElementIndices, const Buffer<Buffer<uint>>& refElements,
	std::string foldername, bool circumcentric, bool minkowskiMetric, bool loadIntegralsFromFile) const {
	if (!mesh_old_ptr)
		return; //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	const BuilderMesh& mesh = *mesh_ptr;
	star.resize(mesh.getEdgeSize());
	const uint dualEdgeCount = 3 * faceMultiIndices.size();

	for (uint ref = 0; ref < refElementIndices.size(); ++ref) {
		Buffer<VectorN> hodgeIntegralsFaceEdges(dualEdgeCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge0Edges(dualEdgeCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge1Edges(dualEdgeCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge2Edges(dualEdgeCount);

		if (loadIntegralsFromFile) { //load the integrals from file
			loadHodgeIntegrals1Forms(foldername, circumcentric, hodgeIntegralsFaceEdges, hodgeIntegralsEdge0Edges, hodgeIntegralsEdge1Edges,
				hodgeIntegralsEdge2Edges, ref);
		}

		else { //compute the integrals in the reference element
			uint f = refElementIndices[ref];
			const Buffer<uint> face_f_Nodes = getFaceNodesCustom(f);
			Buffer<uint> face_f_Edges(3);
			face_f_Edges[0] = mesh_old.findEdge(face_f_Nodes[0], face_f_Nodes[1]);
			face_f_Edges[1] = mesh_old.findEdge(face_f_Nodes[0], face_f_Nodes[2]);
			face_f_Edges[2] = mesh_old.findEdge(face_f_Nodes[1], face_f_Nodes[2]);

			//save node positions, and outward normal vectors of the opposite edges, and gradients of barycentric functions
			Vector2 p0 = mesh_old.getNodePosition2(face_f_Nodes[0]);
			Vector2 p1 = mesh_old.getNodePosition2(face_f_Nodes[1]);
			Vector2 p2 = mesh_old.getNodePosition2(face_f_Nodes[2]);
			Vector2 p0OppositeNormal = (p2 - p1).dual();
			if (p0OppositeNormal.dot(p1 - p0) < 0)
				p0OppositeNormal *= -1.0;
			Vector2 p1OppositeNormal = (p2 - p0).dual();
			if (p1OppositeNormal.dot(p0 - p1) < 0)
				p1OppositeNormal *= -1.0;
			Vector2 p2OppositeNormal = (p1 - p0).dual();
			if (p2OppositeNormal.dot(p0 - p2) < 0)
				p2OppositeNormal *= -1.0;
			Vector2 gradLambda0 = mesh_old.getEdgeDeviation(face_f_Edges[2], mesh_old.getNodePosition(face_f_Nodes[0])).toVector2();
			gradLambda0 /= gradLambda0.lensq();
			Vector2 gradLambda1 = mesh_old.getEdgeDeviation(face_f_Edges[1], mesh_old.getNodePosition(face_f_Nodes[1])).toVector2();
			gradLambda1 /= gradLambda1.lensq();
			Vector2 gradLambda2 = mesh_old.getEdgeDeviation(face_f_Edges[0], mesh_old.getNodePosition(face_f_Nodes[2])).toVector2();
			gradLambda2 /= gradLambda2.lensq();

			//form the list of dual edges
			struct DualEdgePart
			{
				Vector2 p0;
				Vector2 p1;
				Vector2 normalVector;
			};
			std::vector<std::vector<DualEdgePart>> dualEdgeList;
			dualEdgeList.resize(dualEdgeCount);
			uint deIndex = 0;
			for (uint e = 0; e < 3; ++e) { //dual edges of small edges on edges
				for (uint i = 0; i < edgeMultiIndices.size(); ++i) {
					mi_t& mi = edgeMultiIndices[i];
					uint n0Index, n1Index;
					if (e == 0) {
						n0Index = getSmallNodeIndex({ mi[0], mi[1], 0 }, 0, face_f_Nodes);
						n1Index = getSmallNodeIndex({ mi[0], mi[1], 0 }, 1, face_f_Nodes);
					}
					else if (e == 1) {
						n0Index = getSmallNodeIndex({ mi[0], 0, mi[1] }, 0, face_f_Nodes);
						n1Index = getSmallNodeIndex({ mi[0], 0, mi[1] }, 2, face_f_Nodes);
					}
					else if (e == 2) {
						n0Index = getSmallNodeIndex({ 0, mi[0], mi[1] }, 1, face_f_Nodes);
						n1Index = getSmallNodeIndex({ 0, mi[0], mi[1] }, 2, face_f_Nodes);
					}
					uint edgeIndex = mesh.findEdge(n0Index, n1Index);
					const Buffer<uint>& smallFaces = mesh.getEdgeFaces(edgeIndex);
					uint sf = 0;
					while (!isInsideTriangle(mesh.getFaceAverage2(smallFaces[sf]), p0, p1, p2, p0OppositeNormal, p1OppositeNormal, p2OppositeNormal))
						++sf;
					Vector2 sf_bc = (circumcentric ? mesh.getFacePosition2(smallFaces[sf]) : mesh.getFaceAverage2(smallFaces[sf]));
					Vector2 se_bc = mesh.getEdgeAverage2(edgeIndex);
					Vector2 normalVector(sf_bc.y - se_bc.y, se_bc.x - sf_bc.x);
					if (minkowskiMetric)
						normalVector.y = -normalVector.y;
					if (e == 1)
						normalVector *= -1.0;
					dualEdgeList[deIndex].push_back({ sf_bc, se_bc, normalVector });
					++deIndex;
				}
			}
			for (uint e = 0; e < 3; ++e) { //dual edges of small edges in the interior
				for (uint i = 0; i < faceEdgeHoles.size(); ++i) {
					mi_t mi;
					uint n0Index;
					uint n1Index;
					if (e == 0) {
						mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
						n0Index = getSmallNodeIndex(mi, 0, face_f_Nodes);
						n1Index = getSmallNodeIndex(mi, 1, face_f_Nodes);
					}
					else if (e == 1) {
						mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
						n0Index = getSmallNodeIndex(mi, 0, face_f_Nodes);
						n1Index = getSmallNodeIndex(mi, 2, face_f_Nodes);
					}
					else if (e == 2) {
						mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
						n0Index = getSmallNodeIndex(mi, 1, face_f_Nodes);
						n1Index = getSmallNodeIndex(mi, 2, face_f_Nodes);
					}
					uint edgeIndex = mesh.findEdge(n0Index, n1Index);
					const Buffer<uint>& smallFaces = mesh.getEdgeFaces(edgeIndex);
					Vector2 sf0_bc = (circumcentric ? mesh.getFacePosition2(smallFaces[0]) : mesh.getFaceAverage2(smallFaces[0]));
					Vector2 sf1_bc = (circumcentric ? mesh.getFacePosition2(smallFaces[1]) : mesh.getFaceAverage2(smallFaces[1]));
					Vector2 se_bc = mesh.getEdgeAverage2(edgeIndex);
					Vector2 normalVector0(sf0_bc.y - se_bc.y, se_bc.x - sf0_bc.x);
					if (minkowskiMetric)
						normalVector0.y = -normalVector0.y;
					if (mesh.getFaceIncidence(smallFaces[0], edgeIndex) * mesh.getFaceDualVector2(smallFaces[0]) < 0)
						normalVector0 *= -1.0;
					Vector2 normalVector1(sf1_bc.y - se_bc.y, se_bc.x - sf1_bc.x);
					if (minkowskiMetric)
						normalVector1.y = -normalVector1.y;
					if (mesh.getFaceIncidence(smallFaces[1], edgeIndex) * mesh.getFaceDualVector2(smallFaces[1]) < 0)
						normalVector1 *= -1.0;
					dualEdgeList[deIndex].push_back({ sf0_bc, se_bc, normalVector0 });
					dualEdgeList[deIndex].push_back({ sf1_bc, se_bc, normalVector1 });

					++deIndex;
				}
			}

			//compute the integrals of basis functions
			for (uint i = 0; i < dualEdgeCount; ++i) {
				hodgeIntegralsFaceEdges[i].toVectorN(activeSmallEdgesInFaces.size());
				hodgeIntegralsEdge0Edges[i].resize(edgePermutations.size());
				hodgeIntegralsEdge1Edges[i].resize(edgePermutations.size());
				hodgeIntegralsEdge2Edges[i].resize(edgePermutations.size());
				for (uint j = 0; j < edgePermutations.size(); ++j) {
					hodgeIntegralsEdge0Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
					hodgeIntegralsEdge1Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
					hodgeIntegralsEdge2Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
				}
			}

			//first integrate face basis functions
			VectorN faceCochain(activeSmallEdgesInFaces.size(), 0.0);
			VectorN faceCoefficients(activeSmallEdgesInFaces.size(), 0.0);
			for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
				faceCochain[i] = 1.0;
				solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceCochain, faceCoefficients);
				for (uint j = 0; j < dualEdgeCount; ++j) {
					for (uint k = 0; k < dualEdgeList[j].size(); ++k) {
						std::function<double(Vector2)> interpolant{
							[&](Vector2 p) -> double {
								return evaluate1FormWithFaceCoefficients(p, faceCoefficients, face_f_Nodes, face_f_Edges,
									gradLambda0, gradLambda1, gradLambda2).dot(dualEdgeList[j][k].normalVector);
							}
						};
						hodgeIntegralsFaceEdges[j][i] += gfd::integralAverage(interpolant, dualEdgeList[j][k].p0, dualEdgeList[j][k].p1);
					}

				}
				faceCoefficients.val.fill(0.0);
				faceCochain[i] = 0.0;
			}

			//next integrate edge basis functions for each edge and orientation
			VectorN edgeCochain(activeSmallEdgesInEdges.size(), 0.0);
			VectorN edgeCoefficients(activeSmallEdgesInEdges.size(), 0.0);
			VectorN faceValues(activeSmallEdgesInFaces.size(), 0.0);
			for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
				edgeCochain[i] = 1.0;
				solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeCochain, edgeCoefficients);
				for (uint j = 0; j < edgePermutations.size(); ++j) {
					//edge0
					for (uint k = 0; k < activeSmallEdgesInFaces.size(); ++k) {
						faceValues[k] += faceEdgeValuesEdge0[j][k].dot(edgeCoefficients);
					}
					for (uint k = 0; k < dualEdgeCount; ++k) {
						for (uint l = 0; l < dualEdgeList[k].size(); ++l) {
							std::function<double(Vector2)> interpolant{
								[&](Vector2 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients, 0, j, face_f_Nodes, face_f_Edges,
										gradLambda0, gradLambda1, gradLambda2).dot(dualEdgeList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge0Edges[k][j][i] += gfd::integralAverage(interpolant, dualEdgeList[k][l].p0, dualEdgeList[k][l].p1);
						}
						hodgeIntegralsEdge0Edges[k][j][i] += faceValues.dot(hodgeIntegralsFaceEdges[k]);
					}
					faceValues.val.fill(0.0);
					//edge1
					for (uint k = 0; k < activeSmallEdgesInFaces.size(); ++k) {
						faceValues[k] += faceEdgeValuesEdge1[j][k].dot(edgeCoefficients);
					}
					for (uint k = 0; k < dualEdgeCount; ++k) {
						for (uint l = 0; l < dualEdgeList[k].size(); ++l) {
							std::function<double(Vector2)> interpolant{
								[&](Vector2 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients, 1, j, face_f_Nodes, face_f_Edges,
										gradLambda0, gradLambda1, gradLambda2).dot(dualEdgeList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge1Edges[k][j][i] += gfd::integralAverage(interpolant, dualEdgeList[k][l].p0, dualEdgeList[k][l].p1);
						}
						hodgeIntegralsEdge1Edges[k][j][i] += faceValues.dot(hodgeIntegralsFaceEdges[k]);
					}
					faceValues.val.fill(0.0);
					//edge2
					for (uint k = 0; k < activeSmallEdgesInFaces.size(); ++k) {
						faceValues[k] += faceEdgeValuesEdge2[j][k].dot(edgeCoefficients);
					}
					for (uint k = 0; k < dualEdgeCount; ++k) {
						for (uint l = 0; l < dualEdgeList[k].size(); ++l) {
							std::function<double(Vector2)> interpolant{
								[&](Vector2 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients, 2, j, face_f_Nodes, face_f_Edges,
										gradLambda0, gradLambda1, gradLambda2).dot(dualEdgeList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge2Edges[k][j][i] += gfd::integralAverage(interpolant, dualEdgeList[k][l].p0, dualEdgeList[k][l].p1);
						}
						hodgeIntegralsEdge2Edges[k][j][i] += faceValues.dot(hodgeIntegralsFaceEdges[k]);
					}
					faceValues.val.fill(0.0);
				}
				edgeCoefficients.val.fill(0.0);
				edgeCochain[i] = 0.0;
			}
		}

		if (!loadIntegralsFromFile) { //save the integrals to file
			saveHodgeIntegrals1Forms(foldername, circumcentric, hodgeIntegralsFaceEdges, hodgeIntegralsEdge0Edges, hodgeIntegralsEdge1Edges,
				hodgeIntegralsEdge2Edges, ref);
		}

		//form Hodge matrix from the precomputed integrals
		for (uint face = 0; face < refElements[ref].size(); ++face) {
			uint f = refElements[ref][face];
			//find the subsimplices and their orientations
			const Buffer<uint> nodes = getFaceNodesCustom(f);
			Buffer<uint> edges(3);
			edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
			edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
			edges[2] = mesh_old.findEdge(nodes[1], nodes[2]);
			uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
			uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
			uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[1], nodes[2] });

			//save node positions, outward normal vectors of the opposite edges, and the scaling factor in current face
			Vector2 p0 = mesh_old.getNodePosition2(nodes[0]);
			Vector2 p1 = mesh_old.getNodePosition2(nodes[1]);
			Vector2 p2 = mesh_old.getNodePosition2(nodes[2]);
			Vector2 p0OppositeNormal = (p2 - p1).dual();
			if (p0OppositeNormal.dot(p1 - p0) < 0)
				p0OppositeNormal *= -1.0;
			Vector2 p1OppositeNormal = (p2 - p0).dual();
			if (p1OppositeNormal.dot(p0 - p1) < 0)
				p1OppositeNormal *= -1.0;
			Vector2 p2OppositeNormal = (p1 - p0).dual();
			if (p2OppositeNormal.dot(p0 - p2) < 0)
				p2OppositeNormal *= -1.0;
			double scalingConstant = 1.0;

			//find the indices of all simplices in current body and their orientations with respect to what was precomputed
			Buffer<uint> indices(dualEdgeCount);
			Buffer<double> incidences(dualEdgeCount);
			uint deIndex = 0;
			for (uint e = 0; e < 3; ++e) { //small edges on edges
				for (uint i = 0; i < edgeMultiIndices.size(); ++i) {
					mi_t& mi = edgeMultiIndices[i];
					uint n0Index, n1Index;
					if (e == 0) {
						n0Index = getSmallNodeIndex({ mi[0], mi[1], 0 }, 0, nodes);
						n1Index = getSmallNodeIndex({ mi[0], mi[1], 0 }, 1, nodes);
					}
					else if (e == 1) {
						n0Index = getSmallNodeIndex({ mi[0], 0, mi[1] }, 0, nodes);
						n1Index = getSmallNodeIndex({ mi[0], 0, mi[1] }, 2, nodes);
					}
					else if (e == 2) {
						n0Index = getSmallNodeIndex({ 0, mi[0], mi[1] }, 1, nodes);
						n1Index = getSmallNodeIndex({ 0, mi[0], mi[1] }, 2, nodes);
					}
					uint edgeIndex = mesh.findEdge(n0Index, n1Index);
					indices[deIndex] = edgeIndex;
					incidences[deIndex] = mesh.getEdgeIncidence(edgeIndex, n1Index) * scalingConstant;
					++deIndex;
				}
			}
			for (uint e = 0; e < 3; ++e) { //small edges in the interior
				for (uint i = 0; i < faceEdgeHoles.size(); ++i) {
					mi_t mi;
					uint n0Index;
					uint n1Index;
					if (e == 0) {
						mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
						n0Index = getSmallNodeIndex(mi, 0, nodes);
						n1Index = getSmallNodeIndex(mi, 1, nodes);
					}
					else if (e == 1) {
						mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
						n0Index = getSmallNodeIndex(mi, 0, nodes);
						n1Index = getSmallNodeIndex(mi, 2, nodes);
					}
					else if (e == 2) {
						mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
						n0Index = getSmallNodeIndex(mi, 1, nodes);
						n1Index = getSmallNodeIndex(mi, 2, nodes);
					}
					uint edgeIndex = mesh.findEdge(n0Index, n1Index);
					indices[deIndex] = edgeIndex;
					incidences[deIndex] = mesh.getEdgeIncidence(edgeIndex, n1Index) * scalingConstant;
					++deIndex;
				}
			}

			//get the integrals of all basis functions
			for (uint i = 0; i < dualEdgeCount; ++i) {
				for (uint j = 0; j < activeSmallEdgesInEdges.size(); ++j) { //edge basis functions
					star[indices[i]][smallEdgesEdgeList(edges[0], j)] += incidences[i] * hodgeIntegralsEdge0Edges[i][permutationEdge0][j];
					star[indices[i]][smallEdgesEdgeList(edges[1], j)] += incidences[i] * hodgeIntegralsEdge1Edges[i][permutationEdge1][j];
					star[indices[i]][smallEdgesEdgeList(edges[2], j)] += incidences[i] * hodgeIntegralsEdge2Edges[i][permutationEdge2][j];
				}
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) { //face basis functions
					star[indices[i]][smallEdgesFaceList(f, j)] += incidences[i] * hodgeIntegralsFaceEdges[i][j];
				}
			}
		}
	}
}

/*
The following functions return the index of the active small simplex j of the big simplex i.
Assumes that smallEdgesFaceOffsets have been stored in face 0 of the initial mesh when it was refined.
Then the indices of small simplices can be recovered from those of the big simplices and need not be explicitly saved in memory.
*/

uint SmallSimplexPartition2D::smallNodesEdgeList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + i * activeSmallNodesInEdges.size() + j;
}
uint SmallSimplexPartition2D::smallEdgesEdgeList(uint i, uint j) const {
	return i * edgeMultiIndices.size() + j;
}
uint SmallSimplexPartition2D::smallNodesFaceList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + mesh_old_ptr->getEdgeSize() * activeSmallNodesInEdges.size() + i * activeSmallNodesInFaces.size() + j;
}
uint SmallSimplexPartition2D::smallEdgesFaceList(uint i, uint j) const {
	return i * 3 * faceEdgeHoles.size() + smallEdgesFaceOffsets[j];
}
uint SmallSimplexPartition2D::smallFacesFaceList(uint i, uint j) const {
	return i * (faceMultiIndices.size() + faceEdgeHoles.size()) + j;
}

//Returns the smallest index of small edges (not necessarily active) in the interior of the given big face
uint SmallSimplexPartition2D::firstEdgeOfFace(uint i) const {
	return mesh_old_ptr->getEdgeSize() * edgeMultiIndices.size() + i * 3 * faceEdgeHoles.size();
}

//returns the multi-index that maps the given edge to the given hole boundary
mi_t* SmallSimplexPartition2D::getHoleEdgeMultiIndex(const mi_t& holeIndex, uint edgeIndex) const {
	mi_t multiIndex{ holeIndex };
	uint i;
	for (i = 0; i < multiIndex.size(); ++i) {
		if (faceEdges[edgeIndex][i] == 0)
			++multiIndex[i];
	}
	return &faceMultiIndices[indexOf(faceMultiIndices, multiIndex)];
}

//returns the index of the small node that is the image of multiIndex when we map the given node
uint SmallSimplexPartition2D::getSmallNodeIndex(mi_t multiIndex, uint nodeIndex, Buffer<uint> bigSimplexNodes) const {
	const BuilderMesh& mesh_old = *mesh_old_ptr;

	//remove redundant nodes from the big simplex
	for (uint i = multiIndex.size(); i-- > 0; ) {
		if (multiIndex[i] == 0 && nodeIndex != i) {
			bigSimplexNodes.erase(i);
			multiIndex.erase(multiIndex.begin() + i);
			if (i < nodeIndex)
				--nodeIndex;
		}
	}

	if (bigSimplexNodes.size() == 1)
		return bigSimplexNodes[0]; //small node coincides with big node

	//small node is now in the interior of the big simplex

	uint bigSimplexIndex; //index of the big simplex (could also use findSimplexWithNodes(bigSimplexNodes, mesh_old) here)
	Buffer<uint> nodesOrdered; //nodes of the big simplex in the correct order
	switch (bigSimplexNodes.size()) {
	case 2:
		bigSimplexIndex = mesh_old.findEdge(bigSimplexNodes[0], bigSimplexNodes[1]);
		nodesOrdered = mesh_old.getEdgeNodes(bigSimplexIndex);
		break;
	default:
		bigSimplexIndex = mesh_old.findFace(findEdges(bigSimplexNodes[0], bigSimplexNodes[1], bigSimplexNodes[2], mesh_old));
		nodesOrdered = getFaceNodesCustom(bigSimplexIndex);
		break;
	}

	mi_t multiIndexPermuted(multiIndex.size()); //permutation of multiIndex that corresponds to the correct order of nodes
	uint nodeIndexPermuted{}; //index of the given node in the correct order
	for (uint i = 0; i < multiIndex.size(); ++i) {
		for (uint j = 0; j < multiIndex.size(); ++j) {
			if (nodesOrdered[i] == bigSimplexNodes[j]) {
				multiIndexPermuted[i] = multiIndex[j];
				if (j == nodeIndex)
					nodeIndexPermuted = i;
			}
		}
	}

	if (nodeIndexPermuted != 0) { //get the small node as image of node 0
		multiIndexPermuted[0] -= 1;
		multiIndexPermuted[nodeIndexPermuted] += 1;
		nodeIndexPermuted = 0;
	}

	switch (multiIndexPermuted.size()) {
	case 2:
		return smallNodesEdgeList(bigSimplexIndex, indexOf(activeSmallNodesInEdges, multiIndexPermuted));
	case 3:
		return smallNodesFaceList(bigSimplexIndex, indexOf(activeSmallNodesInFaces, multiIndexPermuted));
	}
}

//computes the value resulting from the coefficients on the given edge in the element whose nodes and edges are given
Vector2 SmallSimplexPartition2D::evaluate1FormWithEdgeCoefficients(const Vector2& evaluationPoint, const VectorN& edgeCoefficients, uint edgeIndex, uint edgePermutation,
	const Buffer<uint>& nodes, const Buffer<uint>& edges, const Vector2& gradLambda0, const Vector2& gradLambda1, const Vector2& gradLambda2) const {
	if (!mesh_old_ptr)
		return Vector2(0.0, 0.0); //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4

	//find barycentric coordinates
	double lambda0, lambda1, lambda2;
	lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() * gradLambda0.len();
	lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() * gradLambda1.len();
	lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() * gradLambda2.len();

	//compute the values of the lowest order Whitney forms
	Vector2 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector2 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector2 w2 = lambda1 * gradLambda2 - lambda2 * gradLambda1;

	//compute the result using the coefficients solved above
	Buffer<double> d(3, 0.0);
	std::vector<double> barycentricCoords;

	//small simplices of edges
	if (edgeIndex == 0) {
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[0] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[0] = -d[0];
	}
	else if (edgeIndex == 1) {
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[1] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[1] = -d[1];
	}
	else if (edgeIndex == 2) {
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[2] = -d[2];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2;
}

//computes the value resulting from the given coefficients in the element whose nodes and edges are given
Vector2 SmallSimplexPartition2D::evaluate1FormWithCoefficients(const Vector2& evaluationPoint, const VectorN& edgeCoefficients, uint edgeIndex, uint edgePermutation, const VectorN&
	faceCoefficients, uint element, const Buffer<uint>& nodes, const Buffer<uint>& edges, const Vector2& gradLambda0, const Vector2& gradLambda1, const Vector2& gradLambda2) const {
	if (!mesh_old_ptr)
		return Vector2(0.0, 0.0); //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4

	//find barycentric coordinates
	double lambda0, lambda1, lambda2;
	lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() * gradLambda0.len();
	lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() * gradLambda1.len();
	lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() * gradLambda2.len();

	//compute the values of the lowest order Whitney forms
	Vector2 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector2 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector2 w2 = lambda1 * gradLambda2 - lambda2 * gradLambda1;

	//compute the result using the coefficients solved above
	Buffer<double> d(3, 0.0);
	std::vector<double> barycentricCoords;

	//small simplices of edges
	if (edgeIndex == 0) {
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[0] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[0] = -d[0];
	}
	else if (edgeIndex == 1) {
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[1] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[1] = -d[1];
	}
	else if (edgeIndex == 2) {
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[2] = -d[2];
	}

	//small simplices in the interior
	for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * faceCoefficients[i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2;
}

//computes the value resulting from the given coefficients in the element whose nodes and edges are given
Vector2 SmallSimplexPartition2D::evaluate1FormWithFaceCoefficients(const Vector2& evaluationPoint, const VectorN& faceCoefficients,
	const Buffer<uint>& nodes, const Buffer<uint>& edges, const Vector2& gradLambda0, const Vector2& gradLambda1, const Vector2& gradLambda2) const {
	if (!mesh_old_ptr)
		return Vector2(0.0, 0.0); //SmallSimplexPartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0, 0.0); //for routines that require Vector4

	//find barycentric coordinates
	double lambda0, lambda1, lambda2;
	lambda0 = mesh_old.getEdgeDeviation(edges[2], evaluationPoint4).len() * gradLambda0.len();
	lambda1 = mesh_old.getEdgeDeviation(edges[1], evaluationPoint4).len() * gradLambda1.len();
	lambda2 = mesh_old.getEdgeDeviation(edges[0], evaluationPoint4).len() * gradLambda2.len();

	//compute the values of the lowest order Whitney forms
	Vector2 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector2 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector2 w2 = lambda1 * gradLambda2 - lambda2 * gradLambda1;

	//compute the result using the coefficients solved above
	Buffer<double> d(3, 0.0);
	std::vector<double> barycentricCoords;

	//small simplices in the interior
	for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * faceCoefficients[i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2;
}

//precomputes the matrices required in evaluating higher order Whitney forms
void SmallSimplexPartition2D::formMatrices() {

	uint i, j, k;

	//0-forms
	Buffer<Buffer<double>> integrals0FormsEdges(activeSmallNodesInEdges.size());
	Buffer<Buffer<double>> integrals0FormsFaces(activeSmallNodesInFaces.size());
	const mi_t node0EdgeMultiIndex{ {order - 1, 0} };
	const mi_t node1EdgeMultiIndex{ {0, order - 1} };
	const mi_t node0FaceMultiIndex{ {order - 1, 0, 0} };
	const mi_t node1FaceMultiIndex{ {0, order - 1, 0} };
	const mi_t node2FaceMultiIndex{ {0, 0, order - 1} };

	//0-forms in edges
	edgeValuesNode0.toVectorN(activeSmallNodesInEdges.size());
	edgeValuesNode1.toVectorN(activeSmallNodesInEdges.size());
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		edgeValuesNode0[i] = -integrateWhitneyForm(*activeSmallNodesInEdges[i].multiIndex, *activeSmallNodesInEdges[i].nodeIndices,
			node0EdgeMultiIndex, edgeNodes[0]);
		edgeValuesNode1[i] = -integrateWhitneyForm(*activeSmallNodesInEdges[i].multiIndex, *activeSmallNodesInEdges[i].nodeIndices,
			node1EdgeMultiIndex, edgeNodes[1]);
		integrals0FormsEdges[i].resize(activeSmallNodesInEdges.size());
		for (j = 0; j < activeSmallNodesInEdges.size(); ++j) {
			integrals0FormsEdges[i][j] = integrateWhitneyForm(*activeSmallNodesInEdges[i].multiIndex, *activeSmallNodesInEdges[i].nodeIndices,
				*activeSmallNodesInEdges[j].multiIndex, *activeSmallNodesInEdges[j].nodeIndices);
		}
	}

	//0-forms in faces
	faceValuesNode0.toVectorN(activeSmallNodesInFaces.size());
	faceValuesNode1.toVectorN(activeSmallNodesInFaces.size());
	faceValuesNode2.toVectorN(activeSmallNodesInFaces.size());
	faceNodeValuesEdge0.resize(edgePermutations.size());
	faceNodeValuesEdge1.resize(edgePermutations.size());
	faceNodeValuesEdge2.resize(edgePermutations.size());
	for (i = 0; i < edgePermutations.size(); ++i) {
		faceNodeValuesEdge0[i].resize(activeSmallNodesInFaces.size());
		faceNodeValuesEdge1[i].resize(activeSmallNodesInFaces.size());
		faceNodeValuesEdge2[i].resize(activeSmallNodesInFaces.size());
		for (j = 0; j < activeSmallNodesInFaces.size(); ++j) {
			faceNodeValuesEdge0[i][j].toVectorN(activeSmallNodesInEdges.size());
			faceNodeValuesEdge1[i][j].toVectorN(activeSmallNodesInEdges.size());
			faceNodeValuesEdge2[i][j].toVectorN(activeSmallNodesInEdges.size());
		}
	}
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		faceValuesNode0[i] = -integrateWhitneyForm(*activeSmallNodesInFaces[i].multiIndex, *activeSmallNodesInFaces[i].nodeIndices,
			node0FaceMultiIndex, faceNodes[0]);
		faceValuesNode1[i] = -integrateWhitneyForm(*activeSmallNodesInFaces[i].multiIndex, *activeSmallNodesInFaces[i].nodeIndices,
			node1FaceMultiIndex, faceNodes[1]);
		faceValuesNode2[i] = -integrateWhitneyForm(*activeSmallNodesInFaces[i].multiIndex, *activeSmallNodesInFaces[i].nodeIndices,
			node2FaceMultiIndex, faceNodes[2]);
		for (j = 0; j < activeSmallNodesInEdges.size(); ++j) {
			for (k = 0; k < edgePermutations.size(); ++k) {
				faceNodeValuesEdge0[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInFaces[i].multiIndex, *activeSmallNodesInFaces[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 2, 0), insert(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 2, 0));
				faceNodeValuesEdge1[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInFaces[i].multiIndex, *activeSmallNodesInFaces[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 1, 0), insert(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 1, 0));
				faceNodeValuesEdge2[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInFaces[i].multiIndex, *activeSmallNodesInFaces[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 0, 0), insert(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 0, 0));
			}
		}
		integrals0FormsFaces[i].resize(activeSmallNodesInFaces.size());
		for (j = 0; j < activeSmallNodesInFaces.size(); ++j) {
			integrals0FormsFaces[i][j] = integrateWhitneyForm(*activeSmallNodesInFaces[i].multiIndex, *activeSmallNodesInFaces[i].nodeIndices,
				*activeSmallNodesInFaces[j].multiIndex, *activeSmallNodesInFaces[j].nodeIndices);
		}
	}

	//1-forms
	Buffer<Buffer<double>> integrals1FormsEdges(activeSmallEdgesInEdges.size());
	Buffer<Buffer<double>> integrals1FormsFaces(activeSmallEdgesInFaces.size());

	//1-forms in edges
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		integrals1FormsEdges[i].resize(activeSmallEdgesInEdges.size());
		for (j = 0; j < activeSmallEdgesInEdges.size(); ++j) {
			integrals1FormsEdges[i][j] = integrateWhitneyForm(*activeSmallEdgesInEdges[i].multiIndex, *activeSmallEdgesInEdges[i].nodeIndices,
				*activeSmallEdgesInEdges[j].multiIndex, *activeSmallEdgesInEdges[j].nodeIndices);
		}
	}

	//1-forms in faces
	faceEdgeValuesEdge0.resize(edgePermutations.size());
	faceEdgeValuesEdge1.resize(edgePermutations.size());
	faceEdgeValuesEdge2.resize(edgePermutations.size());
	for (i = 0; i < edgePermutations.size(); ++i) {
		faceEdgeValuesEdge0[i].resize(activeSmallEdgesInFaces.size());
		faceEdgeValuesEdge1[i].resize(activeSmallEdgesInFaces.size());
		faceEdgeValuesEdge2[i].resize(activeSmallEdgesInFaces.size());
		for (j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
			faceEdgeValuesEdge0[i][j].toVectorN(activeSmallEdgesInEdges.size());
			faceEdgeValuesEdge1[i][j].toVectorN(activeSmallEdgesInEdges.size());
			faceEdgeValuesEdge2[i][j].toVectorN(activeSmallEdgesInEdges.size());
		}
	}
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		for (j = 0; j < activeSmallEdgesInEdges.size(); ++j) {
			for (k = 0; k < edgePermutations.size(); ++k) {
				faceEdgeValuesEdge0[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInFaces[i].multiIndex, *activeSmallEdgesInFaces[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 2, 0), insert(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 2, 0));
				faceEdgeValuesEdge1[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInFaces[i].multiIndex, *activeSmallEdgesInFaces[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 1, 0), insert(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 1, 0));
				faceEdgeValuesEdge2[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInFaces[i].multiIndex, *activeSmallEdgesInFaces[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 0, 0), insert(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 0, 0));
			}
		}
		integrals1FormsFaces[i].resize(activeSmallEdgesInFaces.size());
		for (j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
			integrals1FormsFaces[i][j] = integrateWhitneyForm(*activeSmallEdgesInFaces[i].multiIndex, *activeSmallEdgesInFaces[i].nodeIndices,
				*activeSmallEdgesInFaces[j].multiIndex, *activeSmallEdgesInFaces[j].nodeIndices);
		}
	}

	//2-forms
	Buffer<Buffer<double>> integrals2FormsFaces(activeSmallFacesInFaces.size());

	//2-forms in faces
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		integrals2FormsFaces[i].resize(activeSmallFacesInFaces.size());
		for (j = 0; j < activeSmallFacesInFaces.size(); ++j) {
			integrals2FormsFaces[i][j] = integrateWhitneyForm(*activeSmallFacesInFaces[i].multiIndex, *activeSmallFacesInFaces[i].nodeIndices,
				*activeSmallFacesInFaces[j].multiIndex, *activeSmallFacesInFaces[j].nodeIndices);
		}
	}

	//precompute LU decompositions
	double tol = 1e-15;
	decomposeLUP(integrals0FormsEdges, matrix0FormsEdges_p, tol);
	decomposeLUP(integrals0FormsFaces, matrix0FormsFaces_p, tol);
	decomposeLUP(integrals1FormsEdges, matrix1FormsEdges_p, tol);
	decomposeLUP(integrals1FormsFaces, matrix1FormsFaces_p, tol);
	decomposeLUP(integrals2FormsFaces, matrix2FormsFaces_p, tol);
	//convert to MatrixN
	convertMatrix(integrals0FormsEdges, matrix0FormsEdges);
	convertMatrix(integrals0FormsFaces, matrix0FormsFaces);
	convertMatrix(integrals1FormsEdges, matrix1FormsEdges);
	convertMatrix(integrals1FormsFaces, matrix1FormsFaces);
	convertMatrix(integrals2FormsFaces, matrix2FormsFaces);
}

//save the matrices that have been precomputed
void SmallSimplexPartition2D::saveMatrices() const {

	if (order <= 3) //not all the matrices exist, and those who do are easily formed
		return;

	Text path;

	path << "Files/2D/order" << order << "/matrix0FormsEdges.dat";
	saveMatrix(matrix0FormsEdges, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix0FormsFaces.dat";
	saveMatrix(matrix0FormsFaces, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix1FormsEdges.dat";
	saveMatrix(matrix1FormsEdges, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix1FormsFaces.dat";
	saveMatrix(matrix1FormsFaces, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix2FormsFaces.dat";
	saveMatrix(matrix2FormsFaces, path.str());
	path.clear();

	path << "Files/2D/order" << order << "/matrix0FormsEdges_p.dat";
	saveBuffer(matrix0FormsEdges_p, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix0FormsFaces_p.dat";
	saveBuffer(matrix0FormsFaces_p, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix1FormsEdges_p.dat";
	saveBuffer(matrix1FormsEdges_p, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix1FormsFaces_p.dat";
	saveBuffer(matrix1FormsFaces_p, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix2FormsFaces_p.dat";
	saveBuffer(matrix2FormsFaces_p, path.str());
	path.clear();

	path << "Files/2D/order" << order << "/edgeValuesNode0.dat";
	saveBuffer(edgeValuesNode0.val, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/edgeValuesNode1.dat";
	saveBuffer(edgeValuesNode1.val, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceValuesNode0.dat";
	saveBuffer(faceValuesNode0.val, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceValuesNode1.dat";
	saveBuffer(faceValuesNode1.val, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceValuesNode2.dat";
	saveBuffer(faceValuesNode2.val, path.str());
	path.clear();

	path << "Files/2D/order" << order << "/faceNodeValuesEdge0Permutation.dat";
	saveMatrix(faceNodeValuesEdge0, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceNodeValuesEdge1Permutation.dat";
	saveMatrix(faceNodeValuesEdge1, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceNodeValuesEdge2Permutation.dat";
	saveMatrix(faceNodeValuesEdge2, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceEdgeValuesEdge0Permutation.dat";
	saveMatrix(faceEdgeValuesEdge0, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceEdgeValuesEdge1Permutation.dat";
	saveMatrix(faceEdgeValuesEdge1, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceEdgeValuesEdge2Permutation.dat";
	saveMatrix(faceEdgeValuesEdge2, path.str());
	path.clear();
}

//load the matrices that have been precomputed
void SmallSimplexPartition2D::loadMatrices() {

	if (order <= 3) { //not all the matrices exist, and those who do are easily formed
		formMatrices();
		return;
	}

	Text path;

	path << "Files/2D/order" << order << "/matrix0FormsEdges.dat";
	loadMatrix(matrix0FormsEdges, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix0FormsFaces.dat";
	loadMatrix(matrix0FormsFaces, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix1FormsEdges.dat";
	loadMatrix(matrix1FormsEdges, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix1FormsFaces.dat";
	loadMatrix(matrix1FormsFaces, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix2FormsFaces.dat";
	loadMatrix(matrix2FormsFaces, path.str());
	path.clear();

	path << "Files/2D/order" << order << "/matrix0FormsEdges_p.dat";
	loadBuffer(matrix0FormsEdges_p, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix0FormsFaces_p.dat";
	loadBuffer(matrix0FormsFaces_p, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix1FormsEdges_p.dat";
	loadBuffer(matrix1FormsEdges_p, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix1FormsFaces_p.dat";
	loadBuffer(matrix1FormsFaces_p, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/matrix2FormsFaces_p.dat";
	loadBuffer(matrix2FormsFaces_p, path.str());
	path.clear();

	path << "Files/2D/order" << order << "/edgeValuesNode0.dat";
	loadBuffer(edgeValuesNode0.val, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/edgeValuesNode1.dat";
	loadBuffer(edgeValuesNode1.val, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceValuesNode0.dat";
	loadBuffer(faceValuesNode0.val, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceValuesNode1.dat";
	loadBuffer(faceValuesNode1.val, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceValuesNode2.dat";
	loadBuffer(faceValuesNode2.val, path.str());
	path.clear();

	path << "Files/2D/order" << order << "/faceNodeValuesEdge0Permutation.dat";
	loadMatrix(faceNodeValuesEdge0, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceNodeValuesEdge1Permutation.dat";
	loadMatrix(faceNodeValuesEdge1, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceNodeValuesEdge2Permutation.dat";
	loadMatrix(faceNodeValuesEdge2, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceEdgeValuesEdge0Permutation.dat";
	loadMatrix(faceEdgeValuesEdge0, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceEdgeValuesEdge1Permutation.dat";
	loadMatrix(faceEdgeValuesEdge1, path.str());
	path.clear();
	path << "Files/2D/order" << order << "/faceEdgeValuesEdge2Permutation.dat";
	loadMatrix(faceEdgeValuesEdge2, path.str());
	path.clear();
}

/*
Save the integrals required in building to Hodge matrix for 1-forms for a specific element type (which should be indicated by foldername).
*/
void SmallSimplexPartition2D::saveHodgeIntegrals1Forms(std::string foldername, bool circumcentric, const Buffer<VectorN>& hodgeIntegralsFaceEdges,
	const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge0Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge1Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge2Edges, int ref) const {
	std::string dual = circumcentric ? "circumcentric" : "barycentric";
	std::string refElement = (ref < 0 ? "" : "Type" + std::to_string(ref));
	Text path;
	path << "Files/2D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsFaceEdges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsFaceEdges, path.str());
	path.clear();
	path << "Files/2D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge0Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge0Edges, path.str());
	path.clear();
	path << "Files/2D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge1Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge1Edges, path.str());
	path.clear();
	path << "Files/2D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge2Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge2Edges, path.str());
	path.clear();
}

/*
Load the integrals required in building to Hodge matrix for 1-forms for a specific element type (and have been stored in a folder indicated by foldername).
*/
void SmallSimplexPartition2D::loadHodgeIntegrals1Forms(std::string foldername, bool circumcentric, Buffer<VectorN>& hodgeIntegralsFaceEdges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge0Edges,
	Buffer<Buffer<VectorN>>& hodgeIntegralsEdge1Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge2Edges, int ref) const {
	std::string dual = circumcentric ? "circumcentric" : "barycentric";
	std::string refElement = (ref < 0 ? "" : "Type" + std::to_string(ref));
	Text path;
	path << "Files/2D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsFaceEdges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsFaceEdges, path.str());
	path.clear();
	path << "Files/2D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge0Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge0Edges, path.str());
	path.clear();
	path << "Files/2D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge1Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge1Edges, path.str());
	path.clear();
	path << "Files/2D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge2Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge2Edges, path.str());
	path.clear();
}