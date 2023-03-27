#include "SmallCubePartition2D.hpp"
#include "WhitneyForm.hpp"
#include "NumericalIntegration.hpp"
#include "LinearAlgebraFunctions.hpp"

using namespace gfd;

SmallCubePartition2D::SmallCubePartition2D(uint order) : order{ order } {
	initialiseSmallCells();
	formMatrices();
}

void SmallCubePartition2D::initialiseSmallCells() {
	//small edges in edges
	smallEdgesInEdges.reserve(order);
	for (uint xi = 0; xi < order; ++xi) {
		mi_t mi{ {xi} };
		smallEdgesInEdges.push_back({ mi, 0 });
	}
	//small edges in faces
	smallEdgesInFaces.reserve(2 * order * (order - 1));
	for (uint yi = 1; yi < order; ++yi) {
		for (uint xi = 0; xi < order; ++xi) {
			smallEdgesInFaces.push_back({ {xi, yi}, 0 });
		}
	}
	for (uint xi = 1; xi < order; ++xi) {
		for (uint yi = 0; yi < order; ++yi) {
			smallEdgesInFaces.push_back({ {xi, yi}, 1 });
		}
	}
}

void SmallCubePartition2D::refineMesh(const BuilderMesh& mesh_old, BuilderMesh& mesh) {
	//SmallCubePartition2D is now associated with the old mesh (mesh_old) and its refinement (mesh)
	mesh_old_ptr = &mesh_old;
	mesh_ptr = &mesh;

	//in the lowest order case, the mesh remains the same
	if (order == 1) {
		mesh.createCopy(mesh_old);
		return;
	}

	//first add the existing nodes
	for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
		mesh.addNode(mesh_old.getNodePosition(i));
	}

	//refine edges
	for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
		const Buffer<uint>& nodes = mesh_old.getEdgeNodes(i);
		Vector2 p0 = mesh.getNodePosition2(nodes[0]);
		Vector2 p1 = mesh.getNodePosition2(nodes[1]);
		//small nodes on edge i
		for (uint xi = 1; xi < order; ++xi) {
			double weight_x = (double)xi / (double)order;
			mesh.addNode(Vector4((1 - weight_x) * p0 + weight_x * p1, 0.0, 0.0));
		}
		//small edges in edge i
		mesh.addEdge(nodes[0], smallNodesEdgeList(i, 0));
		for (uint xi = 1; xi < order - 1; ++xi) {
			mesh.addEdge(smallNodesEdgeList(i, xi - 1), smallNodesEdgeList(i, xi));
		}
		mesh.addEdge(smallNodesEdgeList(i, order - 2), nodes[1]);
	}

	//refine faces
	for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
		const Buffer<uint> nodes = getQuadrilateralNodes(i, mesh_old);
		Vector2 p0 = mesh.getNodePosition2(nodes[0]);
		Vector2 p1 = mesh.getNodePosition2(nodes[1]);
		Vector2 p2 = mesh.getNodePosition2(nodes[2]);
		Vector2 p3 = mesh.getNodePosition2(nodes[3]);
		Vector2 xdiff = p1 - p0;
		Vector2 ydiff = p3 - p0;
		//small nodes in face i
		for (uint yi = 1; yi < order; ++yi) {
			double weight_y = (double)yi / (double)order;
			for (uint xi = 1; xi < order; ++xi) {
				double weight_x = (double)xi / (double)order;
				mesh.addNode(Vector4(p0 + weight_x * xdiff + weight_y * ydiff, 0.0, 0.0));
			}
		}
		//small edges in face i
		uint edge03 = mesh_old.findEdge(nodes[0], nodes[3]);
		uint edge12 = mesh_old.findEdge(nodes[1], nodes[2]);
		uint edge01 = mesh_old.findEdge(nodes[0], nodes[1]);
		uint edge23 = mesh_old.findEdge(nodes[2], nodes[3]);
		for (uint yi = 1; yi < order; ++yi) {
			double weight_y = (double)yi / (double)order;
			uint prev = getSmallNodeIndex(p0 + weight_y * ydiff, edge03, 1);
			for (uint xi = 1; xi < order; ++xi) {
				double weight_x = (double)xi / (double)order;
				uint curr = getSmallNodeIndex(p0 + weight_x * xdiff + weight_y * ydiff, i, 2);
				mesh.addEdge(prev, curr);
				prev = curr;
			}
			mesh.addEdge(prev, getSmallNodeIndex(p1 + weight_y * ydiff, edge12, 1));
		}
		for (uint xi = 1; xi < order; ++xi) {
			double weight_x = (double)xi / (double)order;
			uint prev = getSmallNodeIndex(p0 + weight_x * xdiff, edge01, 1);
			for (uint yi = 1; yi < order; ++yi) {
				double weight_y = (double)yi / (double)order;
				uint curr = getSmallNodeIndex(p0 + weight_x * xdiff + weight_y * ydiff, i, 2);
				mesh.addEdge(prev, curr);
				prev = curr;
			}
			mesh.addEdge(prev, getSmallNodeIndex(p3 + weight_x * xdiff, edge23, 1));
		}
		//small faces in face i
		for (uint yi = 0; yi < order; ++yi) {
			double weight_y_start = (double)yi / (double)order;
			double weight_y_end = (double)(yi + 1) / (double)order;
			for (uint xi = 0; xi < order; ++xi) {
				double weight_x_start = (double)xi / (double)order;
				double weight_x_end = (double)(xi + 1) / (double)order;
				uint node0, node1, node2, node3;
				Vector2 pos0 = p0 + weight_x_start * xdiff + weight_y_start * ydiff;
				Vector2 pos1 = p0 + weight_x_end * xdiff + weight_y_start * ydiff;
				Vector2 pos2 = p0 + weight_x_end * xdiff + weight_y_end * ydiff;
				Vector2 pos3 = p0 + weight_x_start * xdiff + weight_y_end * ydiff;
				if (xi == 0 && yi == 0) {
					node0 = nodes[0];
				}
				else if (xi == 0) {
					node0 = getSmallNodeIndex(pos0, edge03, 1);
				}
				else if (yi == 0) {
					node0 = getSmallNodeIndex(pos0, edge01, 1);
				}
				else {
					node0 = getSmallNodeIndex(pos0, i, 2);
				}
				if (xi == order - 1 && yi == 0) {
					node1 = nodes[1];
				}
				else if (xi == order - 1) {
					node1 = getSmallNodeIndex(pos1, edge12, 1);
				}
				else if (yi == 0) {
					node1 = getSmallNodeIndex(pos1, edge01, 1);
				}
				else {
					node1 = getSmallNodeIndex(pos1, i, 2);
				}
				if (xi == order - 1 && yi == order - 1) {
					node2 = nodes[2];
				}
				else if (xi == order - 1) {
					node2 = getSmallNodeIndex(pos2, edge12, 1);
				}
				else if (yi == order - 1) {
					node2 = getSmallNodeIndex(pos2, edge23, 1);
				}
				else {
					node2 = getSmallNodeIndex(pos2, i, 2);
				}
				if (xi == 0 && yi == order - 1) {
					node3 = nodes[3];
				}
				else if (xi == 0) {
					node3 = getSmallNodeIndex(pos3, edge03, 1);
				}
				else if (yi == order - 1) {
					node3 = getSmallNodeIndex(pos3, edge23, 1);
				}
				else {
					node3 = getSmallNodeIndex(pos3, i, 2);
				}
				Buffer<uint> edges(4);
				edges[0] = mesh.findEdge(node0, node1);
				edges[1] = mesh.findEdge(node1, node2);
				edges[2] = mesh.findEdge(node3, node2);
				edges[3] = mesh.findEdge(node0, node3);
				mesh.addFace(edges);
			}
		}
	}
}

//solves the coefficients of the (higher order) cubical 1-form (interpolant of discreteForm)
void SmallCubePartition2D::solve1FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallCubePartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;

	edgeCoefficients.resize(mesh_old.getEdgeSize()); //coefficients for small edges that are on big edges
	for (uint i = 0; i < edgeCoefficients.size(); ++i) {
		VectorN edgeValues(smallEdgesInEdges.size(), 0.0);
		for (uint j = 0; j < smallEdgesInEdges.size(); ++j) {
			edgeValues[j] = discreteForm[smallEdgesEdgeList(i, j)];
		}
		solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeValues, edgeCoefficients[i]);
	}

	faceCoefficients.resize(mesh_old.getFaceSize()); //coefficients for small edges in the interior
	for (uint i = 0; i < faceCoefficients.size(); ++i) {
		const Buffer<uint> faceNodes = getQuadrilateralNodes(i, mesh_old);
		Buffer<uint> faceEdges(4);
		faceEdges[0] = mesh_old.findEdge(faceNodes[0], faceNodes[1]);
		faceEdges[1] = mesh_old.findEdge(faceNodes[1], faceNodes[2]);
		faceEdges[2] = mesh_old.findEdge(faceNodes[2], faceNodes[3]);
		faceEdges[3] = mesh_old.findEdge(faceNodes[0], faceNodes[3]);
		VectorN faceValues(smallEdgesInFaces.size(), 0.0);
		for (uint j = 0; j < smallEdgesInFaces.size(); ++j) {
			for (uint k = 0; k < faceEdgeValues1Forms.size(); ++k) {
				faceValues[j] -= faceEdgeValues1Forms[k][j].dot(edgeCoefficients[faceEdges[k]]);
			}
			faceValues[j] += discreteForm[smallEdgesFaceList(i, j)];
		}
		solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues, faceCoefficients[i]);
	}
}

//computes the value of the (higher order) cubical 1-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known
Vector2 SmallCubePartition2D::evaluate1FormWithCoefficients(const Vector2& evaluationPoint, const Buffer<VectorN>& edgeCoefficients, const Buffer<VectorN>& faceCoefficients, uint element) const {
	if (!mesh_old_ptr)
		return Vector2(0.0, 0.0); //SmallCubePartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;

	const Buffer<uint> nodes = getQuadrilateralNodes(element, mesh_old);
	Buffer<uint> edges(4);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[1], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[2], nodes[3]);
	edges[3] = mesh_old.findEdge(nodes[0], nodes[3]);

	Vector2 result(0.0, 0.0);
	Vector2 p0 = mesh_old.getNodePosition2(nodes[0]);
	Vector2 p1 = mesh_old.getNodePosition2(nodes[2]);
	double lenx = p1.x - p0.x;
	double leny = p1.y - p0.y;
	double x = (evaluationPoint.x - p0.x) / lenx;
	double y = (evaluationPoint.y - p0.y) / leny;
	double x_d = (p1.x - evaluationPoint.x) / lenx;
	double y_d = (p1.y - evaluationPoint.y) / leny;
	double sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge0
		sum += edgeCoefficients[edges[0]][i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(y_d, order);
	result.x += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge1
		sum += edgeCoefficients[edges[1]][i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x, order);
	result.y += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge2
		sum += edgeCoefficients[edges[2]][i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(y, order);
	result.x += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge3
		sum += edgeCoefficients[edges[3]][i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x_d, order);
	result.y += sum;

	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) { //face
		const SmallEdge& se = smallEdgesInFaces[i];
		if (se.edge == 0) {
			result.x += faceCoefficients[element][i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]);
		}
		else {
			result.y += faceCoefficients[element][i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - 1 - se.mi[1]) * std::pow(y, se.mi[1]);
		}
	}
	result.x /= lenx;
	result.y /= leny;
	return result;
}

/*
Computes the discrete Hodge for higher order cubical 1-forms to a sparse matrix (of type Buffer<std::unordered_map<uint, double>>).
Integrals over dual edges are computed once in a single element, so all elements are assumed to have the same shape.
*/
void SmallCubePartition2D::formHodgeMatrix1Forms(Buffer<std::unordered_map<uint, double>>& star, bool minkowskiMetric) const {
	if (!mesh_old_ptr)
		return; //SmallCubePartition2D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	const BuilderMesh& mesh = *mesh_ptr;
	star.resize(mesh.getEdgeSize());

	//perform computations in face 0
	uint refElement = 0;
	const Buffer<uint> refElementNodes = getQuadrilateralNodes(0, mesh_old);
	Buffer<Vector2> refElementPos(4);
	for (uint i = 0; i < refElementPos.size(); ++i)
		refElementPos[i] = mesh_old.getNodePosition2(refElementNodes[i]);
	Buffer<uint> refElementEdges(4);
	refElementEdges[0] = mesh_old.findEdge(refElementNodes[0], refElementNodes[1]);
	refElementEdges[1] = mesh_old.findEdge(refElementNodes[1], refElementNodes[2]);
	refElementEdges[2] = mesh_old.findEdge(refElementNodes[2], refElementNodes[3]);
	refElementEdges[3] = mesh_old.findEdge(refElementNodes[0], refElementNodes[3]);
	Vector2 refElementEdgeX = refElementPos[1] - refElementPos[0];
	Vector2 refElementEdgeY = refElementPos[3] - refElementPos[0];
	double refElementLenX = refElementEdgeX.len() / order;
	double refElementLenY = refElementEdgeY.len() / order;

	//form the list of dual edges
	uint dualEdgeCount = 4 * order + 2 * order * (order - 1);
	struct DualEdgePart
	{
		Vector2 p0;
		Vector2 p1;
		uint edge;
		double length;
	};
	std::vector<DualEdgePart> dualEdgeList;
	dualEdgeList.resize(dualEdgeCount);
	uint deIndex = 0;
	//edge 0
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector2 p0 = mesh.getEdgeAverage2(smallEdgesEdgeList(refElementEdges[0], i));
		Vector2 p1 = p0 + refElementEdgeY / (2 * order);
		dualEdgeList[deIndex++] = { p0, p1, 0, refElementLenY / 2 };
	}
	//edge 1
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector2 p1 = mesh.getEdgeAverage2(smallEdgesEdgeList(refElementEdges[1], i));
		Vector2 p0 = p1 - refElementEdgeX / (2 * order);
		dualEdgeList[deIndex++] = { p0, p1, 1, refElementLenX / 2 };
	}
	//edge 2
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector2 p1 = mesh.getEdgeAverage2(smallEdgesEdgeList(refElementEdges[2], i));
		Vector2 p0 = p1 - refElementEdgeY / (2 * order);
		dualEdgeList[deIndex++] = { p0, p1, 0, refElementLenY / 2 };
	}
	//edge 3
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector2 p0 = mesh.getEdgeAverage2(smallEdgesEdgeList(refElementEdges[3], i));
		Vector2 p1 = p0 + refElementEdgeX / (2 * order);
		dualEdgeList[deIndex++] = { p0, p1, 1, refElementLenX / 2 };
	}
	//face
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		Vector2 edgeMidpoint = mesh.getEdgeAverage2(smallEdgesFaceList(0, i));
		uint edge = smallEdgesInFaces[i].edge;
		if (edge == 0) {
			Vector2 p0 = edgeMidpoint - refElementEdgeY / (2 * order);
			Vector2 p1 = edgeMidpoint + refElementEdgeY / (2 * order);
			dualEdgeList[deIndex++] = { p0, p1, 0, refElementLenY };
		}
		else {
			Vector2 p0 = edgeMidpoint - refElementEdgeX / (2 * order);
			Vector2 p1 = edgeMidpoint + refElementEdgeX / (2 * order);
			dualEdgeList[deIndex++] = { p0, p1, 1, refElementLenX };
		}
	}

	if (minkowskiMetric) {
		for (uint i = 0; i < deIndex; ++i) {
			if (dualEdgeList[i].edge == 1)
				dualEdgeList[i].length *= -1.0;
		}
	}

	//compute the integrals of basis functions
	Buffer<VectorN> hodgeIntegralsFaceEdges(dualEdgeCount);
	Buffer<Buffer<VectorN>> hodgeIntegralsEdgeEdges(4);
	for (uint i = 0; i < dualEdgeCount; ++i) {
		hodgeIntegralsFaceEdges[i].toVectorN(smallEdgesInFaces.size());
	}
	for (uint i = 0; i < 4; ++i) {
		hodgeIntegralsEdgeEdges[i].resize(dualEdgeCount);
		for (uint j = 0; j < dualEdgeCount; ++j) {
			hodgeIntegralsEdgeEdges[i][j].toVectorN(smallEdgesInEdges.size());
		}
	}

	//first integrate face basis functions
	VectorN faceCochain(smallEdgesInFaces.size(), 0.0);
	VectorN faceCoefficients(smallEdgesInFaces.size(), 0.0);
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		faceCochain[i] = 1.0;
		solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceCochain, faceCoefficients);
		for (uint j = 0; j < dualEdgeCount; ++j) {
			uint edge = dualEdgeList[j].edge;
			if (smallEdgesInFaces[i].edge != edge)
				continue;
			std::function<double(Vector2)> interpolant;
			if (edge == 0)
				interpolant = [&](Vector2 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients).x; };
			else
				interpolant = [&](Vector2 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients).y; };
			hodgeIntegralsFaceEdges[j][i] += dualEdgeList[j].length * integralAverage(interpolant, dualEdgeList[j].p0, dualEdgeList[j].p1);

		}
		faceCoefficients.val.fill(0.0);
		faceCochain[i] = 0.0;
	}

	//then integrate edge basis functions
	VectorN edgeCochain(smallEdgesInEdges.size(), 0.0);
	VectorN edgeCoefficients(smallEdgesInEdges.size(), 0.0);
	VectorN faceValues(smallEdgesInFaces.size(), 0.0);
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		edgeCochain[i] = 1.0;
		solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeCochain, edgeCoefficients);
		//edge0
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			faceValues[k] += faceEdgeValues1Forms[0][k].dot(edgeCoefficients);
		}
		for (uint j = 0; j < dualEdgeCount; ++j) {
			uint edge = dualEdgeList[j].edge;
			if (edge != 0)
				continue;
			std::function<double(Vector2)> interpolant = [&](Vector2 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 0).x; };
			hodgeIntegralsEdgeEdges[0][j][i] += dualEdgeList[j].length * integralAverage(interpolant, dualEdgeList[j].p0, dualEdgeList[j].p1);
			hodgeIntegralsEdgeEdges[0][j][i] -= faceValues.dot(hodgeIntegralsFaceEdges[j]);
		}
		faceValues.val.fill(0.0);
		//edge1
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			faceValues[k] += faceEdgeValues1Forms[1][k].dot(edgeCoefficients);
		}
		for (uint j = 0; j < dualEdgeCount; ++j) {
			uint edge = dualEdgeList[j].edge;
			if (edge != 1)
				continue;
			std::function<double(Vector2)> interpolant = [&](Vector2 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 1).y; };
			hodgeIntegralsEdgeEdges[1][j][i] += dualEdgeList[j].length * integralAverage(interpolant, dualEdgeList[j].p0, dualEdgeList[j].p1);
			hodgeIntegralsEdgeEdges[1][j][i] -= faceValues.dot(hodgeIntegralsFaceEdges[j]);
		}
		faceValues.val.fill(0.0);
		//edge2
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			faceValues[k] += faceEdgeValues1Forms[2][k].dot(edgeCoefficients);
		}
		for (uint j = 0; j < dualEdgeCount; ++j) {
			uint edge = dualEdgeList[j].edge;
			if (edge != 0)
				continue;
			std::function<double(Vector2)> interpolant = [&](Vector2 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 2).x; };
			hodgeIntegralsEdgeEdges[2][j][i] += dualEdgeList[j].length * integralAverage(interpolant, dualEdgeList[j].p0, dualEdgeList[j].p1);
			hodgeIntegralsEdgeEdges[2][j][i] -= faceValues.dot(hodgeIntegralsFaceEdges[j]);
		}
		faceValues.val.fill(0.0);
		//edge3
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			faceValues[k] += faceEdgeValues1Forms[3][k].dot(edgeCoefficients);
		}
		for (uint j = 0; j < dualEdgeCount; ++j) {
			uint edge = dualEdgeList[j].edge;
			if (edge != 1)
				continue;
			std::function<double(Vector2)> interpolant = [&](Vector2 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 3).y; };
			hodgeIntegralsEdgeEdges[3][j][i] += dualEdgeList[j].length * integralAverage(interpolant, dualEdgeList[j].p0, dualEdgeList[j].p1);
			hodgeIntegralsEdgeEdges[3][j][i] -= faceValues.dot(hodgeIntegralsFaceEdges[j]);
		}
		faceValues.val.fill(0.0);
		edgeCoefficients.val.fill(0.0);
		edgeCochain[i] = 0.0;
	}

	//form Hodge matrix from the precomputed integrals
	for (uint f = 0; f < mesh_old.getFaceSize(); ++f) {
		const Buffer<uint> nodes = getQuadrilateralNodes(f, mesh_old);
		Buffer<Vector2> pos(4);
		for (uint i = 0; i < pos.size(); ++i)
			pos[i] = mesh_old.getNodePosition2(nodes[i]);
		Buffer<uint> edges(4);
		edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
		edges[1] = mesh_old.findEdge(nodes[1], nodes[2]);
		edges[2] = mesh_old.findEdge(nodes[2], nodes[3]);
		edges[3] = mesh_old.findEdge(nodes[0], nodes[3]);
		Vector2 edgeX = pos[1] - pos[0];
		Vector2 edgeY = pos[3] - pos[0];
		double lenX = edgeX.len() / order;
		double lenY = edgeY.len() / order;

		//find the indices
		Buffer<uint> indices(dualEdgeCount);
		deIndex = 0;
		//edge 0
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[deIndex++] = smallEdgesEdgeList(edges[0], i);
		}
		//edge 1
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[deIndex++] = smallEdgesEdgeList(edges[1], i);
		}
		//edge 2
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[deIndex++] = smallEdgesEdgeList(edges[2], i);
		}
		//edge 3
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[deIndex++] = smallEdgesEdgeList(edges[3], i);
		}
		//face
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			indices[deIndex++] = smallEdgesFaceList(f, i);
		}

		for (uint i = 0; i < dualEdgeCount; ++i) { //edge basis functions
			for (uint j = 0; j < smallEdgesInEdges.size(); ++j) {
				star[indices[i]][smallEdgesEdgeList(edges[0], j)] += hodgeIntegralsEdgeEdges[0][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[1], j)] += hodgeIntegralsEdgeEdges[1][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[2], j)] += hodgeIntegralsEdgeEdges[2][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[3], j)] += hodgeIntegralsEdgeEdges[3][i][j];
			}
		}
		for (uint i = 0; i < dualEdgeCount; ++i) { //face basis functions
			for (uint j = 0; j < smallEdgesInFaces.size(); ++j) {
				star[indices[i]][smallEdgesFaceList(f, j)] += hodgeIntegralsFaceEdges[i][j];
			}
		}
	}
}

//The following functions return the index of the small cell j of the big cell i.

uint SmallCubePartition2D::smallNodesEdgeList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + i * (order - 1) + j;
}
uint SmallCubePartition2D::smallEdgesEdgeList(uint i, uint j) const {
	return i * order + j;
}
uint SmallCubePartition2D::smallNodesFaceList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + mesh_old_ptr->getEdgeSize() * (order - 1) + i * (order - 1) * (order - 1) + j;
}
uint SmallCubePartition2D::smallEdgesFaceList(uint i, uint j) const {
	return mesh_old_ptr->getEdgeSize() * order + i * 2 * order * (order - 1) + j;
}
uint SmallCubePartition2D::smallFacesFaceList(uint i, uint j) const {
	return i * order * order + j;
}

//returns the index of the small node located at pos in the element of given dimension
uint SmallCubePartition2D::getSmallNodeIndex(const Vector2& pos, uint element, uint dim) {
	const BuilderMesh& mesh = *mesh_ptr;
	if (dim == 1) {
		for (uint i = 0; i < order - 1; ++i) {
			uint node = smallNodesEdgeList(element, i);
			if ((mesh.getNodePosition2(node) - pos).lensq() < 1e-24)
				return node;
		}
	}
	else if (dim == 2) {
		for (uint i = 0; i < (order - 1) * (order - 1); ++i) {
			uint node = smallNodesFaceList(element, i);
			if ((mesh.getNodePosition2(node) - pos).lensq() < 1e-24)
				return node;
		}
	}
	return NONE;
}

//computes the value resulting from the coefficients on the given edge in the element whose nodes are given in pos
Vector2 SmallCubePartition2D::evaluate1FormWithEdgeCoefficients(const Vector2& evaluationPoint, const Buffer<Vector2>& pos, const VectorN& edgeCoefficients, uint edge) const {
	Vector2 result(0.0, 0.0);
	double lenx = pos[1].x - pos[0].x;
	double leny = pos[3].y - pos[0].y;
	double x = (evaluationPoint.x - pos[0].x) / lenx;
	double y = (evaluationPoint.y - pos[0].y) / leny;
	double x_d = (pos[1].x - evaluationPoint.x) / lenx;
	double y_d = (pos[3].y - evaluationPoint.y) / leny;
	double sum = 0.0;
	if (edge == 0) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(y_d, order);
		result.x += sum;
		sum = 0.0;
	}
	else if (edge == 1) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x, order);
		result.y += sum;
		sum = 0.0;
	}
	else if (edge == 2) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(y, order);
		result.x += sum;
		sum = 0.0;
	}
	else if (edge == 3) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x_d, order);
		result.y += sum;
		sum = 0.0;
	}
	result.x /= lenx;
	result.y /= leny;
	return result;
}

//computes the value resulting from the coefficients in the element whose nodes are given in pos
Vector2 SmallCubePartition2D::evaluate1FormWithFaceCoefficients(const Vector2& evaluationPoint, const Buffer<Vector2>& pos, const VectorN& faceCoefficients) const {
	Vector2 result(0.0, 0.0);
	double lenx = pos[1].x - pos[0].x;
	double leny = pos[3].y - pos[0].y;
	double x = (evaluationPoint.x - pos[0].x) / lenx;
	double y = (evaluationPoint.y - pos[0].y) / leny;
	double x_d = (pos[1].x - evaluationPoint.x) / lenx;
	double y_d = (pos[3].y - evaluationPoint.y) / leny;
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		const SmallEdge& se = smallEdgesInFaces[i];
		if (se.edge == 0) {
			result.x += faceCoefficients[i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]);
		}
		else {
			result.y += faceCoefficients[i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - 1 - se.mi[1]) * std::pow(y, se.mi[1]);
		}
	}
	result.x /= lenx;
	result.y /= leny;
	return result;
}

//precomputes the matrices required in evaluating higher order cubical forms
void SmallCubePartition2D::formMatrices() {
	integrals1Forms.resize(smallEdgesInEdges.size());
	Buffer<Buffer<double>> integrals1FormsEdges(smallEdgesInEdges.size());
	Buffer<Buffer<double>> integrals1FormsFaces(smallEdgesInFaces.size());

	//1-forms in edges
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		integrals1FormsEdges[i].resize(smallEdgesInEdges.size());
		integrals1Forms[i].resize(smallEdgesInEdges.size());
		const SmallEdge& se_i = smallEdgesInEdges[i]; //small edge i is the domain of integration
		for (uint j = 0; j < smallEdgesInEdges.size(); ++j) {
			const SmallEdge& se_j = smallEdgesInEdges[j]; //the form corresponding to small edge j is the integrand
			double x0 = (double)(se_i.mi[0]) / (double)(order);
			double x1 = (double)(se_i.mi[0] + 1) / (double)(order);
			integrals1Forms[i][j] = integralAverage([=](double x) -> double {return std::pow(1 - x, order - 1 - se_j.mi[0]) * std::pow(x, se_j.mi[0]); }, x0, x1) / order;
			integrals1FormsEdges[i][j] = integrals1Forms[i][j];
		}
	}

	//1-forms in faces
	faceEdgeValues1Forms.resize(4);
	for (uint k = 0; k < 4; ++k) {
		faceEdgeValues1Forms[k].resize(smallEdgesInFaces.size());
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			faceEdgeValues1Forms[k][i].toVectorN(smallEdgesInEdges.size());
			const SmallEdge& se_i = smallEdgesInFaces[i]; //small edge i is the domain of integration
			for (uint j = 0; j < smallEdgesInEdges.size(); ++j) {
				const SmallEdge& se_j = smallEdgesInEdges[j]; //the form corresponding to small edge j (in edge k) is the integrand
				if ((se_i.edge == 0 && (k == 1 || k == 3)) || (se_i.edge == 1 && (k == 0 || k == 2))) {
					faceEdgeValues1Forms[k][i][j] = 0.0;
				}
				else {
					if (k == 0) {
						double y = (double)(se_i.mi[1]) / (double)(order);
						faceEdgeValues1Forms[k][i][j] = std::pow(1 - y, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
					}
					else if (k == 2) {
						double y = (double)(se_i.mi[1]) / (double)(order);
						faceEdgeValues1Forms[k][i][j] = std::pow(y, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
					}
					else if (k == 3) {
						double x = (double)(se_i.mi[0]) / (double)(order);
						faceEdgeValues1Forms[k][i][j] = std::pow(1 - x, order) * integrals1Forms[se_i.mi[1]][se_j.mi[0]];
					}
					else { //k == 1
						double x = (double)(se_i.mi[0]) / (double)(order);
						faceEdgeValues1Forms[k][i][j] = std::pow(x, order) * integrals1Forms[se_i.mi[1]][se_j.mi[0]];
					}
				}
			}
		}
	}
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		integrals1FormsFaces[i].resize(smallEdgesInFaces.size());
		const SmallEdge& se_i = smallEdgesInFaces[i]; //small edge i is the domain of integration
		for (uint j = 0; j < smallEdgesInFaces.size(); ++j) {
			const SmallEdge& se_j = smallEdgesInFaces[j]; //the form corresponding to small edge j is the integrand
			if (se_i.edge != se_j.edge)
				integrals1FormsFaces[i][j] = 0.0;
			else {
				if (se_i.edge == 0) { //dx-component
					double y = (double)(se_i.mi[1]) / (double)(order);
					integrals1FormsFaces[i][j] = std::pow(1 - y, order - se_j.mi[1]) * std::pow(y, se_j.mi[1]) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
				}
				else { //dy-component
					double x = (double)(se_i.mi[0]) / (double)(order);
					integrals1FormsFaces[i][j] = std::pow(1 - x, order - se_j.mi[0]) * std::pow(x, se_j.mi[0]) * integrals1Forms[se_i.mi[1]][se_j.mi[1]];
				}
			}
		}
	}

	//precompute LU decompositions
	double tol = 1e-15;
	decomposeLUP(integrals1FormsEdges, matrix1FormsEdges_p, tol);
	decomposeLUP(integrals1FormsFaces, matrix1FormsFaces_p, tol);
	//convert to MatrixN
	convertMatrix(integrals1FormsEdges, matrix1FormsEdges);
	convertMatrix(integrals1FormsFaces, matrix1FormsFaces);
}
