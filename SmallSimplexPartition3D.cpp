#include "SmallSimplexPartition3D.hpp"
#include "WhitneyForm.hpp"
#include "NumericalIntegration.hpp"
#include "LinearAlgebraFunctions.hpp"

using namespace gfd;

SmallSimplexPartition3D::SmallSimplexPartition3D(uint order, bool loadMatricesFromFile) : SmallSimplexPartition(order) {
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
void SmallSimplexPartition3D::initialiseMultiIndices() {
	//initialise the multi-indices in edges, faces, and bodies
	getMultiIndices(2, order - 1, edgeMultiIndices);
	getMultiIndices(3, order - 1, faceMultiIndices);
	getMultiIndices(4, order - 1, bodyMultiIndices);
	//initialise the hole indices in faces and bodies
	if (order >= 2) {
		getMultiIndices(3, order - 2, faceEdgeHoles);
		getMultiIndices(4, order - 2, bodyFaceHoles);
	}
	if(order >= 3)
		getMultiIndices(4, order - 3, bodyEdgeHoles);
}

//index the small simplices that are in the interior of a big simplex, choosing and omitting redundant ones
void SmallSimplexPartition3D::formActiveSmallSimplices() {
	
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

	//body as big simplex
	activeSmallBodiesInBodies.resize(bodyMultiIndices.size());
	for (i = 0; i < bodyMultiIndices.size(); ++i) { //small nodes and small bodies
		if (mapsToInterior(bodyMultiIndices[i], bodyNodes[0])) {
			activeSmallNodesInBodies.push_back({ &bodyMultiIndices[i], &bodyNodes[0] });
		}
		activeSmallBodiesInBodies[i].multiIndex = &bodyMultiIndices[i];
		activeSmallBodiesInBodies[i].nodeIndices = &bodyBody;
	}
	if (order >= 2) { //small faces
		for (i = 0; i < bodyFaceHoles.size(); ++i) {
			for (j = 0; j < bodyFaces.size(); ++j) {
				if (j == 0) { // omit small faces in the direction of face 0
					continue;
				}
				//if (i % bodyFaces.size() == j) { // alternative: omit a similar number of faces in each direction
				//	continue;
				//}
				activeSmallFacesInBodies.push_back({ getHoleFaceMultiIndex(bodyFaceHoles[i], j), &bodyFaces[j] });
			}
		}
	}
	if (order >= 3) { //small edges
		for (i = 0; i < bodyEdgeHoles.size(); ++i) { // in each inverted tetrahedron choose three of the six edges as active
			uint edge1, edge2, edge3;
			if (true) { // choose edges 1, 3, and 5
				edge1 = 1;
				edge2 = 3;
				edge3 = 5;
			}
			//if (i % 2 == 0) { // alternative: choose a similar number of edges in each direction
			//	edge1 = 0;
			//	edge2 = 3;
			//	edge3 = 5;
			//}
			//else {
			//	edge1 = 1;
			//	edge2 = 2;
			//	edge3 = 4;
			//}
			activeSmallEdgesInBodies.push_back({ getHoleEdgeMultiIndex(bodyEdgeHoles[i], edge1), &bodyEdges[edge1] });
			activeSmallEdgesInBodies.push_back({ getHoleEdgeMultiIndex(bodyEdgeHoles[i], edge2), &bodyEdges[edge2] });
			activeSmallEdgesInBodies.push_back({ getHoleEdgeMultiIndex(bodyEdgeHoles[i], edge3), &bodyEdges[edge3] });
		}
	}
}

//refine mesh and save the required information that implies indices of small simplices
void SmallSimplexPartition3D::refineMesh(const BuilderMesh& mesh_old, BuilderMesh& mesh) {

	uint i, j, k;
	//SmallSimplexPartition3D is now associated with the old mesh (mesh_old) and its refinement (mesh)
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
		const Buffer<uint> nodes = mesh_old.getFaceNodes(i);
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

	//refine bodies
	smallEdgesBodyOffsets.resize(activeSmallEdgesInBodies.size());
	smallFacesBodyOffsets.resize(activeSmallFacesInBodies.size());
	for (i = 0; i < mesh_old.getBodySize(); ++i) {
		const Buffer<uint> nodes = getBodyNodesCustom(i);
		//small nodes in body i
		for (j = 0; j < activeSmallNodesInBodies.size(); ++j) {
			mesh.addNode(getNodeImage(*activeSmallNodesInBodies[j].multiIndex, mesh_old, nodes, *activeSmallNodesInBodies[j].nodeIndices));
		}
		//small edges, small faces, and small bodies in body i
		for (j = 0; j < bodyMultiIndices.size(); ++j) {
			uint node0 = getSmallNodeIndex(bodyMultiIndices[j], 0, nodes);
			uint node1 = getSmallNodeIndex(bodyMultiIndices[j], 1, nodes);
			uint node2 = getSmallNodeIndex(bodyMultiIndices[j], 2, nodes);
			uint node3 = getSmallNodeIndex(bodyMultiIndices[j], 3, nodes);
			Buffer<uint> edges(6);
			edges[0] = mesh.addEdge(node0, node1);
			edges[1] = mesh.addEdge(node0, node2);
			edges[2] = mesh.addEdge(node0, node3);
			edges[3] = mesh.addEdge(node1, node2);
			edges[4] = mesh.addEdge(node1, node3);
			edges[5] = mesh.addEdge(node2, node3);
			if (i == 0) { //the information implying the indices of small simplices is collected in body 0
				for (k = 0; k < edges.size(); ++k) {
					if (mapsToInterior(bodyMultiIndices[j], bodyEdges[k])) {
						SmallSimplex smallEdge{ &bodyMultiIndices[j], &bodyEdges[k] };
						uint smallEdgeIndex = indexOf(activeSmallEdgesInBodies, smallEdge);
						if (smallEdgeIndex != NONE)
							smallEdgesBodyOffsets[smallEdgeIndex] = edges[k];
					}
				}
			}
			Buffer<uint> faces(4);
			faces[0] = mesh.addFace(findEdges(node0, node1, node2, mesh));
			faces[1] = mesh.addFace(findEdges(node0, node1, node3, mesh));
			faces[2] = mesh.addFace(findEdges(node0, node2, node3, mesh));
			faces[3] = mesh.addFace(findEdges(node1, node2, node3, mesh));
			if (i == 0) { //the information implying the indices of small simplices is collected in body 0
				for (k = 0; k < faces.size(); ++k) {
					if (mapsToInterior(bodyMultiIndices[j], bodyFaces[k])) {
						SmallSimplex smallFace{ &bodyMultiIndices[j], &bodyFaces[k] };
						uint smallFaceIndex = indexOf(activeSmallFacesInBodies, smallFace);
						if (smallFaceIndex != NONE)
							smallFacesBodyOffsets[smallFaceIndex] = faces[k];
					}
				}
			}
			mesh.addBody(faces);
		}
		//fill the octahedra in body i
		enum class ClosestOppositeEdges { //denotes which of the opposite edges are closest to each other in body i
			edges05, //edges 0 and 5 are closest opposite edges
			edges14,
			edges23,
			edges05and14, //edges 0 and 5 have same distance as edges 1 and 4, both pairs are closest
			edges05and23,
			edges14and23,
			edgesAll //all opposite edges have same distance
		};
		ClosestOppositeEdges closestOppositeEdges;
		uint edge0 = mesh_old.findEdge(nodes[0], nodes[1]);
		uint edge1 = mesh_old.findEdge(nodes[0], nodes[2]);
		uint edge2 = mesh_old.findEdge(nodes[0], nodes[3]);
		uint edge3 = mesh_old.findEdge(nodes[1], nodes[2]);
		uint edge4 = mesh_old.findEdge(nodes[1], nodes[3]);
		uint edge5 = mesh_old.findEdge(nodes[2], nodes[3]);
		double dist05 = (mesh_old.getEdgePosition(edge0) - mesh_old.getEdgePosition(edge5)).len();
		double dist14 = (mesh_old.getEdgePosition(edge1) - mesh_old.getEdgePosition(edge4)).len();
		double dist23 = (mesh_old.getEdgePosition(edge2) - mesh_old.getEdgePosition(edge3)).len();
		if (dist05 < dist14 && dist05 < dist23) //add tolerance if you prefer
			closestOppositeEdges = ClosestOppositeEdges::edges05;
		else if (dist14 < dist05 && dist14 < dist23)
			closestOppositeEdges = ClosestOppositeEdges::edges14;
		else /*if (dist23 < dist05&& dist23 < dist14)*/
			closestOppositeEdges = ClosestOppositeEdges::edges23; //currently always divides into four tetrahedra; remove these comments for the possibility of combining into pyramids or octahedra
		/*else if (dist05 < dist23)
			closestOppositeEdges = ClosestOppositeEdges::edges05and14;
		else if (dist05 < dist14)
			closestOppositeEdges = ClosestOppositeEdges::edges05and23;
		else if (dist14 < dist05)
			closestOppositeEdges = ClosestOppositeEdges::edges14and23;
		else
			closestOppositeEdges = ClosestOppositeEdges::edgesAll;*/
		for (j = 0; j < bodyFaceHoles.size(); ++j) {
			Buffer<uint> holeNodes(6);
			for (k = 0; k < 3; ++k) {
				mi_t multiIndex{ bodyFaceHoles[j] };
				multiIndex[k] += 1;
				holeNodes[k] = getSmallNodeIndex(multiIndex, (k + 1) % 3, nodes);
			}
			uint holeNode2 = holeNodes[1]; //
			holeNodes[1] = holeNodes[2]; //
			holeNodes[2] = holeNode2; //order such that holeNodes[i] is in edge parallel to edge i
			mi_t multiIndex{ bodyFaceHoles[j] };
			multiIndex[3] += 1;
			for (k = 0; k < 3; ++k) {
				holeNodes[3 + k] = getSmallNodeIndex(multiIndex, k, nodes);
			}
			mesh.addFace(findEdges(holeNodes[0], holeNodes[1], holeNodes[2], mesh));
			mesh.addFace(findEdges(holeNodes[0], holeNodes[3], holeNodes[4], mesh));
			mesh.addFace(findEdges(holeNodes[1], holeNodes[3], holeNodes[5], mesh));
			mesh.addFace(findEdges(holeNodes[2], holeNodes[4], holeNodes[5], mesh));
			switch (closestOppositeEdges) {
			case ClosestOppositeEdges::edges05:
				mesh.addEdge(holeNodes[0], holeNodes[5]);
				mesh.addFace(findEdges(holeNodes[0], holeNodes[5], holeNodes[1], mesh));
				mesh.addFace(findEdges(holeNodes[0], holeNodes[5], holeNodes[2], mesh));
				mesh.addFace(findEdges(holeNodes[0], holeNodes[5], holeNodes[3], mesh));
				mesh.addFace(findEdges(holeNodes[0], holeNodes[5], holeNodes[4], mesh));
				mesh.addBody(findFaces(holeNodes[0], holeNodes[5], holeNodes[1], holeNodes[2], mesh));
				mesh.addBody(findFaces(holeNodes[0], holeNodes[5], holeNodes[1], holeNodes[3], mesh));
				mesh.addBody(findFaces(holeNodes[0], holeNodes[5], holeNodes[2], holeNodes[4], mesh));
				mesh.addBody(findFaces(holeNodes[0], holeNodes[5], holeNodes[3], holeNodes[4], mesh));
				break;
			case ClosestOppositeEdges::edges14:
				mesh.addEdge(holeNodes[1], holeNodes[4]);
				mesh.addFace(findEdges(holeNodes[1], holeNodes[4], holeNodes[0], mesh));
				mesh.addFace(findEdges(holeNodes[1], holeNodes[4], holeNodes[2], mesh));
				mesh.addFace(findEdges(holeNodes[1], holeNodes[4], holeNodes[3], mesh));
				mesh.addFace(findEdges(holeNodes[1], holeNodes[4], holeNodes[5], mesh));
				mesh.addBody(findFaces(holeNodes[1], holeNodes[4], holeNodes[0], holeNodes[2], mesh));
				mesh.addBody(findFaces(holeNodes[1], holeNodes[4], holeNodes[0], holeNodes[3], mesh));
				mesh.addBody(findFaces(holeNodes[1], holeNodes[4], holeNodes[2], holeNodes[5], mesh));
				mesh.addBody(findFaces(holeNodes[1], holeNodes[4], holeNodes[3], holeNodes[5], mesh));
				break;
			case ClosestOppositeEdges::edges23:
				mesh.addEdge(holeNodes[2], holeNodes[3]);
				mesh.addFace(findEdges(holeNodes[2], holeNodes[3], holeNodes[0], mesh));
				mesh.addFace(findEdges(holeNodes[2], holeNodes[3], holeNodes[1], mesh));
				mesh.addFace(findEdges(holeNodes[2], holeNodes[3], holeNodes[4], mesh));
				mesh.addFace(findEdges(holeNodes[2], holeNodes[3], holeNodes[5], mesh));
				mesh.addBody(findFaces(holeNodes[2], holeNodes[3], holeNodes[0], holeNodes[1], mesh));
				mesh.addBody(findFaces(holeNodes[2], holeNodes[3], holeNodes[0], holeNodes[4], mesh));
				mesh.addBody(findFaces(holeNodes[2], holeNodes[3], holeNodes[1], holeNodes[5], mesh));
				mesh.addBody(findFaces(holeNodes[2], holeNodes[3], holeNodes[4], holeNodes[5], mesh));
				break;
			case ClosestOppositeEdges::edges05and14: {
				Buffer<uint> baseEdges(4);
				baseEdges[0] = mesh.findEdge(holeNodes[0], holeNodes[1]);
				baseEdges[1] = mesh.findEdge(holeNodes[0], holeNodes[4]);
				baseEdges[2] = mesh.findEdge(holeNodes[1], holeNodes[5]);
				baseEdges[3] = mesh.findEdge(holeNodes[4], holeNodes[5]);
				uint base = mesh.addFace(baseEdges);
				Buffer<uint> pyramidFaces1(5);
				pyramidFaces1[0] = base;
				pyramidFaces1[1] = mesh.findFace(findEdges(holeNodes[2], holeNodes[0], holeNodes[1], mesh));
				pyramidFaces1[2] = mesh.findFace(findEdges(holeNodes[2], holeNodes[0], holeNodes[4], mesh));
				pyramidFaces1[3] = mesh.findFace(findEdges(holeNodes[2], holeNodes[1], holeNodes[5], mesh));
				pyramidFaces1[4] = mesh.findFace(findEdges(holeNodes[2], holeNodes[4], holeNodes[5], mesh));
				mesh.addBody(pyramidFaces1);
				Buffer<uint> pyramidFaces2(5);
				pyramidFaces2[0] = base;
				pyramidFaces2[1] = mesh.findFace(findEdges(holeNodes[3], holeNodes[0], holeNodes[1], mesh));
				pyramidFaces2[2] = mesh.findFace(findEdges(holeNodes[3], holeNodes[0], holeNodes[4], mesh));
				pyramidFaces2[3] = mesh.findFace(findEdges(holeNodes[3], holeNodes[1], holeNodes[5], mesh));
				pyramidFaces2[4] = mesh.findFace(findEdges(holeNodes[3], holeNodes[4], holeNodes[5], mesh));
				mesh.addBody(pyramidFaces2);
				break;
			}
			case ClosestOppositeEdges::edges05and23: {
				Buffer<uint> baseEdges(4);
				baseEdges[0] = mesh.findEdge(holeNodes[0], holeNodes[2]);
				baseEdges[1] = mesh.findEdge(holeNodes[0], holeNodes[3]);
				baseEdges[2] = mesh.findEdge(holeNodes[2], holeNodes[5]);
				baseEdges[3] = mesh.findEdge(holeNodes[3], holeNodes[5]);
				uint base = mesh.addFace(baseEdges);
				Buffer<uint> pyramidFaces1(5);
				pyramidFaces1[0] = base;
				pyramidFaces1[1] = mesh.findFace(findEdges(holeNodes[1], holeNodes[0], holeNodes[2], mesh));
				pyramidFaces1[2] = mesh.findFace(findEdges(holeNodes[1], holeNodes[0], holeNodes[3], mesh));
				pyramidFaces1[3] = mesh.findFace(findEdges(holeNodes[1], holeNodes[2], holeNodes[5], mesh));
				pyramidFaces1[4] = mesh.findFace(findEdges(holeNodes[1], holeNodes[3], holeNodes[5], mesh));
				mesh.addBody(pyramidFaces1);
				Buffer<uint> pyramidFaces2(5);
				pyramidFaces2[0] = base;
				pyramidFaces2[1] = mesh.findFace(findEdges(holeNodes[2], holeNodes[0], holeNodes[2], mesh));
				pyramidFaces2[2] = mesh.findFace(findEdges(holeNodes[2], holeNodes[0], holeNodes[3], mesh));
				pyramidFaces2[3] = mesh.findFace(findEdges(holeNodes[2], holeNodes[2], holeNodes[5], mesh));
				pyramidFaces2[4] = mesh.findFace(findEdges(holeNodes[2], holeNodes[3], holeNodes[5], mesh));
				mesh.addBody(pyramidFaces2);
				break;
			}
			case ClosestOppositeEdges::edges14and23: {
				Buffer<uint> baseEdges(4);
				baseEdges[0] = mesh.findEdge(holeNodes[1], holeNodes[2]);
				baseEdges[1] = mesh.findEdge(holeNodes[1], holeNodes[3]);
				baseEdges[2] = mesh.findEdge(holeNodes[2], holeNodes[4]);
				baseEdges[3] = mesh.findEdge(holeNodes[3], holeNodes[4]);
				uint base = mesh.addFace(baseEdges);
				Buffer<uint> pyramidFaces1(5);
				pyramidFaces1[0] = base;
				pyramidFaces1[1] = mesh.findFace(findEdges(holeNodes[0], holeNodes[1], holeNodes[2], mesh));
				pyramidFaces1[2] = mesh.findFace(findEdges(holeNodes[0], holeNodes[1], holeNodes[3], mesh));
				pyramidFaces1[3] = mesh.findFace(findEdges(holeNodes[0], holeNodes[2], holeNodes[4], mesh));
				pyramidFaces1[4] = mesh.findFace(findEdges(holeNodes[0], holeNodes[3], holeNodes[4], mesh));
				mesh.addBody(pyramidFaces1);
				Buffer<uint> pyramidFaces2(5);
				pyramidFaces2[0] = base;
				pyramidFaces2[1] = mesh.findFace(findEdges(holeNodes[5], holeNodes[1], holeNodes[2], mesh));
				pyramidFaces2[2] = mesh.findFace(findEdges(holeNodes[5], holeNodes[1], holeNodes[3], mesh));
				pyramidFaces2[3] = mesh.findFace(findEdges(holeNodes[5], holeNodes[2], holeNodes[4], mesh));
				pyramidFaces2[4] = mesh.findFace(findEdges(holeNodes[5], holeNodes[3], holeNodes[4], mesh));
				mesh.addBody(pyramidFaces2);
				break;
			}
			default: //edgesAll
				Buffer<uint> octahedronFaces(8);
				octahedronFaces[0] = mesh.findFace(findEdges(holeNodes[0], holeNodes[1], holeNodes[2], mesh));
				octahedronFaces[1] = mesh.findFace(findEdges(holeNodes[0], holeNodes[1], holeNodes[3], mesh));
				octahedronFaces[2] = mesh.findFace(findEdges(holeNodes[0], holeNodes[2], holeNodes[4], mesh));
				octahedronFaces[3] = mesh.findFace(findEdges(holeNodes[0], holeNodes[3], holeNodes[4], mesh));
				octahedronFaces[4] = mesh.findFace(findEdges(holeNodes[1], holeNodes[2], holeNodes[5], mesh));
				octahedronFaces[5] = mesh.findFace(findEdges(holeNodes[1], holeNodes[3], holeNodes[5], mesh));
				octahedronFaces[6] = mesh.findFace(findEdges(holeNodes[2], holeNodes[4], holeNodes[5], mesh));
				octahedronFaces[7] = mesh.findFace(findEdges(holeNodes[3], holeNodes[4], holeNodes[5], mesh));
				mesh.addBody(octahedronFaces);
				break;
			}
		}
		//fill inverted tetrahedra in body i
		for (j = 0; j < bodyEdgeHoles.size(); ++j) {
			Buffer<uint> holeNodes(4);
			for (k = 0; k < holeNodes.size(); ++k) {
				mi_t multiIndex{ bodyEdgeHoles[j] };
				multiIndex[k] += 1;
				multiIndex[(k + 1) % holeNodes.size()] += 1;
				holeNodes[k] = getSmallNodeIndex(multiIndex, (k + 2) % holeNodes.size(), nodes);
			}
			mesh.addBody(findFaces(holeNodes[0], holeNodes[1], holeNodes[2], holeNodes[3], mesh));
		}
	}
}

/*
Computes the value of the (higher order) Whitney 0-form (interpolant of discreteForm) at evaluationPoint.
When evaluating many times, it is more efficient to first solve the coefficients with solve0FormCoefficients and then evaluate with evaluate0FormWithCoefficients.
*/
double SmallSimplexPartition3D::evaluate0Form(const Buffer<double>& discreteForm, const Vector3& evaluationPoint) const {

	if (!mesh_old_ptr)
		return 0.0; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find the element containing evaluationPoint
	static uint element = 0;
	if (!mesh_old.findElement(evaluationPoint4, element))
		return 0.0; //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getBodyNodesCustom(element);
	Buffer<double> nodeCoefficients(nodes.size()); //coefficients for small nodes that are also big nodes are obtained at once
	uint i, j;
	for (i = 0; i < nodeCoefficients.size(); ++i) {
		nodeCoefficients[i] = discreteForm[nodes[i]];
	}

	Buffer<uint> edges(6);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[0], nodes[3]);
	edges[3] = mesh_old.findEdge(nodes[1], nodes[2]);
	edges[4] = mesh_old.findEdge(nodes[1], nodes[3]);
	edges[5] = mesh_old.findEdge(nodes[2], nodes[3]);
	Buffer<VectorN> edgeCoefficients(edges.size()); //coefficients for small nodes that are on big edges
	for (i = 0; i < edgeCoefficients.size(); ++i) {
		Buffer<uint> edgeNodes = mesh_old.getEdgeNodes(edges[i]);
		VectorN edgeValues = edgeValuesNode0 * discreteForm[edgeNodes[0]] + edgeValuesNode1 * discreteForm[edgeNodes[1]];
		for (j = 0; j < activeSmallNodesInEdges.size(); ++j) {
			edgeValues[j] += discreteForm[smallNodesEdgeList(edges[i], j)];
		}
		solveLUP(matrix0FormsEdges, matrix0FormsEdges_p, edgeValues, edgeCoefficients[i]);
	}

	Buffer<uint> faces(4);
	faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
	faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
	faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
	faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));
	Buffer<VectorN> faceCoefficients(faces.size()); //coefficients for small nodes that are on big faces
	for (i = 0; i < faceCoefficients.size(); ++i) {
		const Buffer<uint> faceNodes = mesh_old.getFaceNodes(faces[i]);
		Buffer<uint> faceEdges(3);
		faceEdges[0] = mesh_old.findEdge(faceNodes[0], faceNodes[1]);
		faceEdges[1] = mesh_old.findEdge(faceNodes[0], faceNodes[2]);
		faceEdges[2] = mesh_old.findEdge(faceNodes[1], faceNodes[2]);
		VectorN faceValues = faceValuesNode0 * discreteForm[faceNodes[0]] + faceValuesNode1 * discreteForm[faceNodes[1]] + faceValuesNode2 * discreteForm[faceNodes[2]];
		uint permutationIndexEdge0 = getPermutation(mesh_old.getEdgeNodes(faceEdges[0]), { faceNodes[0], faceNodes[1] });
		uint permutationIndexEdge1 = getPermutation(mesh_old.getEdgeNodes(faceEdges[1]), { faceNodes[0], faceNodes[2] });
		uint permutationIndexEdge2 = getPermutation(mesh_old.getEdgeNodes(faceEdges[2]), { faceNodes[1], faceNodes[2] });
		for (j = 0; j < activeSmallNodesInFaces.size(); ++j) {
			faceValues[j] += faceNodeValuesEdge0[permutationIndexEdge0][j].dot(edgeCoefficients[edges.findFirst(faceEdges[0])]);
			faceValues[j] += faceNodeValuesEdge1[permutationIndexEdge1][j].dot(edgeCoefficients[edges.findFirst(faceEdges[1])]);
			faceValues[j] += faceNodeValuesEdge2[permutationIndexEdge2][j].dot(edgeCoefficients[edges.findFirst(faceEdges[2])]);
			faceValues[j] += discreteForm[smallNodesFaceList(faces[i], j)];
		}
		solveLUP(matrix0FormsFaces, matrix0FormsFaces_p, faceValues, faceCoefficients[i]);
	}

	VectorN bodyValues = bodyValuesNode0 * discreteForm[nodes[0]] + bodyValuesNode1 * discreteForm[nodes[1]]
		+ bodyValuesNode2 * discreteForm[nodes[2]] + bodyValuesNode3 * discreteForm[nodes[3]];
	uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
	uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
	uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[0], nodes[3] });
	uint permutationEdge3 = getPermutation(mesh_old.getEdgeNodes(edges[3]), { nodes[1], nodes[2] });
	uint permutationEdge4 = getPermutation(mesh_old.getEdgeNodes(edges[4]), { nodes[1], nodes[3] });
	uint permutationEdge5 = getPermutation(mesh_old.getEdgeNodes(edges[5]), { nodes[2], nodes[3] });
	uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
	uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
	uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
	uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });
	for (i = 0; i < activeSmallNodesInBodies.size(); ++i) {
		bodyValues[i] += bodyNodeValuesEdge0[permutationEdge0][i].dot(edgeCoefficients[0]);
		bodyValues[i] += bodyNodeValuesEdge1[permutationEdge1][i].dot(edgeCoefficients[1]);
		bodyValues[i] += bodyNodeValuesEdge2[permutationEdge2][i].dot(edgeCoefficients[2]);
		bodyValues[i] += bodyNodeValuesEdge3[permutationEdge3][i].dot(edgeCoefficients[3]);
		bodyValues[i] += bodyNodeValuesEdge4[permutationEdge4][i].dot(edgeCoefficients[4]);
		bodyValues[i] += bodyNodeValuesEdge5[permutationEdge5][i].dot(edgeCoefficients[5]);
		bodyValues[i] += bodyNodeValuesFace0[permutationFace0][i].dot(faceCoefficients[0]);
		bodyValues[i] += bodyNodeValuesFace1[permutationFace1][i].dot(faceCoefficients[1]);
		bodyValues[i] += bodyNodeValuesFace2[permutationFace2][i].dot(faceCoefficients[2]);
		bodyValues[i] += bodyNodeValuesFace3[permutationFace3][i].dot(faceCoefficients[3]);
		bodyValues[i] += discreteForm[smallNodesBodyList(element, i)];
	}
	VectorN bodyCoefficients; //coefficients for small nodes in the interior of the tetrahedron
	solveLUP(matrix0FormsBodies, matrix0FormsBodies_p, bodyValues, bodyCoefficients);

	//find barycentric coordinates of evaluationPoint
	double h0 = mesh_old.getFaceDeviation(faces[3], mesh_old.getNodePosition(nodes[0])).len();
	double lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() / h0;
	double h1 = mesh_old.getFaceDeviation(faces[2], mesh_old.getNodePosition(nodes[1])).len();
	double lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() / h1;
	double h2 = mesh_old.getFaceDeviation(faces[1], mesh_old.getNodePosition(nodes[2])).len();
	double lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() / h2;
	double h3 = mesh_old.getFaceDeviation(faces[0], mesh_old.getNodePosition(nodes[3])).len();
	double lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() / h3;

	//compute the result using the coefficients solved above
	std::vector<double> barycentricCoords;

	//contribution from small nodes that are also big nodes
	double result = std::pow(lambda0, order) * nodeCoefficients[0] + std::pow(lambda1, order) * nodeCoefficients[1]
		+ std::pow(lambda2, order) * nodeCoefficients[2] + std::pow(lambda3, order) * nodeCoefficients[3];

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
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda3 }, edgePermutations[permutationEdge2]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[2][i];
	}
	//edge3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[permutationEdge3]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[3][i];
	}
	//edge4
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda3 }, edgePermutations[permutationEdge4]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[4][i];
	}
	//edge5
	barycentricCoords = permute(std::vector<double>{ lambda2, lambda3 }, edgePermutations[permutationEdge5]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[5][i];
	}

	//small nodes on faces
	//face0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[permutationFace0]);
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[0][i];
	}
	//face1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[permutationFace1]);
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[1][i];
	}
	//face2
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[permutationFace2]);
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[2][i];
	}
	//face3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[permutationFace3]);
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[3][i];
	}

	//small nodes in the interior
	for (i = 0; i < activeSmallNodesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInBodies[i].multiIndex;
		result += std::pow(lambda0, mi[0] + 1) * std::pow(lambda1, mi[1]) * std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[i];
	}

	return result;
}

/*
Computes the value of the (higher order) Whitney 1-form (interpolant of discreteForm) at evaluationPoint.
When evaluating many times, it is more efficient to first solve the coefficients with solve1FormCoefficients and then evaluate with evaluate1FormWithCoefficients.
*/
Vector3 SmallSimplexPartition3D::evaluate1Form(const Buffer<double>& discreteForm, const Vector3& evaluationPoint) const {

	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find the element containing evaluationPoint
	static uint element = 0;
	if (!mesh_old.findElement(evaluationPoint4, element))
		return Vector3(0.0, 0.0, 0.0); //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getBodyNodesCustom(element);
	uint i, j;

	Buffer<uint> edges(6);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[0], nodes[3]);
	edges[3] = mesh_old.findEdge(nodes[1], nodes[2]);
	edges[4] = mesh_old.findEdge(nodes[1], nodes[3]);
	edges[5] = mesh_old.findEdge(nodes[2], nodes[3]);
	Buffer<VectorN> edgeCoefficients(edges.size()); //coefficients for small edges that are on big edges
	for (i = 0; i < edgeCoefficients.size(); ++i) {
		VectorN edgeValues(activeSmallEdgesInEdges.size(), 0.0);
		for (j = 0; j < activeSmallEdgesInEdges.size(); ++j) {
			edgeValues[j] = discreteForm[smallEdgesEdgeList(edges[i], j)];
		}
		solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeValues, edgeCoefficients[i]);
	}

	Buffer<uint> faces(4);
	faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
	faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
	faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
	faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));
	Buffer<VectorN> faceCoefficients(faces.size()); //coefficients for small edges that are on big faces
	for (i = 0; i < faceCoefficients.size(); ++i) {
		const Buffer<uint> faceNodes = mesh_old.getFaceNodes(faces[i]);
		Buffer<uint> faceEdges(3);
		faceEdges[0] = mesh_old.findEdge(faceNodes[0], faceNodes[1]);
		faceEdges[1] = mesh_old.findEdge(faceNodes[0], faceNodes[2]);
		faceEdges[2] = mesh_old.findEdge(faceNodes[1], faceNodes[2]);
		VectorN faceValues(activeSmallEdgesInFaces.size(), 0.0);
		uint permutationIndexEdge0 = getPermutation(mesh_old.getEdgeNodes(faceEdges[0]), { faceNodes[0], faceNodes[1] });
		uint permutationIndexEdge1 = getPermutation(mesh_old.getEdgeNodes(faceEdges[1]), { faceNodes[0], faceNodes[2] });
		uint permutationIndexEdge2 = getPermutation(mesh_old.getEdgeNodes(faceEdges[2]), { faceNodes[1], faceNodes[2] });
		for (j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
			faceValues[j] += faceEdgeValuesEdge0[permutationIndexEdge0][j].dot(edgeCoefficients[edges.findFirst(faceEdges[0])]);
			faceValues[j] += faceEdgeValuesEdge1[permutationIndexEdge1][j].dot(edgeCoefficients[edges.findFirst(faceEdges[1])]);
			faceValues[j] += faceEdgeValuesEdge2[permutationIndexEdge2][j].dot(edgeCoefficients[edges.findFirst(faceEdges[2])]);
			faceValues[j] += discreteForm[smallEdgesFaceList(faces[i], j)];
		}
		solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues, faceCoefficients[i]);
	}

	VectorN bodyValues(activeSmallEdgesInBodies.size(), 0.0);
	uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
	uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
	uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[0], nodes[3] });
	uint permutationEdge3 = getPermutation(mesh_old.getEdgeNodes(edges[3]), { nodes[1], nodes[2] });
	uint permutationEdge4 = getPermutation(mesh_old.getEdgeNodes(edges[4]), { nodes[1], nodes[3] });
	uint permutationEdge5 = getPermutation(mesh_old.getEdgeNodes(edges[5]), { nodes[2], nodes[3] });
	uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
	uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
	uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
	uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });
	for (i = 0; i < activeSmallEdgesInBodies.size(); ++i) {
		bodyValues[i] += bodyEdgeValuesEdge0[permutationEdge0][i].dot(edgeCoefficients[0]);
		bodyValues[i] += bodyEdgeValuesEdge1[permutationEdge1][i].dot(edgeCoefficients[1]);
		bodyValues[i] += bodyEdgeValuesEdge2[permutationEdge2][i].dot(edgeCoefficients[2]);
		bodyValues[i] += bodyEdgeValuesEdge3[permutationEdge3][i].dot(edgeCoefficients[3]);
		bodyValues[i] += bodyEdgeValuesEdge4[permutationEdge4][i].dot(edgeCoefficients[4]);
		bodyValues[i] += bodyEdgeValuesEdge5[permutationEdge5][i].dot(edgeCoefficients[5]);
		bodyValues[i] += bodyEdgeValuesFace0[permutationFace0][i].dot(faceCoefficients[0]);
		bodyValues[i] += bodyEdgeValuesFace1[permutationFace1][i].dot(faceCoefficients[1]);
		bodyValues[i] += bodyEdgeValuesFace2[permutationFace2][i].dot(faceCoefficients[2]);
		bodyValues[i] += bodyEdgeValuesFace3[permutationFace3][i].dot(faceCoefficients[3]);
		bodyValues[i] += discreteForm[smallEdgesBodyList(element, i)];
	}
	VectorN bodyCoefficients;
	solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);

	//find barycentric coordinates and their gradients at evaluationPoint
	double lambda0, lambda1, lambda2, lambda3;
	Vector3 gradLambda0 = mesh_old.getFaceDeviation(faces[3], mesh_old.getNodePosition(nodes[0])).toVector3();
	double lensq0 = gradLambda0.lensq();
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() / std::sqrt(lensq0);
	gradLambda0 /= lensq0;
	Vector3 gradLambda1 = mesh_old.getFaceDeviation(faces[2], mesh_old.getNodePosition(nodes[1])).toVector3();
	double lensq1 = gradLambda1.lensq();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() / std::sqrt(lensq1);
	gradLambda1 /= lensq1;
	Vector3 gradLambda2 = mesh_old.getFaceDeviation(faces[1], mesh_old.getNodePosition(nodes[2])).toVector3();
	double lensq2 = gradLambda2.lensq();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() / std::sqrt(lensq2);
	gradLambda2 /= lensq2;
	Vector3 gradLambda3 = mesh_old.getFaceDeviation(faces[0], mesh_old.getNodePosition(nodes[3])).toVector3();
	double lensq3 = gradLambda3.lensq();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() / std::sqrt(lensq3);
	gradLambda3 /= lensq3;

	//compute the values of the lowest order Whitney forms
	Vector3 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector3 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector3 w2 = lambda0 * gradLambda3 - lambda3 * gradLambda0;
	Vector3 w3 = lambda1 * gradLambda2 - lambda2 * gradLambda1;
	Vector3 w4 = lambda1 * gradLambda3 - lambda3 * gradLambda1;
	Vector3 w5 = lambda2 * gradLambda3 - lambda3 * gradLambda2;

	//compute the result using the coefficients solved above
	Buffer<double> d(6, 0.0);
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
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda3 }, edgePermutations[permutationEdge2]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[2][i];
	}
	//edge3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[permutationEdge3]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[3] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[3][i];
	}
	//edge4
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda3 }, edgePermutations[permutationEdge4]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[4] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[4][i];
	}
	//edge5
	barycentricCoords = permute(std::vector<double>{ lambda2, lambda3 }, edgePermutations[permutationEdge5]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[5] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[5][i];
	}
	//take orientations into account
	if (getEdgePermutationSign(permutationEdge0) < 0)
		d[0] = -d[0];
	if (getEdgePermutationSign(permutationEdge1) < 0)
		d[1] = -d[1];
	if (getEdgePermutationSign(permutationEdge2) < 0)
		d[2] = -d[2];
	if (getEdgePermutationSign(permutationEdge3) < 0)
		d[3] = -d[3];
	if (getEdgePermutationSign(permutationEdge4) < 0)
		d[4] = -d[4];
	if (getEdgePermutationSign(permutationEdge5) < 0)
		d[5] = -d[5];

	//small simplices of faces
	Buffer<double> faceContributions(3, 0.0); //contributions to the edges of the faces
	//face 0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[permutationFace0]);
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
			* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[0][i];
	}
	//distribute the contributions for the correct edges
	switch (permutationFace0) {
	case 0:
		d[0] += faceContributions[0];
		d[1] += faceContributions[1];
		d[3] += faceContributions[2];
		break;
	case 1:
		d[0] += faceContributions[1];
		d[1] += faceContributions[0];
		d[3] += -faceContributions[2];
		break;
	case 2:
		d[0] += -faceContributions[0];
		d[1] += faceContributions[2];
		d[3] += faceContributions[1];
		break;
	case 3:
		d[0] += -faceContributions[1];
		d[1] += -faceContributions[2];
		d[3] += faceContributions[0];
		break;
	case 4:
		d[0] += faceContributions[2];
		d[1] += -faceContributions[0];
		d[3] += -faceContributions[1];
		break;
	case 5:
		d[0] += -faceContributions[2];
		d[1] += -faceContributions[1];
		d[3] += -faceContributions[0];
		break;
	}
	faceContributions.fill(0.0);
	//face1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[permutationFace1]);
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
			* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[1][i];
	}
	switch (permutationFace1) {
	case 0:
		d[0] += faceContributions[0];
		d[2] += faceContributions[1];
		d[4] += faceContributions[2];
		break;
	case 1:
		d[0] += faceContributions[1];
		d[2] += faceContributions[0];
		d[4] += -faceContributions[2];
		break;
	case 2:
		d[0] += -faceContributions[0];
		d[2] += faceContributions[2];
		d[4] += faceContributions[1];
		break;
	case 3:
		d[0] += -faceContributions[1];
		d[2] += -faceContributions[2];
		d[4] += faceContributions[0];
		break;
	case 4:
		d[0] += faceContributions[2];
		d[2] += -faceContributions[0];
		d[4] += -faceContributions[1];
		break;
	case 5:
		d[0] += -faceContributions[2];
		d[2] += -faceContributions[1];
		d[4] += -faceContributions[0];
		break;
	}
	faceContributions.fill(0.0);
	//face2
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[permutationFace2]);
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
			* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[2][i];
	}
	switch (permutationFace2) {
	case 0:
		d[1] += faceContributions[0];
		d[2] += faceContributions[1];
		d[5] += faceContributions[2];
		break;
	case 1:
		d[1] += faceContributions[1];
		d[2] += faceContributions[0];
		d[5] += -faceContributions[2];
		break;
	case 2:
		d[1] += -faceContributions[0];
		d[2] += faceContributions[2];
		d[5] += faceContributions[1];
		break;
	case 3:
		d[1] += -faceContributions[1];
		d[2] += -faceContributions[2];
		d[5] += faceContributions[0];
		break;
	case 4:
		d[1] += faceContributions[2];
		d[2] += -faceContributions[0];
		d[5] += -faceContributions[1];
		break;
	case 5:
		d[1] += -faceContributions[2];
		d[2] += -faceContributions[1];
		d[5] += -faceContributions[0];
		break;
	}
	faceContributions.fill(0.0);
	//face3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[permutationFace3]);
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
			* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[3][i];
	}
	switch (permutationFace3) {
	case 0:
		d[3] += faceContributions[0];
		d[4] += faceContributions[1];
		d[5] += faceContributions[2];
		break;
	case 1:
		d[3] += faceContributions[1];
		d[4] += faceContributions[0];
		d[5] += -faceContributions[2];
		break;
	case 2:
		d[3] += -faceContributions[0];
		d[4] += faceContributions[2];
		d[5] += faceContributions[1];
		break;
	case 3:
		d[3] += -faceContributions[1];
		d[4] += -faceContributions[2];
		d[5] += faceContributions[0];
		break;
	case 4:
		d[3] += faceContributions[2];
		d[4] += -faceContributions[0];
		d[5] += -faceContributions[1];
		break;
	case 5:
		d[3] += -faceContributions[2];
		d[4] += -faceContributions[1];
		d[5] += -faceContributions[0];
		break;
	}

	//small simplices in the interior
	for (i = 0; i < activeSmallEdgesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInBodies[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInBodies[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3 + d[4] * w4 + d[5] * w5;
}

/*
Computes the value of the (higher order) Whitney 2-form (interpolant of discreteForm) at evaluationPoint.
When evaluating many times, it is more efficient to first solve the coefficients with solve2FormCoefficients and then evaluate with evaluate2FormWithCoefficients.
*/
Vector3 SmallSimplexPartition3D::evaluate2Form(const Buffer<double>& discreteForm, const Vector3& evaluationPoint) const {

	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find the element containing evaluationPoint
	static uint element = 0;
	if (!mesh_old.findElement(evaluationPoint4, element))
		return Vector3(0.0, 0.0, 0.0); //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getBodyNodesCustom(element);
	uint i, j;

	Buffer<uint> faces(4);
	faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
	faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
	faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
	faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));
	Buffer<VectorN> faceCoefficients(faces.size()); //coefficients for small faces that are on big faces
	for (i = 0; i < faceCoefficients.size(); ++i) {
		VectorN faceValues(activeSmallFacesInFaces.size(), 0.0);
		for (j = 0; j < activeSmallFacesInFaces.size(); ++j) {
			faceValues[j] = discreteForm[smallFacesFaceList(faces[i], j)];
		}
		solveLUP(matrix2FormsFaces, matrix2FormsFaces_p, faceValues, faceCoefficients[i]);
	}

	VectorN bodyValues(activeSmallFacesInBodies.size(), 0.0);
	uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
	uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
	uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
	uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });
	for (i = 0; i < activeSmallFacesInBodies.size(); ++i) {
		bodyValues[i] += bodyFaceValuesFace0[permutationFace0][i].dot(faceCoefficients[0]);
		bodyValues[i] += bodyFaceValuesFace1[permutationFace1][i].dot(faceCoefficients[1]);
		bodyValues[i] += bodyFaceValuesFace2[permutationFace2][i].dot(faceCoefficients[2]);
		bodyValues[i] += bodyFaceValuesFace3[permutationFace3][i].dot(faceCoefficients[3]);
		bodyValues[i] += discreteForm[smallFacesBodyList(element, i)];
	}
	VectorN bodyCoefficients; //coefficients for small faces in the interior of the tetrahedron
	solveLUP(matrix2FormsBodies, matrix2FormsBodies_p, bodyValues, bodyCoefficients);

	//find barycentric coordinates and their gradients at evaluationPoint
	double lambda0, lambda1, lambda2, lambda3;
	Vector3 gradLambda0 = mesh_old.getFaceDeviation(faces[3], mesh_old.getNodePosition(nodes[0])).toVector3();
	double lensq0 = gradLambda0.lensq();
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() / std::sqrt(lensq0);
	gradLambda0 /= lensq0;
	Vector3 gradLambda1 = mesh_old.getFaceDeviation(faces[2], mesh_old.getNodePosition(nodes[1])).toVector3();
	double lensq1 = gradLambda1.lensq();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() / std::sqrt(lensq1);
	gradLambda1 /= lensq1;
	Vector3 gradLambda2 = mesh_old.getFaceDeviation(faces[1], mesh_old.getNodePosition(nodes[2])).toVector3();
	double lensq2 = gradLambda2.lensq();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() / std::sqrt(lensq2);
	gradLambda2 /= lensq2;
	Vector3 gradLambda3 = mesh_old.getFaceDeviation(faces[0], mesh_old.getNodePosition(nodes[3])).toVector3();
	double lensq3 = gradLambda3.lensq();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() / std::sqrt(lensq3);
	gradLambda3 /= lensq3;

	//compute the values of the lowest order Whitney forms
	double multiplier = 2.0 / ThreeVector3(mesh_old.getNodePosition3(nodes[1]) - mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[2]) - mesh_old.getNodePosition3(nodes[0]),
		mesh_old.getNodePosition3(nodes[3]) - mesh_old.getNodePosition3(nodes[0])).xyz;
	Vector3 w0 = multiplier * (mesh_old.getNodePosition3(nodes[3]) - evaluationPoint);
	Vector3 w1 = -multiplier * (mesh_old.getNodePosition3(nodes[2]) - evaluationPoint);
	Vector3 w2 = multiplier * (mesh_old.getNodePosition3(nodes[1]) - evaluationPoint);
	Vector3 w3 = -multiplier * (mesh_old.getNodePosition3(nodes[0]) - evaluationPoint);

	//compute the result using the coefficients solved above
	Buffer<double> d(4, 0.0);
	std::vector<double> barycentricCoords;

	//small simplices of faces
	//face0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[permutationFace0]);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		d[0] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[0][i];
	}
	//face1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[permutationFace1]);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		d[1] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[1][i];
	}
	//face2
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[permutationFace2]);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[2][i];
	}
	//face3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[permutationFace3]);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		d[3] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[3][i];
	}
	//take orientations into account
	if (getFacePermutationSign(permutationFace0) < 0)
		d[0] = -d[0];
	if (getFacePermutationSign(permutationFace1) < 0)
		d[1] = -d[1];
	if (getFacePermutationSign(permutationFace2) < 0)
		d[2] = -d[2];
	if (getFacePermutationSign(permutationFace3) < 0)
		d[3] = -d[3];

	//small simplices in the interior
	for (i = 0; i < activeSmallFacesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInBodies[i].multiIndex;
		d[getFaceIndex(*activeSmallFacesInBodies[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3;
}

/*
Computes the value of the (higher order) Whitney 3-form (interpolant of discreteForm) at evaluationPoint.
When evaluating many times, it is more efficient to first solve the coefficients with solve3FormCoefficients and then evaluate with evaluate3FormWithCoefficients.
N.B. discreteForm[i] is interpreted with the positive sign as if getBodyVector(i).xyz > 0.
*/
double SmallSimplexPartition3D::evaluate3Form(const Buffer<double>& discreteForm, const Vector3& evaluationPoint) const {

	if (!mesh_old_ptr)
		return 0.0; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find the element containing evaluationPoint
	static uint element = 0;
	if (!mesh_old.findElement(evaluationPoint4, element))
		return 0.0; //evaluationPoint is not contained in mesh_old

	uint i;
	//solve the coefficients
	VectorN bodyValues(activeSmallBodiesInBodies.size(), 0.0);
	for (i = 0; i < activeSmallBodiesInBodies.size(); ++i) {
		bodyValues[i] = discreteForm[smallBodiesBodyList(element, i)];
	}
	VectorN bodyCoefficients;
	solveLUP(matrix3Forms, matrix3Forms_p, bodyValues, bodyCoefficients);

	Buffer<uint> nodes = getBodyNodesCustom(element);
	//find barycentric coordinates of evaluationPoint
	double lambda0, lambda1, lambda2, lambda3;
	findBarycentricCoordinates(mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[1]), mesh_old.getNodePosition3(nodes[2]), mesh_old.getNodePosition3(nodes[3]),
		evaluationPoint4.toVector3(), lambda0, lambda1, lambda2, lambda3);

	double result = 0.0; //compute the result using the coefficients solved above

	for (i = 0; i < activeSmallBodiesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallBodiesInBodies[i].multiIndex;
		result += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1]) * std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[i];
	}
	return result / std::abs(mesh_old.getBodyVector3(element).determinant());
}

//solves the coefficients of the (higher order) Whitney 0-form (interpolant of discreteForm)
void SmallSimplexPartition3D::solve0FormCoefficients(const Buffer<double>& discreteForm, Buffer<double>& nodeCoefficients, Buffer<VectorN>& edgeCoefficients,
	Buffer<VectorN>& faceCoefficients, Buffer<VectorN>& bodyCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallSimplexPartition3D has not been associated with a mesh

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
		const Buffer<uint> faceNodes = mesh_old.getFaceNodes(i);
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

	bodyCoefficients.resize(mesh_old.getBodySize()); //coefficients for small nodes in the interior of the tetrahedron
	for (i = 0; i < bodyCoefficients.size(); ++i) {
		const Buffer<uint> nodes = getBodyNodesCustom(i);
		Buffer<uint> edges(6);
		edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
		edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
		edges[2] = mesh_old.findEdge(nodes[0], nodes[3]);
		edges[3] = mesh_old.findEdge(nodes[1], nodes[2]);
		edges[4] = mesh_old.findEdge(nodes[1], nodes[3]);
		edges[5] = mesh_old.findEdge(nodes[2], nodes[3]);
		Buffer<uint> faces(4);
		faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
		faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
		faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
		faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));
		VectorN bodyValues = bodyValuesNode0 * discreteForm[nodes[0]] + bodyValuesNode1 * discreteForm[nodes[1]]
			+ bodyValuesNode2 * discreteForm[nodes[2]] + bodyValuesNode3 * discreteForm[nodes[3]];
		uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
		uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
		uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[0], nodes[3] });
		uint permutationEdge3 = getPermutation(mesh_old.getEdgeNodes(edges[3]), { nodes[1], nodes[2] });
		uint permutationEdge4 = getPermutation(mesh_old.getEdgeNodes(edges[4]), { nodes[1], nodes[3] });
		uint permutationEdge5 = getPermutation(mesh_old.getEdgeNodes(edges[5]), { nodes[2], nodes[3] });
		uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
		uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
		uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
		uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });
		for (j = 0; j < activeSmallNodesInBodies.size(); ++j) {
			bodyValues[j] += bodyNodeValuesEdge0[permutationEdge0][j].dot(edgeCoefficients[edges[0]]);
			bodyValues[j] += bodyNodeValuesEdge1[permutationEdge1][j].dot(edgeCoefficients[edges[1]]);
			bodyValues[j] += bodyNodeValuesEdge2[permutationEdge2][j].dot(edgeCoefficients[edges[2]]);
			bodyValues[j] += bodyNodeValuesEdge3[permutationEdge3][j].dot(edgeCoefficients[edges[3]]);
			bodyValues[j] += bodyNodeValuesEdge4[permutationEdge4][j].dot(edgeCoefficients[edges[4]]);
			bodyValues[j] += bodyNodeValuesEdge5[permutationEdge5][j].dot(edgeCoefficients[edges[5]]);
			bodyValues[j] += bodyNodeValuesFace0[permutationFace0][j].dot(faceCoefficients[faces[0]]);
			bodyValues[j] += bodyNodeValuesFace1[permutationFace1][j].dot(faceCoefficients[faces[1]]);
			bodyValues[j] += bodyNodeValuesFace2[permutationFace2][j].dot(faceCoefficients[faces[2]]);
			bodyValues[j] += bodyNodeValuesFace3[permutationFace3][j].dot(faceCoefficients[faces[3]]);
			bodyValues[j] += discreteForm[smallNodesBodyList(i, j)];
		}
		solveLUP(matrix0FormsBodies, matrix0FormsBodies_p, bodyValues, bodyCoefficients[i]);
	}
}

//solves the coefficients of the (higher order) Whitney 1-form (interpolant of discreteForm)
void SmallSimplexPartition3D::solve1FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients,
	Buffer<VectorN>& bodyCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallSimplexPartition3D has not been associated with a mesh

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
		const Buffer<uint> faceNodes = mesh_old.getFaceNodes(i);
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

	bodyCoefficients.resize(mesh_old.getBodySize()); //coefficients for small edges in the interior of the tetrahedron
	for (i = 0; i < bodyCoefficients.size(); ++i) {
		const Buffer<uint> nodes = getBodyNodesCustom(i);
		Buffer<uint> edges(6);
		edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
		edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
		edges[2] = mesh_old.findEdge(nodes[0], nodes[3]);
		edges[3] = mesh_old.findEdge(nodes[1], nodes[2]);
		edges[4] = mesh_old.findEdge(nodes[1], nodes[3]);
		edges[5] = mesh_old.findEdge(nodes[2], nodes[3]);
		Buffer<uint> faces(4);
		faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
		faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
		faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
		faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));
		VectorN bodyValues(activeSmallEdgesInBodies.size(), 0.0);
		uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
		uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
		uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[0], nodes[3] });
		uint permutationEdge3 = getPermutation(mesh_old.getEdgeNodes(edges[3]), { nodes[1], nodes[2] });
		uint permutationEdge4 = getPermutation(mesh_old.getEdgeNodes(edges[4]), { nodes[1], nodes[3] });
		uint permutationEdge5 = getPermutation(mesh_old.getEdgeNodes(edges[5]), { nodes[2], nodes[3] });
		uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
		uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
		uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
		uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });
		for (j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
			bodyValues[j] += bodyEdgeValuesEdge0[permutationEdge0][j].dot(edgeCoefficients[edges[0]]);
			bodyValues[j] += bodyEdgeValuesEdge1[permutationEdge1][j].dot(edgeCoefficients[edges[1]]);
			bodyValues[j] += bodyEdgeValuesEdge2[permutationEdge2][j].dot(edgeCoefficients[edges[2]]);
			bodyValues[j] += bodyEdgeValuesEdge3[permutationEdge3][j].dot(edgeCoefficients[edges[3]]);
			bodyValues[j] += bodyEdgeValuesEdge4[permutationEdge4][j].dot(edgeCoefficients[edges[4]]);
			bodyValues[j] += bodyEdgeValuesEdge5[permutationEdge5][j].dot(edgeCoefficients[edges[5]]);
			bodyValues[j] += bodyEdgeValuesFace0[permutationFace0][j].dot(faceCoefficients[faces[0]]);
			bodyValues[j] += bodyEdgeValuesFace1[permutationFace1][j].dot(faceCoefficients[faces[1]]);
			bodyValues[j] += bodyEdgeValuesFace2[permutationFace2][j].dot(faceCoefficients[faces[2]]);
			bodyValues[j] += bodyEdgeValuesFace3[permutationFace3][j].dot(faceCoefficients[faces[3]]);
			bodyValues[j] += discreteForm[smallEdgesBodyList(i, j)];
		}
		solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients[i]);
	}
}

//solves the coefficients of the (higher order) Whitney 2-form (interpolant of discreteForm)
void SmallSimplexPartition3D::solve2FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& faceCoefficients, Buffer<VectorN>& bodyCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallSimplexPartition3D has not been associated with a mesh

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

	bodyCoefficients.resize(mesh_old.getBodySize()); //coefficients for small faces in the interior of the tetrahedron
	for (i = 0; i < bodyCoefficients.size(); ++i) {
		const Buffer<uint> nodes = getBodyNodesCustom(i);
		Buffer<uint> faces(4);
		faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
		faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
		faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
		faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));
		VectorN bodyValues(activeSmallFacesInBodies.size(), 0.0);
		uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
		uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
		uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
		uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });
		for (j = 0; j < activeSmallFacesInBodies.size(); ++j) {
			bodyValues[j] += bodyFaceValuesFace0[permutationFace0][j].dot(faceCoefficients[faces[0]]);
			bodyValues[j] += bodyFaceValuesFace1[permutationFace1][j].dot(faceCoefficients[faces[1]]);
			bodyValues[j] += bodyFaceValuesFace2[permutationFace2][j].dot(faceCoefficients[faces[2]]);
			bodyValues[j] += bodyFaceValuesFace3[permutationFace3][j].dot(faceCoefficients[faces[3]]);
			bodyValues[j] += discreteForm[smallFacesBodyList(i, j)];
		}
		solveLUP(matrix2FormsBodies, matrix2FormsBodies_p, bodyValues, bodyCoefficients[i]);
	}
}

//solves the coefficients of the (higher order) Whitney 3-form (interpolant of discreteForm)
void SmallSimplexPartition3D::solve3FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& bodyCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	uint i, j;

	bodyCoefficients.resize(mesh_old.getBodySize()); //solve the coefficients
	for (i = 0; i < bodyCoefficients.size(); ++i) {
		VectorN bodyValues(activeSmallBodiesInBodies.size(), 0.0);
		for (j = 0; j < activeSmallBodiesInBodies.size(); ++j) {
			bodyValues[j] = discreteForm[smallBodiesBodyList(i, j)];
		}
		solveLUP(matrix3Forms, matrix3Forms_p, bodyValues, bodyCoefficients[i]);
	}
}

//computes the value of the (higher order) Whitney 0-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known
double SmallSimplexPartition3D::evaluate0FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<double>& nodeCoefficients, const Buffer<VectorN>& edgeCoefficients,
	const Buffer<VectorN>& faceCoefficients, const Buffer<VectorN>& bodyCoefficients, uint element) const {

	if (!mesh_old_ptr)
		return 0.0; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4
	uint i;

	//find the element containing evaluationPoint
	if (!mesh_old.findElement(evaluationPoint4, element))
		return 0.0; //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getBodyNodesCustom(element);
	Buffer<uint> edges(6);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[0], nodes[3]);
	edges[3] = mesh_old.findEdge(nodes[1], nodes[2]);
	edges[4] = mesh_old.findEdge(nodes[1], nodes[3]);
	edges[5] = mesh_old.findEdge(nodes[2], nodes[3]);
	Buffer<uint> faces(4);
	faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
	faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
	faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
	faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));

	uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
	uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
	uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[0], nodes[3] });
	uint permutationEdge3 = getPermutation(mesh_old.getEdgeNodes(edges[3]), { nodes[1], nodes[2] });
	uint permutationEdge4 = getPermutation(mesh_old.getEdgeNodes(edges[4]), { nodes[1], nodes[3] });
	uint permutationEdge5 = getPermutation(mesh_old.getEdgeNodes(edges[5]), { nodes[2], nodes[3] });
	uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
	uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
	uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
	uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });

	//find barycentric coordinates of evaluationPoint
	double h0 = mesh_old.getFaceDeviation(faces[3], mesh_old.getNodePosition(nodes[0])).len();
	double lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() / h0;
	double h1 = mesh_old.getFaceDeviation(faces[2], mesh_old.getNodePosition(nodes[1])).len();
	double lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() / h1;
	double h2 = mesh_old.getFaceDeviation(faces[1], mesh_old.getNodePosition(nodes[2])).len();
	double lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() / h2;
	double h3 = mesh_old.getFaceDeviation(faces[0], mesh_old.getNodePosition(nodes[3])).len();
	double lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() / h3;

	//compute the result using the given coefficients
	std::vector<double> barycentricCoords;

	//contribution from small nodes that are also big nodes
	double result = std::pow(lambda0, order) * nodeCoefficients[nodes[0]] + std::pow(lambda1, order) * nodeCoefficients[nodes[1]]
		+ std::pow(lambda2, order) * nodeCoefficients[nodes[2]] + std::pow(lambda3, order) * nodeCoefficients[nodes[3]];

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
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda3 }, edgePermutations[permutationEdge2]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[2]][i];
	}
	//edge3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[permutationEdge3]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[3]][i];
	}
	//edge4
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda3 }, edgePermutations[permutationEdge4]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[4]][i];
	}
	//edge5
	barycentricCoords = permute(std::vector<double>{ lambda2, lambda3 }, edgePermutations[permutationEdge5]);
	for (i = 0; i < activeSmallNodesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInEdges[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[5]][i];
	}

	//small nodes on faces
	//face0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[permutationFace0]);
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[0]][i];
	}
	//face1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[permutationFace1]);
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[1]][i];
	}
	//face2
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[permutationFace2]);
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[2]][i];
	}
	//face3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[permutationFace3]);
	for (i = 0; i < activeSmallNodesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInFaces[i].multiIndex;
		result += std::pow(barycentricCoords[0], mi[0] + 1) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[3]][i];
	}

	//small nodes in the interior
	for (i = 0; i < activeSmallNodesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallNodesInBodies[i].multiIndex;
		result += std::pow(lambda0, mi[0] + 1) * std::pow(lambda1, mi[1]) * std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[element][i];
	}

	return result;
}

//computes the value of the (higher order) Whitney 1-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known
Vector3 SmallSimplexPartition3D::evaluate1FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<VectorN>& edgeCoefficients,
	const Buffer<VectorN>& faceCoefficients, const Buffer<VectorN>& bodyCoefficients, uint element) const {

	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4
	uint i;

	//find the element containing evaluationPoint
	if (!mesh_old.findElement(evaluationPoint4, element))
		return Vector3(0.0, 0.0, 0.0); //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getBodyNodesCustom(element);
	Buffer<uint> edges(6);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[0], nodes[3]);
	edges[3] = mesh_old.findEdge(nodes[1], nodes[2]);
	edges[4] = mesh_old.findEdge(nodes[1], nodes[3]);
	edges[5] = mesh_old.findEdge(nodes[2], nodes[3]);
	Buffer<uint> faces(4);
	faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
	faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
	faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
	faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));

	uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
	uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
	uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[0], nodes[3] });
	uint permutationEdge3 = getPermutation(mesh_old.getEdgeNodes(edges[3]), { nodes[1], nodes[2] });
	uint permutationEdge4 = getPermutation(mesh_old.getEdgeNodes(edges[4]), { nodes[1], nodes[3] });
	uint permutationEdge5 = getPermutation(mesh_old.getEdgeNodes(edges[5]), { nodes[2], nodes[3] });
	uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
	uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
	uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
	uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });

	//find barycentric coordinates and their gradients at evaluationPoint
	double lambda0, lambda1, lambda2, lambda3;
	Vector3 gradLambda0 = mesh_old.getFaceDeviation(faces[3], mesh_old.getNodePosition(nodes[0])).toVector3();
	double lensq0 = gradLambda0.lensq();
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() / std::sqrt(lensq0);
	gradLambda0 /= lensq0;
	Vector3 gradLambda1 = mesh_old.getFaceDeviation(faces[2], mesh_old.getNodePosition(nodes[1])).toVector3();
	double lensq1 = gradLambda1.lensq();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() / std::sqrt(lensq1);
	gradLambda1 /= lensq1;
	Vector3 gradLambda2 = mesh_old.getFaceDeviation(faces[1], mesh_old.getNodePosition(nodes[2])).toVector3();
	double lensq2 = gradLambda2.lensq();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() / std::sqrt(lensq2);
	gradLambda2 /= lensq2;
	Vector3 gradLambda3 = mesh_old.getFaceDeviation(faces[0], mesh_old.getNodePosition(nodes[3])).toVector3();
	double lensq3 = gradLambda3.lensq();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() / std::sqrt(lensq3);
	gradLambda3 /= lensq3;

	//compute the values of the lowest order Whitney forms
	Vector3 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector3 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector3 w2 = lambda0 * gradLambda3 - lambda3 * gradLambda0;
	Vector3 w3 = lambda1 * gradLambda2 - lambda2 * gradLambda1;
	Vector3 w4 = lambda1 * gradLambda3 - lambda3 * gradLambda1;
	Vector3 w5 = lambda2 * gradLambda3 - lambda3 * gradLambda2;

	//compute the result using the coefficients solved above
	Buffer<double> d(6, 0.0);
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
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda3 }, edgePermutations[permutationEdge2]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[2]][i];
	}
	//edge3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[permutationEdge3]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[3] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[3]][i];
	}
	//edge4
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda3 }, edgePermutations[permutationEdge4]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[4] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[4]][i];
	}
	//edge5
	barycentricCoords = permute(std::vector<double>{ lambda2, lambda3 }, edgePermutations[permutationEdge5]);
	for (i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
		d[5] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[edges[5]][i];
	}
	//take orientations into account
	if (getEdgePermutationSign(permutationEdge0) < 0)
		d[0] = -d[0];
	if (getEdgePermutationSign(permutationEdge1) < 0)
		d[1] = -d[1];
	if (getEdgePermutationSign(permutationEdge2) < 0)
		d[2] = -d[2];
	if (getEdgePermutationSign(permutationEdge3) < 0)
		d[3] = -d[3];
	if (getEdgePermutationSign(permutationEdge4) < 0)
		d[4] = -d[4];
	if (getEdgePermutationSign(permutationEdge5) < 0)
		d[5] = -d[5];

	//small simplices of faces
	Buffer<double> faceContributions(3, 0.0); //contributions to the edges of the faces
	//face 0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[permutationFace0]);
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
			* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[0]][i];
	}
	//distribute the contributions for the correct edges
	switch (permutationFace0) {
	case 0:
		d[0] += faceContributions[0];
		d[1] += faceContributions[1];
		d[3] += faceContributions[2];
		break;
	case 1:
		d[0] += faceContributions[1];
		d[1] += faceContributions[0];
		d[3] += -faceContributions[2];
		break;
	case 2:
		d[0] += -faceContributions[0];
		d[1] += faceContributions[2];
		d[3] += faceContributions[1];
		break;
	case 3:
		d[0] += -faceContributions[1];
		d[1] += -faceContributions[2];
		d[3] += faceContributions[0];
		break;
	case 4:
		d[0] += faceContributions[2];
		d[1] += -faceContributions[0];
		d[3] += -faceContributions[1];
		break;
	case 5:
		d[0] += -faceContributions[2];
		d[1] += -faceContributions[1];
		d[3] += -faceContributions[0];
		break;
	}
	faceContributions.fill(0.0);
	//face1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[permutationFace1]);
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
			* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[1]][i];
	}
	switch (permutationFace1) {
	case 0:
		d[0] += faceContributions[0];
		d[2] += faceContributions[1];
		d[4] += faceContributions[2];
		break;
	case 1:
		d[0] += faceContributions[1];
		d[2] += faceContributions[0];
		d[4] += -faceContributions[2];
		break;
	case 2:
		d[0] += -faceContributions[0];
		d[2] += faceContributions[2];
		d[4] += faceContributions[1];
		break;
	case 3:
		d[0] += -faceContributions[1];
		d[2] += -faceContributions[2];
		d[4] += faceContributions[0];
		break;
	case 4:
		d[0] += faceContributions[2];
		d[2] += -faceContributions[0];
		d[4] += -faceContributions[1];
		break;
	case 5:
		d[0] += -faceContributions[2];
		d[2] += -faceContributions[1];
		d[4] += -faceContributions[0];
		break;
	}
	faceContributions.fill(0.0);
	//face2
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[permutationFace2]);
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
			* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[2]][i];
	}
	switch (permutationFace2) {
	case 0:
		d[1] += faceContributions[0];
		d[2] += faceContributions[1];
		d[5] += faceContributions[2];
		break;
	case 1:
		d[1] += faceContributions[1];
		d[2] += faceContributions[0];
		d[5] += -faceContributions[2];
		break;
	case 2:
		d[1] += -faceContributions[0];
		d[2] += faceContributions[2];
		d[5] += faceContributions[1];
		break;
	case 3:
		d[1] += -faceContributions[1];
		d[2] += -faceContributions[2];
		d[5] += faceContributions[0];
		break;
	case 4:
		d[1] += faceContributions[2];
		d[2] += -faceContributions[0];
		d[5] += -faceContributions[1];
		break;
	case 5:
		d[1] += -faceContributions[2];
		d[2] += -faceContributions[1];
		d[5] += -faceContributions[0];
		break;
	}
	faceContributions.fill(0.0);
	//face3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[permutationFace3]);
	for (i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
		faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
			* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[3]][i];
	}
	switch (permutationFace3) {
	case 0:
		d[3] += faceContributions[0];
		d[4] += faceContributions[1];
		d[5] += faceContributions[2];
		break;
	case 1:
		d[3] += faceContributions[1];
		d[4] += faceContributions[0];
		d[5] += -faceContributions[2];
		break;
	case 2:
		d[3] += -faceContributions[0];
		d[4] += faceContributions[2];
		d[5] += faceContributions[1];
		break;
	case 3:
		d[3] += -faceContributions[1];
		d[4] += -faceContributions[2];
		d[5] += faceContributions[0];
		break;
	case 4:
		d[3] += faceContributions[2];
		d[4] += -faceContributions[0];
		d[5] += -faceContributions[1];
		break;
	case 5:
		d[3] += -faceContributions[2];
		d[4] += -faceContributions[1];
		d[5] += -faceContributions[0];
		break;
	}

	//small simplices in the interior
	for (i = 0; i < activeSmallEdgesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInBodies[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInBodies[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[element][i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3 + d[4] * w4 + d[5] * w5;
}

//computes the value of the (higher order) Whitney 2-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known
Vector3 SmallSimplexPartition3D::evaluate2FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<VectorN>& faceCoefficients, const Buffer<VectorN>& bodyCoefficients, uint element) const {

	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4
	uint i;

	//find the element containing evaluationPoint
	if (!mesh_old.findElement(evaluationPoint4, element))
		return Vector3(0.0, 0.0, 0.0); //evaluationPoint is not contained in mesh_old

	const Buffer<uint> nodes = getBodyNodesCustom(element);
	Buffer<uint> faces(4);
	faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
	faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
	faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
	faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));

	uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
	uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
	uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
	uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });

	//find barycentric coordinates and their gradients at evaluationPoint
	double lambda0, lambda1, lambda2, lambda3;
	Vector3 gradLambda0 = mesh_old.getFaceDeviation(faces[3], mesh_old.getNodePosition(nodes[0])).toVector3();
	double lensq0 = gradLambda0.lensq();
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() / std::sqrt(lensq0);
	gradLambda0 /= lensq0;
	Vector3 gradLambda1 = mesh_old.getFaceDeviation(faces[2], mesh_old.getNodePosition(nodes[1])).toVector3();
	double lensq1 = gradLambda1.lensq();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() / std::sqrt(lensq1);
	gradLambda1 /= lensq1;
	Vector3 gradLambda2 = mesh_old.getFaceDeviation(faces[1], mesh_old.getNodePosition(nodes[2])).toVector3();
	double lensq2 = gradLambda2.lensq();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() / std::sqrt(lensq2);
	gradLambda2 /= lensq2;
	Vector3 gradLambda3 = mesh_old.getFaceDeviation(faces[0], mesh_old.getNodePosition(nodes[3])).toVector3();
	double lensq3 = gradLambda3.lensq();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() / std::sqrt(lensq3);
	gradLambda3 /= lensq3;

	//compute the values of the lowest order Whitney forms
	double multiplier = 2.0 / ThreeVector3(mesh_old.getNodePosition3(nodes[1]) - mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[2]) - mesh_old.getNodePosition3(nodes[0]),
		mesh_old.getNodePosition3(nodes[3]) - mesh_old.getNodePosition3(nodes[0])).xyz;
	Vector3 w0 = multiplier * (mesh_old.getNodePosition3(nodes[3]) - evaluationPoint);
	Vector3 w1 = -multiplier * (mesh_old.getNodePosition3(nodes[2]) - evaluationPoint);
	Vector3 w2 = multiplier * (mesh_old.getNodePosition3(nodes[1]) - evaluationPoint);
	Vector3 w3 = -multiplier * (mesh_old.getNodePosition3(nodes[0]) - evaluationPoint);

	//compute the result using the coefficients solved above
	Buffer<double> d(4, 0.0);
	std::vector<double> barycentricCoords;

	//small simplices of faces
	//face0
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[permutationFace0]);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		d[0] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[0]][i];
	}
	//face1
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[permutationFace1]);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		d[1] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[1]][i];
	}
	//face2
	barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[permutationFace2]);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[2]][i];
	}
	//face3
	barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[permutationFace3]);
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInFaces[i].multiIndex;
		d[3] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[faces[3]][i];
	}
	//take orientations into account
	if (getFacePermutationSign(permutationFace0) < 0)
		d[0] = -d[0];
	if (getFacePermutationSign(permutationFace1) < 0)
		d[1] = -d[1];
	if (getFacePermutationSign(permutationFace2) < 0)
		d[2] = -d[2];
	if (getFacePermutationSign(permutationFace3) < 0)
		d[3] = -d[3];

	//small simplices in the interior
	for (i = 0; i < activeSmallFacesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallFacesInBodies[i].multiIndex;
		d[getFaceIndex(*activeSmallFacesInBodies[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[element][i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3;
}

/*
Computes the value of the (higher order) Whitney 3-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known.
N.B. discreteForm[i] is interpreted with the positive sign as if getBodyVector(i).xyz > 0.
*/
double SmallSimplexPartition3D::evaluate3FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<VectorN>& bodyCoefficients, uint element) const {

	if (!mesh_old_ptr)
		return 0.0; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4
	uint i;

	//find the element containing evaluationPoint
	if (!mesh_old.findElement(evaluationPoint4, element))
		return 0.0; //evaluationPoint is not contained in mesh_old

	Buffer<uint> nodes = getBodyNodesCustom(element);
	//find barycentric coordinates of evaluationPoint
	double lambda0, lambda1, lambda2, lambda3;
	findBarycentricCoordinates(mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[1]), mesh_old.getNodePosition3(nodes[2]), mesh_old.getNodePosition3(nodes[3]),
		evaluationPoint4.toVector3(), lambda0, lambda1, lambda2, lambda3);

	double result = 0.0; //compute the result using the given coefficients

	for (i = 0; i < activeSmallBodiesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallBodiesInBodies[i].multiIndex;
		result += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1]) * std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[element][i];
	}
	return result / std::abs(mesh_old.getBodyVector3(element).determinant());
}

//Gives a strict upper bound for the number of nonzero elements in the higher order Hodge matrix for 1-forms. The actual number is smaller only if some integrals evaluate exactly to zero.
uint SmallSimplexPartition3D::estimateNonzeros() const {
	if (!mesh_old_ptr)
		return 0; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	uint a_e = activeSmallEdgesInEdges.size();
	uint a_f = activeSmallEdgesInFaces.size();
	uint a_b = activeSmallEdgesInBodies.size();
	uint smallEdgesInBodies = 6 * bodyEdgeHoles.size() + bodyFaceHoles.size();
	uint smallEdgesInFaces = 3 * faceEdgeHoles.size();
	uint smallEdgesInEdges = activeSmallEdgesInEdges.size();
	uint estimate = mesh_old.getBodySize() * smallEdgesInBodies * (6 * a_e + 4 * a_f + a_b);
	for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
		if (mesh_old.getFaceBodies(i).size() == 1) { //boundary face
			estimate += smallEdgesInFaces * (6 * a_e + 4 * a_f + a_b);
		}
		else {
			estimate += smallEdgesInFaces * (9 * a_e + 7 * a_f + 2 * a_b);
		}
	}
	for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
		const Buffer<uint> bodies = mesh_old.getEdgeBodies(i);
		Buffer<uint> faces;
		Buffer<uint> edges;
		for (uint b = 0; b < bodies.size(); ++b) {
			faces = mesh_old.getBodyFaces(bodies[b]).getUnion(faces);
			edges = mesh_old.getBodyEdges(bodies[b]).getUnion(edges);
		}
		estimate += smallEdgesInEdges * (edges.size() * a_e + faces.size() * a_f + bodies.size() * a_b);
	}
	return estimate;
}

/*
Computes the discrete Hodge for higher order Whitney 1-forms to a sparse matrix (of type Buffer<std::unordered_map<uint, double>>).
Integrals over dual faces are computed in each element separately. If all elements have same shape, use formHodgeMatrix1FormsCustom for efficiency.
*/
void SmallSimplexPartition3D::formHodgeMatrix1Forms(Buffer<std::unordered_map<uint, double>>& star, bool circumcentric, bool minkowskiMetric) const {
	if (!mesh_old_ptr)
		return; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	const BuilderMesh& mesh = *mesh_ptr;
	star.resize(mesh.getEdgeSize());

	struct DualFacePart
	{
		uint index;
		Vector3 p0;
		Vector3 p1;
		Vector3 p2;
		Vector3 normalVector;
	};

	for (uint b = 0; b < mesh_old.getBodySize(); ++b) {
		//find the subsimplices and their orientations
		const Buffer<uint> nodes = getBodyNodesCustom(b);
		Buffer<uint> edges(6);
		edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
		edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
		edges[2] = mesh_old.findEdge(nodes[0], nodes[3]);
		edges[3] = mesh_old.findEdge(nodes[1], nodes[2]);
		edges[4] = mesh_old.findEdge(nodes[1], nodes[3]);
		edges[5] = mesh_old.findEdge(nodes[2], nodes[3]);
		Buffer<uint> faces(4);
		faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
		faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
		faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
		faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));
		Buffer<uint> faceNodes0 = mesh_old.getFaceNodes(faces[0]);
		Buffer<uint> faceNodes1 = mesh_old.getFaceNodes(faces[1]);
		Buffer<uint> faceNodes2 = mesh_old.getFaceNodes(faces[2]);
		Buffer<uint> faceNodes3 = mesh_old.getFaceNodes(faces[3]);
		Buffer<uint> faceEdges0(3);
		faceEdges0[0] = mesh_old.findEdge(faceNodes0[0], faceNodes0[1]);
		faceEdges0[1] = mesh_old.findEdge(faceNodes0[0], faceNodes0[2]);
		faceEdges0[2] = mesh_old.findEdge(faceNodes0[1], faceNodes0[2]);
		uint permutationIndexEdge0Face0 = getPermutation(mesh_old.getEdgeNodes(faceEdges0[0]), { faceNodes0[0], faceNodes0[1] });
		uint permutationIndexEdge1Face0 = getPermutation(mesh_old.getEdgeNodes(faceEdges0[1]), { faceNodes0[0], faceNodes0[2] });
		uint permutationIndexEdge2Face0 = getPermutation(mesh_old.getEdgeNodes(faceEdges0[2]), { faceNodes0[1], faceNodes0[2] });
		Buffer<uint> faceEdges1(3);
		faceEdges1[0] = mesh_old.findEdge(faceNodes1[0], faceNodes1[1]);
		faceEdges1[1] = mesh_old.findEdge(faceNodes1[0], faceNodes1[2]);
		faceEdges1[2] = mesh_old.findEdge(faceNodes1[1], faceNodes1[2]);
		uint permutationIndexEdge0Face1 = getPermutation(mesh_old.getEdgeNodes(faceEdges1[0]), { faceNodes1[0], faceNodes1[1] });
		uint permutationIndexEdge1Face1 = getPermutation(mesh_old.getEdgeNodes(faceEdges1[1]), { faceNodes1[0], faceNodes1[2] });
		uint permutationIndexEdge2Face1 = getPermutation(mesh_old.getEdgeNodes(faceEdges1[2]), { faceNodes1[1], faceNodes1[2] });
		Buffer<uint> faceEdges2(3);
		faceEdges2[0] = mesh_old.findEdge(faceNodes2[0], faceNodes2[1]);
		faceEdges2[1] = mesh_old.findEdge(faceNodes2[0], faceNodes2[2]);
		faceEdges2[2] = mesh_old.findEdge(faceNodes2[1], faceNodes2[2]);
		uint permutationIndexEdge0Face2 = getPermutation(mesh_old.getEdgeNodes(faceEdges2[0]), { faceNodes2[0], faceNodes2[1] });
		uint permutationIndexEdge1Face2 = getPermutation(mesh_old.getEdgeNodes(faceEdges2[1]), { faceNodes2[0], faceNodes2[2] });
		uint permutationIndexEdge2Face2 = getPermutation(mesh_old.getEdgeNodes(faceEdges2[2]), { faceNodes2[1], faceNodes2[2] });
		Buffer<uint> faceEdges3(3);
		faceEdges3[0] = mesh_old.findEdge(faceNodes3[0], faceNodes3[1]);
		faceEdges3[1] = mesh_old.findEdge(faceNodes3[0], faceNodes3[2]);
		faceEdges3[2] = mesh_old.findEdge(faceNodes3[1], faceNodes3[2]);
		uint permutationIndexEdge0Face3 = getPermutation(mesh_old.getEdgeNodes(faceEdges3[0]), { faceNodes3[0], faceNodes3[1] });
		uint permutationIndexEdge1Face3 = getPermutation(mesh_old.getEdgeNodes(faceEdges3[1]), { faceNodes3[0], faceNodes3[2] });
		uint permutationIndexEdge2Face3 = getPermutation(mesh_old.getEdgeNodes(faceEdges3[2]), { faceNodes3[1], faceNodes3[2] });
		uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
		uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
		uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[0], nodes[3] });
		uint permutationEdge3 = getPermutation(mesh_old.getEdgeNodes(edges[3]), { nodes[1], nodes[2] });
		uint permutationEdge4 = getPermutation(mesh_old.getEdgeNodes(edges[4]), { nodes[1], nodes[3] });
		uint permutationEdge5 = getPermutation(mesh_old.getEdgeNodes(edges[5]), { nodes[2], nodes[3] });
		uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
		uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
		uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
		uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });

		//save node positions, and outward normal vectors of the opposite faces, and gradients of barycentric functions
		Vector3 p0 = mesh_old.getNodePosition3(nodes[0]);
		Vector3 p1 = mesh_old.getNodePosition3(nodes[1]);
		Vector3 p2 = mesh_old.getNodePosition3(nodes[2]);
		Vector3 p3 = mesh_old.getNodePosition3(nodes[3]);
		Vector3 p0OppositeNormal = TwoVector3(p2 - p1, p3 - p1).dual();
		if (p0OppositeNormal.dot(p1 - p0) < 0)
			p0OppositeNormal *= -1.0;
		Vector3 p1OppositeNormal = TwoVector3(p2 - p0, p3 - p0).dual();
		if (p1OppositeNormal.dot(p0 - p1) < 0)
			p1OppositeNormal *= -1.0;
		Vector3 p2OppositeNormal = TwoVector3(p1 - p0, p3 - p0).dual();
		if (p2OppositeNormal.dot(p0 - p2) < 0)
			p2OppositeNormal *= -1.0;
		Vector3 p3OppositeNormal = TwoVector3(p1 - p0, p2 - p0).dual();
		if (p3OppositeNormal.dot(p0 - p3) < 0)
			p3OppositeNormal *= -1.0;
		Vector3 gradLambda0 = mesh_old.getFaceDeviation(faces[3], mesh_old.getNodePosition(nodes[0])).toVector3();
		gradLambda0 /= gradLambda0.lensq();
		Vector3 gradLambda1 = mesh_old.getFaceDeviation(faces[2], mesh_old.getNodePosition(nodes[1])).toVector3();
		gradLambda1 /= gradLambda1.lensq();
		Vector3 gradLambda2 = mesh_old.getFaceDeviation(faces[1], mesh_old.getNodePosition(nodes[2])).toVector3();
		gradLambda2 /= gradLambda2.lensq();
		Vector3 gradLambda3 = mesh_old.getFaceDeviation(faces[0], mesh_old.getNodePosition(nodes[3])).toVector3();
		gradLambda3 /= gradLambda3.lensq();

		//form the list of dual face parts contained in body b
		std::vector<DualFacePart> triangleList;
		triangleList.reserve(60 * bodyMultiIndices.size() + 10 * bodyFaceHoles.size());
		for (uint e = 0; e < edges.size(); ++e) {
			for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
				const Buffer<uint> smallBodies = mesh.getEdgeBodies(smallEdgesEdgeList(edges[e], i));
				const Buffer<uint>& smallNodes = mesh.getEdgeNodes(smallEdgesEdgeList(edges[e], i));
				Vector3 n0 = mesh.getNodePosition3(smallNodes[0]);
				Vector3 n1 = mesh.getNodePosition3(smallNodes[1]);
				uint sb = 0;
				while (!isInsideTetrahedron(mesh.getBodyAverage3(smallBodies[sb]), p0, p1, p2, p3, p0OppositeNormal, p1OppositeNormal, p2OppositeNormal, p3OppositeNormal))
					++sb;
				const Buffer<uint> smallFaces = mesh.getEdgeFaces(smallEdgesEdgeList(edges[e], i)).getIntersection(mesh.getBodyFaces(smallBodies[sb]));
				Vector3 sf0_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[0]) : mesh.getFaceAverage3(smallFaces[0]));
				Vector3 sf1_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[1]) : mesh.getFaceAverage3(smallFaces[1]));
				Vector3 sb_bc = (circumcentric ? mesh.getBodyPosition3(smallBodies[sb]) : mesh.getBodyAverage3(smallBodies[sb]));
				Vector3 se_bc = 0.5 * (n0 + n1);
				Vector3 normalVector0 = 0.5 * TwoVector3(sf0_bc - se_bc, sb_bc - se_bc).dual();
				if (normalVector0.dot(n1 - n0) < 0)
					normalVector0 *= -1.0;
				if (minkowskiMetric)
					normalVector0.z = -normalVector0.z;
				Vector3 normalVector1 = 0.5 * TwoVector3(sf1_bc - se_bc, sb_bc - se_bc).dual();
				if (normalVector1.dot(n1 - n0) < 0)
					normalVector1 *= -1.0;
				if (minkowskiMetric)
					normalVector1.z = -normalVector1.z;
				triangleList.push_back({ smallEdgesEdgeList(edges[e], i), sf0_bc, sb_bc, se_bc, normalVector0 });
				triangleList.push_back({ smallEdgesEdgeList(edges[e], i), sf1_bc, sb_bc, se_bc, normalVector1 });
			}
		}
		for (uint f = 0; f < faces.size(); ++f) {
			for (uint i = firstEdgeOfFace(faces[f]); i < firstEdgeOfFace(faces[f]) + 3 * faceEdgeHoles.size(); ++i) {
				const Buffer<uint> smallBodies = mesh.getEdgeBodies(i);
				const Buffer<uint>& smallNodes = mesh.getEdgeNodes(i);
				Vector3 n0 = mesh.getNodePosition3(smallNodes[0]);
				Vector3 n1 = mesh.getNodePosition3(smallNodes[1]);
				for (uint sb = 0; sb < smallBodies.size(); ++sb) {
					if (!isInsideTetrahedron(mesh.getBodyAverage3(smallBodies[sb]), p0, p1, p2, p3, p0OppositeNormal, p1OppositeNormal, p2OppositeNormal, p3OppositeNormal))
						continue;
					const Buffer<uint> smallFaces = mesh.getEdgeFaces(i).getIntersection(mesh.getBodyFaces(smallBodies[sb]));
					Vector3 sf0_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[0]) : mesh.getFaceAverage3(smallFaces[0]));
					Vector3 sf1_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[1]) : mesh.getFaceAverage3(smallFaces[1]));
					Vector3 sb_bc = (circumcentric ? mesh.getBodyPosition3(smallBodies[sb]) : mesh.getBodyAverage3(smallBodies[sb]));
					Vector3 se_bc = 0.5 * (n0 + n1);
					Vector3 normalVector0 = 0.5 * TwoVector3(sf0_bc - se_bc, sb_bc - se_bc).dual();
					if (normalVector0.dot(n1 - n0) < 0)
						normalVector0 *= -1.0;
					if (minkowskiMetric)
						normalVector0.z = -normalVector0.z;
					Vector3 normalVector1 = 0.5 * TwoVector3(sf1_bc - se_bc, sb_bc - se_bc).dual();
					if (normalVector1.dot(n1 - n0) < 0)
						normalVector1 *= -1.0;
					if (minkowskiMetric)
						normalVector1.z = -normalVector1.z;
					triangleList.push_back({ i, sf0_bc, sb_bc, se_bc, normalVector0 });
					triangleList.push_back({ i, sf1_bc, sb_bc, se_bc, normalVector1 });
				}
			}
		}
		for (uint i = firstEdgeOfBody(b); i < firstEdgeOfBody(b) + 6 * bodyEdgeHoles.size() + bodyFaceHoles.size(); ++i) {
			const Buffer<uint> smallBodies = mesh.getEdgeBodies(i);
			const Buffer<uint>& smallNodes = mesh.getEdgeNodes(i);
			Vector3 n0 = mesh.getNodePosition3(smallNodes[0]);
			Vector3 n1 = mesh.getNodePosition3(smallNodes[1]);
			for (uint sb = 0; sb < smallBodies.size(); ++sb) {
				const Buffer<uint> smallFaces = mesh.getEdgeFaces(i).getIntersection(mesh.getBodyFaces(smallBodies[sb]));
				Vector3 sf0_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[0]) : mesh.getFaceAverage3(smallFaces[0]));
				Vector3 sf1_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[1]) : mesh.getFaceAverage3(smallFaces[1]));
				Vector3 sb_bc = (circumcentric ? mesh.getBodyPosition3(smallBodies[sb]) : mesh.getBodyAverage3(smallBodies[sb]));
				Vector3 se_bc = 0.5 * (n0 + n1);
				Vector3 normalVector0 = 0.5 * TwoVector3(sf0_bc - se_bc, sb_bc - se_bc).dual();
				if (normalVector0.dot(n1 - n0) < 0)
					normalVector0 *= -1.0;
				if (minkowskiMetric)
					normalVector0.z = -normalVector0.z;
				Vector3 normalVector1 = 0.5 * TwoVector3(sf1_bc - se_bc, sb_bc - se_bc).dual();
				if (normalVector1.dot(n1 - n0) < 0)
					normalVector1 *= -1.0;
				if (minkowskiMetric)
					normalVector1.z = -normalVector1.z;
				triangleList.push_back({ i, sf0_bc, sb_bc, se_bc, normalVector0 });
				triangleList.push_back({ i, sf1_bc, sb_bc, se_bc, normalVector1 });
			}
		}

		//integrate all basis functions over the triangles in triangleList and add the values in Hodge matrix
		VectorN edgeCochain(activeSmallEdgesInEdges.size(), 0.0);
		VectorN faceCochain(activeSmallEdgesInFaces.size(), 0.0);
		VectorN bodyCochain(activeSmallEdgesInBodies.size(), 0.0);
		VectorN edgeCoefficients(activeSmallEdgesInEdges.size(), 0.0);
		VectorN faceValues0(activeSmallEdgesInFaces.size(), 0.0);
		VectorN faceValues1(activeSmallEdgesInFaces.size(), 0.0);
		VectorN faceCoefficients0(activeSmallEdgesInFaces.size(), 0.0);
		VectorN faceCoefficients1(activeSmallEdgesInFaces.size(), 0.0);
		VectorN faceCoefficients(activeSmallEdgesInFaces.size(), 0.0);
		VectorN bodyValues(activeSmallEdgesInBodies.size(), 0.0);
		VectorN bodyCoefficients(activeSmallEdgesInBodies.size(), 0.0);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) { //edge basis functions
			edgeCochain[i] = 1.0;
			solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeCochain, edgeCoefficients);
			//edge0 (belongs to face0 and face1)
			if (edges[0] == faceEdges0[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face0][j].dot(edgeCoefficients);
				}
			else if (edges[0] == faceEdges0[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face0][j].dot(edgeCoefficients);
				}
			else if (edges[0] == faceEdges0[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face0][j].dot(edgeCoefficients);
				}
			if (edges[0] == faceEdges1[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face1][j].dot(edgeCoefficients);
				}
			else if (edges[0] == faceEdges1[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face1][j].dot(edgeCoefficients);
				}
			else if (edges[0] == faceEdges1[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face1][j].dot(edgeCoefficients);
				}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues0, faceCoefficients0);
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues1, faceCoefficients1);
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesEdge0[permutationEdge0][j].dot(edgeCoefficients);
				bodyValues[j] += bodyEdgeValuesFace0[permutationFace0][j].dot(faceCoefficients0);
				bodyValues[j] += bodyEdgeValuesFace1[permutationFace1][j].dot(faceCoefficients1);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 0, permutationEdge0, faceCoefficients0, 0, permutationFace0,
								faceCoefficients1, 1, permutationFace1, bodyCoefficients, nodes, edges, faces,
								gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesEdgeList(edges[0], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			faceValues0.val.fill(0.0);
			faceValues1.val.fill(0.0);
			faceCoefficients0.val.fill(0.0);
			faceCoefficients1.val.fill(0.0);
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			//edge1 (belongs to face0 and face2)
			if (edges[1] == faceEdges0[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face0][j].dot(edgeCoefficients);
				}
			else if (edges[1] == faceEdges0[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face0][j].dot(edgeCoefficients);
				}
			else if (edges[1] == faceEdges0[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face0][j].dot(edgeCoefficients);
				}
			if (edges[1] == faceEdges2[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face2][j].dot(edgeCoefficients);
				}
			else if (edges[1] == faceEdges2[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face2][j].dot(edgeCoefficients);
				}
			else if (edges[1] == faceEdges2[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face2][j].dot(edgeCoefficients);
				}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues0, faceCoefficients0);
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues1, faceCoefficients1);
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesEdge1[permutationEdge1][j].dot(edgeCoefficients);
				bodyValues[j] += bodyEdgeValuesFace0[permutationFace0][j].dot(faceCoefficients0);
				bodyValues[j] += bodyEdgeValuesFace2[permutationFace2][j].dot(faceCoefficients1);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 1, permutationEdge1, faceCoefficients0, 0, permutationFace0,
								faceCoefficients1, 2, permutationFace2, bodyCoefficients, nodes, edges, faces,
								gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesEdgeList(edges[1], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			faceValues0.val.fill(0.0);
			faceValues1.val.fill(0.0);
			faceCoefficients0.val.fill(0.0);
			faceCoefficients1.val.fill(0.0);
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			//edge2 (belongs to face1 and face2)
			if (edges[2] == faceEdges1[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face1][j].dot(edgeCoefficients);
				}
			else if (edges[2] == faceEdges1[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face1][j].dot(edgeCoefficients);
				}
			else if (edges[2] == faceEdges1[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face1][j].dot(edgeCoefficients);
				}
			if (edges[2] == faceEdges2[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face2][j].dot(edgeCoefficients);
				}
			else if (edges[2] == faceEdges2[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face2][j].dot(edgeCoefficients);
				}
			else if (edges[2] == faceEdges2[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face2][j].dot(edgeCoefficients);
				}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues0, faceCoefficients0);
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues1, faceCoefficients1);
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesEdge2[permutationEdge2][j].dot(edgeCoefficients);
				bodyValues[j] += bodyEdgeValuesFace1[permutationFace1][j].dot(faceCoefficients0);
				bodyValues[j] += bodyEdgeValuesFace2[permutationFace2][j].dot(faceCoefficients1);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 2, permutationEdge2, faceCoefficients0, 1, permutationFace1,
								faceCoefficients1, 2, permutationFace2, bodyCoefficients, nodes, edges, faces,
								gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesEdgeList(edges[2], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			faceValues0.val.fill(0.0);
			faceValues1.val.fill(0.0);
			faceCoefficients0.val.fill(0.0);
			faceCoefficients1.val.fill(0.0);
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			//edge3 (belongs to face0 and face3)
			if (edges[3] == faceEdges0[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face0][j].dot(edgeCoefficients);
				}
			else if (edges[3] == faceEdges0[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face0][j].dot(edgeCoefficients);
				}
			else if (edges[3] == faceEdges0[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face0][j].dot(edgeCoefficients);
				}
			if (edges[3] == faceEdges3[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face3][j].dot(edgeCoefficients);
				}
			else if (edges[3] == faceEdges3[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face3][j].dot(edgeCoefficients);
				}
			else if (edges[3] == faceEdges3[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face3][j].dot(edgeCoefficients);
				}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues0, faceCoefficients0);
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues1, faceCoefficients1);
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesEdge3[permutationEdge3][j].dot(edgeCoefficients);
				bodyValues[j] += bodyEdgeValuesFace0[permutationFace0][j].dot(faceCoefficients0);
				bodyValues[j] += bodyEdgeValuesFace3[permutationFace3][j].dot(faceCoefficients1);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 3, permutationEdge3, faceCoefficients0, 0, permutationFace0,
								faceCoefficients1, 3, permutationFace3, bodyCoefficients, nodes, edges, faces,
								gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesEdgeList(edges[3], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			faceValues0.val.fill(0.0);
			faceValues1.val.fill(0.0);
			faceCoefficients0.val.fill(0.0);
			faceCoefficients1.val.fill(0.0);
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			//edge4 (belongs to face1 and face3)
			if (edges[4] == faceEdges1[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face1][j].dot(edgeCoefficients);
				}
			else if (edges[4] == faceEdges1[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face1][j].dot(edgeCoefficients);
				}
			else if (edges[4] == faceEdges1[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face1][j].dot(edgeCoefficients);
				}
			if (edges[4] == faceEdges3[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face3][j].dot(edgeCoefficients);
				}
			else if (edges[4] == faceEdges3[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face3][j].dot(edgeCoefficients);
				}
			else if (edges[4] == faceEdges3[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face3][j].dot(edgeCoefficients);
				}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues0, faceCoefficients0);
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues1, faceCoefficients1);
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesEdge4[permutationEdge4][j].dot(edgeCoefficients);
				bodyValues[j] += bodyEdgeValuesFace1[permutationFace1][j].dot(faceCoefficients0);
				bodyValues[j] += bodyEdgeValuesFace3[permutationFace3][j].dot(faceCoefficients1);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 4, permutationEdge4, faceCoefficients0, 1, permutationFace1,
								faceCoefficients1, 3, permutationFace3, bodyCoefficients, nodes, edges, faces,
								gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesEdgeList(edges[4], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			faceValues0.val.fill(0.0);
			faceValues1.val.fill(0.0);
			faceCoefficients0.val.fill(0.0);
			faceCoefficients1.val.fill(0.0);
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			//edge5 (belongs to face2 and face3)
			if (edges[5] == faceEdges2[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face2][j].dot(edgeCoefficients);
				}
			else if (edges[5] == faceEdges2[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face2][j].dot(edgeCoefficients);
				}
			else if (edges[5] == faceEdges2[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face2][j].dot(edgeCoefficients);
				}
			if (edges[5] == faceEdges3[0])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face3][j].dot(edgeCoefficients);
				}
			else if (edges[5] == faceEdges3[1])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face3][j].dot(edgeCoefficients);
				}
			else if (edges[5] == faceEdges3[2])
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face3][j].dot(edgeCoefficients);
				}
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues0, faceCoefficients0);
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceValues1, faceCoefficients1);
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesEdge5[permutationEdge5][j].dot(edgeCoefficients);
				bodyValues[j] += bodyEdgeValuesFace2[permutationFace2][j].dot(faceCoefficients0);
				bodyValues[j] += bodyEdgeValuesFace3[permutationFace3][j].dot(faceCoefficients1);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, edgeCoefficients, 5, permutationEdge5, faceCoefficients0, 2, permutationFace2,
								faceCoefficients1, 3, permutationFace3, bodyCoefficients, nodes, edges, faces,
								gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesEdgeList(edges[5], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			edgeCoefficients.val.fill(0.0);
			faceValues0.val.fill(0.0);
			faceValues1.val.fill(0.0);
			faceCoefficients0.val.fill(0.0);
			faceCoefficients1.val.fill(0.0);
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			edgeCochain[i] = 0.0;
		}
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) { //face basis functions
			faceCochain[i] = 1.0;
			solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceCochain, faceCoefficients);
			//face0
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesFace0[permutationFace0][j].dot(faceCoefficients);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, faceCoefficients, 0, permutationFace0,
								bodyCoefficients, nodes, edges, faces, gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesFaceList(faces[0], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			//face1
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesFace1[permutationFace1][j].dot(faceCoefficients);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, faceCoefficients, 1, permutationFace1,
								bodyCoefficients, nodes, edges, faces, gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesFaceList(faces[1], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			//face2
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesFace2[permutationFace2][j].dot(faceCoefficients);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, faceCoefficients, 2, permutationFace2,
								bodyCoefficients, nodes, edges, faces, gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesFaceList(faces[2], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			//face3
			for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
				bodyValues[j] += bodyEdgeValuesFace3[permutationFace3][j].dot(faceCoefficients);
			}
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithCoefficients(p, faceCoefficients, 3, permutationFace3,
								bodyCoefficients, nodes, edges, faces, gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesFaceList(faces[3], i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			faceCoefficients.val.fill(0.0);
			bodyValues.val.fill(0.0);
			bodyCoefficients.val.fill(0.0);
			faceCochain[i] = 0.0;
		}
		for (uint i = 0; i < activeSmallEdgesInBodies.size(); ++i) { //body basis functions
			bodyCochain[i] = 1.0;
			solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyCochain, bodyCoefficients);
			for (uint j = 0; j < triangleList.size(); ++j) {
				std::function<double(Vector3)> interpolant{
						[&](Vector3 p) -> double {
							return evaluate1FormWithBodyCoefficients(p, bodyCoefficients, nodes, edges, faces,
								gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(triangleList[j].normalVector);
						}
				};
				star[triangleList[j].index][smallEdgesBodyList(b, i)] += gfd::integralAverage(interpolant, triangleList[j].p0, triangleList[j].p1, triangleList[j].p2);
			}
			bodyCoefficients.val.fill(0.0);
			bodyCochain[i] = 0.0;
		}
	}
}

/*
Reimplementation of Mesh::getBodyNodes that defines a definite order of vertices for the specific element type currently in use.
When all elements are congruent, this enables one to compute the integrals for Hodge matrix only once in a single element.
*/
Buffer<uint> SmallSimplexPartition3D::getBodyNodesCustom(uint b) const {
	//this should be properly reimplemented, but the following implemention currently works for "BccTetrahedron" and "CubeTetrahedron"
	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Buffer<uint> result(4);
	const Buffer<uint>& faces = mesh_old.getBodyFaces(b);
	const Buffer<uint>& edges0 = mesh_old.getFaceEdges(faces[0]);
	double maxLensq = 0.0;
	uint maxInd = 0;
	for (uint i = 0; i < edges0.size(); ++i) {
		const Buffer<uint>& nodes = mesh_old.getEdgeNodes(edges0[i]);
		double lensq = (mesh_old.getNodePosition3(nodes[1]) - mesh_old.getNodePosition3(nodes[0])).lensq();
		if (lensq > maxLensq) {
			maxLensq = lensq;
			maxInd = i;
		}
	}
	const Buffer<uint>& longestEdgeNodes = mesh_old.getEdgeNodes(edges0[maxInd]);
	result[0] = longestEdgeNodes[0];
	result[3] = longestEdgeNodes[1];
	const Buffer<uint>& nextEdgeNodes = mesh_old.getEdgeNodes(edges0[(maxInd + 1) % 3]);
	if (nextEdgeNodes[0] == result[0] || nextEdgeNodes[0] == result[3])
		result[1] = nextEdgeNodes[1];
	else
		result[1] = nextEdgeNodes[0];
	const Buffer<uint>& edges1 = mesh_old.getFaceEdges(faces[1]);
	for (uint i = 0; i < edges1.size(); ++i) {
		const Buffer<uint>& nodes = mesh_old.getEdgeNodes(edges1[i]);
		if (!(nodes[0] == result[0] || nodes[0] == result[1] || nodes[0] == result[3])) {
			result[2] = nodes[0];
			break;
		}
		else if (!(nodes[1] == result[0] || nodes[1] == result[1] || nodes[1] == result[3])) {
			result[2] = nodes[1];
			break;
		}
	}
	const Vector3 p0 = mesh_old.getNodePosition3(result[0]);
	const Vector3 p1 = mesh_old.getNodePosition3(result[1]);
	const Vector3 p2 = mesh_old.getNodePosition3(result[2]);
	const Vector3 p3 = mesh_old.getNodePosition3(result[3]);
	if ((p3 - p0).dot(TwoVector3(p1 - p0, p2 - p0).dual()) < 0) {
		uint swap = result[2];
		result[2] = result[1];
		result[1] = swap;
	}
	return result;
}

/*
Computes the discrete Hodge (with respect to Euclidean metric) for higher order Whitney 1-forms to a sparse matrix (of type Buffer<std::unordered_map<uint, double>>).
Integrals over dual faces are computed only once in a single element.
Requires that all elements have same shape and getBodyNodesCustom has been implemented for that element shape.
Optionally the required integrals can be saved to or loaded from files in folder given by foldername.
*/
void SmallSimplexPartition3D::formHodgeMatrix1FormsCustom(Buffer<std::unordered_map<uint, double>>& star, std::string foldername, bool circumcentric, bool loadIntegralsFromFile) const {
	if (!mesh_old_ptr)
		return; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Buffer<uint> refElementIndices(1, 0);
	Buffer<Buffer<uint>> refElements(1);
	refElements[0].resize(mesh_old.getBodySize());
	for (uint i = 0; i < mesh_old.getBodySize(); ++i) {
		refElements[0][i] = i;
	}
	formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, foldername, circumcentric, false, loadIntegralsFromFile);
}

/*
Computes the discrete Hodge for higher order Whitney 1-forms to a sparse matrix (of type Buffer<std::unordered_map<uint, double>>).
Integrals over dual faces are computed once in each reference element whose index is given in refElementIndices.
Requires that every element should be in one refElements buffer; the computations done in the corresponding reference element are used for all elements in the buffer.
Optionally the required integrals can be saved to or loaded from files in folder given by foldername.
*/
void SmallSimplexPartition3D::formHodgeMatrix1FormsCustom(Buffer<std::unordered_map<uint, double>>& star, const Buffer<uint>& refElementIndices, const Buffer<Buffer<uint>>& refElements,
	std::string foldername, bool circumcentric, bool minkowskiMetric, bool loadIntegralsFromFile) const {
	if (!mesh_old_ptr)
		return; //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	const BuilderMesh& mesh = *mesh_ptr;
	star.resize(mesh.getEdgeSize());
	const uint dualFaceCount = 6 * bodyMultiIndices.size() + bodyFaceHoles.size();
	double scalingFactor = 1.0;

	for (uint ref = 0; ref < refElementIndices.size(); ++ref) {
		Buffer<VectorN> hodgeIntegralsBodyEdges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsFace0Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsFace1Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsFace2Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsFace3Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge0Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge1Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge2Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge3Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge4Edges(dualFaceCount);
		Buffer<Buffer<VectorN>> hodgeIntegralsEdge5Edges(dualFaceCount);
		Buffer<VectorN> edgeCoefficients(activeSmallEdgesInEdges.size());

		if (loadIntegralsFromFile) { //load the integrals from file
			scalingFactor = loadHodgeIntegrals1Forms(foldername, circumcentric, hodgeIntegralsBodyEdges, hodgeIntegralsFace0Edges, hodgeIntegralsFace1Edges, hodgeIntegralsFace2Edges,
				hodgeIntegralsFace3Edges, hodgeIntegralsEdge0Edges, hodgeIntegralsEdge1Edges, hodgeIntegralsEdge2Edges, hodgeIntegralsEdge3Edges, hodgeIntegralsEdge4Edges,
				hodgeIntegralsEdge5Edges, ref);
			VectorN edgeCochain(activeSmallEdgesInEdges.size(), 0.0);
			for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
				edgeCochain[i] = 1.0;
				solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeCochain, edgeCoefficients[i]);
				edgeCochain[i] = 0.0;
			}
		}
		else { //compute the integrals in the reference element
			uint b = refElementIndices[ref];
			const Buffer<uint> body_b_Nodes = getBodyNodesCustom(b);
			Buffer<uint> body_b_Edges(6);
			body_b_Edges[0] = mesh_old.findEdge(body_b_Nodes[0], body_b_Nodes[1]);
			body_b_Edges[1] = mesh_old.findEdge(body_b_Nodes[0], body_b_Nodes[2]);
			body_b_Edges[2] = mesh_old.findEdge(body_b_Nodes[0], body_b_Nodes[3]);
			body_b_Edges[3] = mesh_old.findEdge(body_b_Nodes[1], body_b_Nodes[2]);
			body_b_Edges[4] = mesh_old.findEdge(body_b_Nodes[1], body_b_Nodes[3]);
			body_b_Edges[5] = mesh_old.findEdge(body_b_Nodes[2], body_b_Nodes[3]);
			Buffer<uint> body_b_Faces(4);
			body_b_Faces[0] = mesh_old.findFace(findEdges(body_b_Nodes[0], body_b_Nodes[1], body_b_Nodes[2], mesh_old));
			body_b_Faces[1] = mesh_old.findFace(findEdges(body_b_Nodes[0], body_b_Nodes[1], body_b_Nodes[3], mesh_old));
			body_b_Faces[2] = mesh_old.findFace(findEdges(body_b_Nodes[0], body_b_Nodes[2], body_b_Nodes[3], mesh_old));
			body_b_Faces[3] = mesh_old.findFace(findEdges(body_b_Nodes[1], body_b_Nodes[2], body_b_Nodes[3], mesh_old));
			scalingFactor = (mesh_old.getNodePosition3(body_b_Nodes[1]) - mesh_old.getNodePosition3(body_b_Nodes[0])).len();

			//save node positions, and outward normal vectors of the opposite faces, and gradients of barycentric functions
			Vector3 p0 = mesh_old.getNodePosition3(body_b_Nodes[0]);
			Vector3 p1 = mesh_old.getNodePosition3(body_b_Nodes[1]);
			Vector3 p2 = mesh_old.getNodePosition3(body_b_Nodes[2]);
			Vector3 p3 = mesh_old.getNodePosition3(body_b_Nodes[3]);
			Vector3 p0OppositeNormal = TwoVector3(p2 - p1, p3 - p1).dual();
			if (p0OppositeNormal.dot(p1 - p0) < 0)
				p0OppositeNormal *= -1.0;
			Vector3 p1OppositeNormal = TwoVector3(p2 - p0, p3 - p0).dual();
			if (p1OppositeNormal.dot(p0 - p1) < 0)
				p1OppositeNormal *= -1.0;
			Vector3 p2OppositeNormal = TwoVector3(p1 - p0, p3 - p0).dual();
			if (p2OppositeNormal.dot(p0 - p2) < 0)
				p2OppositeNormal *= -1.0;
			Vector3 p3OppositeNormal = TwoVector3(p1 - p0, p2 - p0).dual();
			if (p3OppositeNormal.dot(p0 - p3) < 0)
				p3OppositeNormal *= -1.0;
			Vector3 gradLambda0 = mesh_old.getFaceDeviation(body_b_Faces[3], mesh_old.getNodePosition(body_b_Nodes[0])).toVector3();
			gradLambda0 /= gradLambda0.lensq();
			Vector3 gradLambda1 = mesh_old.getFaceDeviation(body_b_Faces[2], mesh_old.getNodePosition(body_b_Nodes[1])).toVector3();
			gradLambda1 /= gradLambda1.lensq();
			Vector3 gradLambda2 = mesh_old.getFaceDeviation(body_b_Faces[1], mesh_old.getNodePosition(body_b_Nodes[2])).toVector3();
			gradLambda2 /= gradLambda2.lensq();
			Vector3 gradLambda3 = mesh_old.getFaceDeviation(body_b_Faces[0], mesh_old.getNodePosition(body_b_Nodes[3])).toVector3();
			gradLambda3 /= gradLambda3.lensq();

			//form the list of dual faces
			struct DualFacePart
			{
				Vector3 p0;
				Vector3 p1;
				Vector3 p2;
				Vector3 normalVector;
			};
			std::vector<std::vector<DualFacePart>> dualFaceList;
			dualFaceList.resize(dualFaceCount);
			uint dfIndex = 0;
			for (uint e = 0; e < 6; ++e) { //dual faces of small edges on edges
				for (uint i = 0; i < edgeMultiIndices.size(); ++i) {
					mi_t& mi = edgeMultiIndices[i];
					uint n0Index, n1Index;
					if (e == 0) {
						n0Index = getSmallNodeIndex({ mi[0], mi[1], 0, 0 }, 0, body_b_Nodes);
						n1Index = getSmallNodeIndex({ mi[0], mi[1], 0, 0 }, 1, body_b_Nodes);
					}
					else if (e == 1) {
						n0Index = getSmallNodeIndex({ mi[0], 0, mi[1], 0 }, 0, body_b_Nodes);
						n1Index = getSmallNodeIndex({ mi[0], 0, mi[1], 0 }, 2, body_b_Nodes);
					}
					else if (e == 2) {
						n0Index = getSmallNodeIndex({ mi[0], 0, 0, mi[1] }, 0, body_b_Nodes);
						n1Index = getSmallNodeIndex({ mi[0], 0, 0, mi[1] }, 3, body_b_Nodes);
					}
					else if (e == 3) {
						n0Index = getSmallNodeIndex({ 0, mi[0], mi[1], 0 }, 1, body_b_Nodes);
						n1Index = getSmallNodeIndex({ 0, mi[0], mi[1], 0 }, 2, body_b_Nodes);
					}
					else if (e == 4) {
						n0Index = getSmallNodeIndex({ 0, mi[0], 0, mi[1] }, 1, body_b_Nodes);
						n1Index = getSmallNodeIndex({ 0, mi[0], 0, mi[1] }, 3, body_b_Nodes);
					}
					else if (e == 5) {
						n0Index = getSmallNodeIndex({ 0, 0, mi[0], mi[1] }, 2, body_b_Nodes);
						n1Index = getSmallNodeIndex({ 0, 0, mi[0], mi[1] }, 3, body_b_Nodes);
					}
					uint edgeIndex = mesh.findEdge(n0Index, n1Index);
					const Buffer<uint> smallBodies = mesh.getEdgeBodies(edgeIndex);
					Vector3 n0 = mesh.getNodePosition3(n0Index);
					Vector3 n1 = mesh.getNodePosition3(n1Index);
					uint sb = 0;
					while (!isInsideTetrahedron(mesh.getBodyAverage3(smallBodies[sb]), p0, p1, p2, p3, p0OppositeNormal, p1OppositeNormal, p2OppositeNormal, p3OppositeNormal))
						++sb;
					const Buffer<uint> smallFaces = mesh.getEdgeFaces(edgeIndex).getIntersection(mesh.getBodyFaces(smallBodies[sb]));
					Vector3 sf0_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[0]) : mesh.getFaceAverage3(smallFaces[0]));
					Vector3 sf1_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[1]) : mesh.getFaceAverage3(smallFaces[1]));
					Vector3 sb_bc = (circumcentric ? mesh.getBodyPosition3(smallBodies[sb]) : mesh.getBodyAverage3(smallBodies[sb]));
					Vector3 se_bc = 0.5 * (n0 + n1);
					Vector3 normalVector0 = 0.5 * TwoVector3(sf0_bc - se_bc, sb_bc - se_bc).dual();
					if (normalVector0.dot(n1 - n0) < 0)
						normalVector0 *= -1.0;
					if (minkowskiMetric)
						normalVector0.z = -normalVector0.z;
					Vector3 normalVector1 = 0.5 * TwoVector3(sf1_bc - se_bc, sb_bc - se_bc).dual();
					if (normalVector1.dot(n1 - n0) < 0)
						normalVector1 *= -1.0;
					if (minkowskiMetric)
						normalVector1.z = -normalVector1.z;
					dualFaceList[dfIndex].push_back({ sf0_bc, sb_bc, se_bc, normalVector0 });
					dualFaceList[dfIndex].push_back({ sf1_bc, sb_bc, se_bc, normalVector1 });
					++dfIndex;
				}
			}
			for (uint f = 0; f < 4; ++f) { //dual faces of small edges on faces
				for (uint e = 0; e < 3; ++e) {
					for (uint i = 0; i < faceEdgeHoles.size(); ++i) {
						mi_t mi;
						uint n0Index, n1Index;
						if (f == 0) {
							if (e == 0) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 0, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 1, body_b_Nodes);
							}
							else if (e == 1) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 0, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 2, body_b_Nodes);
							}
							else if (e == 2) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 1, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 2, body_b_Nodes);
							}
						}
						else if (f == 1) {
							if (e == 0) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 0, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 1, body_b_Nodes);
							}
							else if (e == 1) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 0, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 3, body_b_Nodes);
							}
							else if (e == 2) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 1, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 3, body_b_Nodes);
							}
						}
						else if (f == 2) {
							if (e == 0) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
								n0Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 0, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 2, body_b_Nodes);
							}
							else if (e == 1) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
								n0Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 0, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 3, body_b_Nodes);
							}
							else if (e == 2) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
								n0Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 2, body_b_Nodes);
								n1Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 3, body_b_Nodes);
							}
						}
						else if (f == 3) {
							if (e == 0) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
								n0Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 1, body_b_Nodes);
								n1Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 2, body_b_Nodes);
							}
							else if (e == 1) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
								n0Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 1, body_b_Nodes);
								n1Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 3, body_b_Nodes);
							}
							else if (e == 2) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
								n0Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 2, body_b_Nodes);
								n1Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 3, body_b_Nodes);
							}
						}
						uint edgeIndex = mesh.findEdge(n0Index, n1Index);
						const Buffer<uint> smallBodies = mesh.getEdgeBodies(edgeIndex);
						Vector3 n0 = mesh.getNodePosition3(n0Index);
						Vector3 n1 = mesh.getNodePosition3(n1Index);
						for (uint sb = 0; sb < smallBodies.size(); ++sb) {
							if (!isInsideTetrahedron(mesh.getBodyAverage3(smallBodies[sb]), p0, p1, p2, p3, p0OppositeNormal, p1OppositeNormal, p2OppositeNormal, p3OppositeNormal))
								continue;
							const Buffer<uint> smallFaces = mesh.getEdgeFaces(edgeIndex).getIntersection(mesh.getBodyFaces(smallBodies[sb]));
							Vector3 sf0_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[0]) : mesh.getFaceAverage3(smallFaces[0]));
							Vector3 sf1_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[1]) : mesh.getFaceAverage3(smallFaces[1]));
							Vector3 sb_bc = (circumcentric ? mesh.getBodyPosition3(smallBodies[sb]) : mesh.getBodyAverage3(smallBodies[sb]));
							Vector3 se_bc = 0.5 * (n0 + n1);
							Vector3 normalVector0 = 0.5 * TwoVector3(sf0_bc - se_bc, sb_bc - se_bc).dual();
							if (normalVector0.dot(n1 - n0) < 0)
								normalVector0 *= -1.0;
							if (minkowskiMetric)
								normalVector0.z = -normalVector0.z;
							Vector3 normalVector1 = 0.5 * TwoVector3(sf1_bc - se_bc, sb_bc - se_bc).dual();
							if (normalVector1.dot(n1 - n0) < 0)
								normalVector1 *= -1.0;
							if (minkowskiMetric)
								normalVector1.z = -normalVector1.z;
							dualFaceList[dfIndex].push_back({ sf0_bc, sb_bc, se_bc, normalVector0 });
							dualFaceList[dfIndex].push_back({ sf1_bc, sb_bc, se_bc, normalVector1 });
						}
						++dfIndex;
					}
				}
			}
			for (uint e = 0; e < 6; ++e) { //dual faces of small edges in the interior
				for (uint i = 0; i < bodyEdgeHoles.size(); ++i) {
					mi_t mi;
					uint n0Index;
					uint n1Index;
					if (e == 0) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 0);
						n0Index = getSmallNodeIndex(mi, 0, body_b_Nodes);
						n1Index = getSmallNodeIndex(mi, 1, body_b_Nodes);
					}
					else if (e == 1) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 1);
						n0Index = getSmallNodeIndex(mi, 0, body_b_Nodes);
						n1Index = getSmallNodeIndex(mi, 2, body_b_Nodes);
					}
					else if (e == 2) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 2);
						n0Index = getSmallNodeIndex(mi, 0, body_b_Nodes);
						n1Index = getSmallNodeIndex(mi, 3, body_b_Nodes);
					}
					else if (e == 3) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 3);
						n0Index = getSmallNodeIndex(mi, 1, body_b_Nodes);
						n1Index = getSmallNodeIndex(mi, 2, body_b_Nodes);
					}
					else if (e == 4) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 4);
						n0Index = getSmallNodeIndex(mi, 1, body_b_Nodes);
						n1Index = getSmallNodeIndex(mi, 3, body_b_Nodes);
					}
					else if (e == 5) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 5);
						n0Index = getSmallNodeIndex(mi, 2, body_b_Nodes);
						n1Index = getSmallNodeIndex(mi, 3, body_b_Nodes);
					}
					uint edgeIndex = mesh.findEdge(n0Index, n1Index);
					const Buffer<uint> smallBodies = mesh.getEdgeBodies(edgeIndex);
					Vector3 n0 = mesh.getNodePosition3(n0Index);
					Vector3 n1 = mesh.getNodePosition3(n1Index);
					for (uint sb = 0; sb < smallBodies.size(); ++sb) {
						const Buffer<uint> smallFaces = mesh.getEdgeFaces(edgeIndex).getIntersection(mesh.getBodyFaces(smallBodies[sb]));
						Vector3 sf0_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[0]) : mesh.getFaceAverage3(smallFaces[0]));
						Vector3 sf1_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[1]) : mesh.getFaceAverage3(smallFaces[1]));
						Vector3 sb_bc = (circumcentric ? mesh.getBodyPosition3(smallBodies[sb]) : mesh.getBodyAverage3(smallBodies[sb]));
						Vector3 se_bc = 0.5 * (n0 + n1);
						Vector3 normalVector0 = 0.5 * TwoVector3(sf0_bc - se_bc, sb_bc - se_bc).dual();
						if (normalVector0.dot(n1 - n0) < 0)
							normalVector0 *= -1.0;
						if (minkowskiMetric)
							normalVector0.z = -normalVector0.z;
						Vector3 normalVector1 = 0.5 * TwoVector3(sf1_bc - se_bc, sb_bc - se_bc).dual();
						if (normalVector1.dot(n1 - n0) < 0)
							normalVector1 *= -1.0;
						if (minkowskiMetric)
							normalVector1.z = -normalVector1.z;
						dualFaceList[dfIndex].push_back({ sf0_bc, sb_bc, se_bc, normalVector0 });
						dualFaceList[dfIndex].push_back({ sf1_bc, sb_bc, se_bc, normalVector1 });
					}
					++dfIndex;
				}
			}
			for (uint i = firstEdgeOfBody(b) + 6 * bodyEdgeHoles.size(); i < firstEdgeOfBody(b) + 6 * bodyEdgeHoles.size() + bodyFaceHoles.size(); ++i) { //dual faces of extra small edges
				const Buffer<uint> smallBodies = mesh.getEdgeBodies(i);
				const Buffer<uint>& smallNodes = mesh.getEdgeNodes(i);
				Vector3 n0 = mesh.getNodePosition3(smallNodes[0]);
				Vector3 n1 = mesh.getNodePosition3(smallNodes[1]);
				for (uint sb = 0; sb < smallBodies.size(); ++sb) {
					const Buffer<uint> smallFaces = mesh.getEdgeFaces(i).getIntersection(mesh.getBodyFaces(smallBodies[sb]));
					Vector3 sf0_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[0]) : mesh.getFaceAverage3(smallFaces[0]));
					Vector3 sf1_bc = (circumcentric ? mesh.getFacePosition3(smallFaces[1]) : mesh.getFaceAverage3(smallFaces[1]));
					Vector3 sb_bc = (circumcentric ? mesh.getBodyPosition3(smallBodies[sb]) : mesh.getBodyAverage3(smallBodies[sb]));
					Vector3 se_bc = 0.5 * (n0 + n1);
					Vector3 normalVector0 = 0.5 * TwoVector3(sf0_bc - se_bc, sb_bc - se_bc).dual();
					if (normalVector0.dot(n1 - n0) < 0)
						normalVector0 *= -1.0;
					if (minkowskiMetric)
						normalVector0.z = -normalVector0.z;
					Vector3 normalVector1 = 0.5 * TwoVector3(sf1_bc - se_bc, sb_bc - se_bc).dual();
					if (normalVector1.dot(n1 - n0) < 0)
						normalVector1 *= -1.0;
					if (minkowskiMetric)
						normalVector1.z = -normalVector1.z;
					dualFaceList[dfIndex].push_back({ sf0_bc, sb_bc, se_bc, normalVector0 });
					dualFaceList[dfIndex].push_back({ sf1_bc, sb_bc, se_bc, normalVector1 });
				}
				++dfIndex;
			}

			//compute the integrals of basis functions
			for (uint i = 0; i < dualFaceCount; ++i) {
				hodgeIntegralsBodyEdges[i].toVectorN(activeSmallEdgesInBodies.size());
				hodgeIntegralsFace0Edges[i].resize(facePermutations.size());
				hodgeIntegralsFace1Edges[i].resize(facePermutations.size());
				hodgeIntegralsFace2Edges[i].resize(facePermutations.size());
				hodgeIntegralsFace3Edges[i].resize(facePermutations.size());
				for (uint j = 0; j < facePermutations.size(); ++j) {
					hodgeIntegralsFace0Edges[i][j].toVectorN(activeSmallEdgesInFaces.size());
					hodgeIntegralsFace1Edges[i][j].toVectorN(activeSmallEdgesInFaces.size());
					hodgeIntegralsFace2Edges[i][j].toVectorN(activeSmallEdgesInFaces.size());
					hodgeIntegralsFace3Edges[i][j].toVectorN(activeSmallEdgesInFaces.size());
				}
				hodgeIntegralsEdge0Edges[i].resize(edgePermutations.size());
				hodgeIntegralsEdge1Edges[i].resize(edgePermutations.size());
				hodgeIntegralsEdge2Edges[i].resize(edgePermutations.size());
				hodgeIntegralsEdge3Edges[i].resize(edgePermutations.size());
				hodgeIntegralsEdge4Edges[i].resize(edgePermutations.size());
				hodgeIntegralsEdge5Edges[i].resize(edgePermutations.size());
				for (uint j = 0; j < edgePermutations.size(); ++j) {
					hodgeIntegralsEdge0Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
					hodgeIntegralsEdge1Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
					hodgeIntegralsEdge2Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
					hodgeIntegralsEdge3Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
					hodgeIntegralsEdge4Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
					hodgeIntegralsEdge5Edges[i][j].toVectorN(activeSmallEdgesInEdges.size());
				}
			}

			//first integrate body basis functions
			VectorN bodyCochain(activeSmallEdgesInBodies.size(), 0.0);
			VectorN bodyCoefficients(activeSmallEdgesInBodies.size(), 0.0);
			for (uint i = 0; i < activeSmallEdgesInBodies.size(); ++i) {
				bodyCochain[i] = 1.0;
				solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyCochain, bodyCoefficients);
				for (uint j = 0; j < dualFaceCount; ++j) {
					for (uint k = 0; k < dualFaceList[j].size(); ++k) {
						std::function<double(Vector3)> interpolant{
							[&](Vector3 p) -> double {
								return evaluate1FormWithBodyCoefficients(p, bodyCoefficients, body_b_Nodes, body_b_Edges, body_b_Faces,
									gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[j][k].normalVector);
							}
						};
						hodgeIntegralsBodyEdges[j][i] += gfd::integralAverage(interpolant, dualFaceList[j][k].p0, dualFaceList[j][k].p1, dualFaceList[j][k].p2);
					}

				}
				bodyCoefficients.val.fill(0.0);
				bodyCochain[i] = 0.0;
			}

			//next integrate face basis functions for each face and orientation
			VectorN faceCochain(activeSmallEdgesInFaces.size(), 0.0);
			VectorN faceCoefficients(activeSmallEdgesInFaces.size(), 0.0);
			VectorN bodyValues(activeSmallEdgesInBodies.size(), 0.0);
			for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
				faceCochain[i] = 1.0;
				solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceCochain, faceCoefficients);
				for (uint j = 0; j < facePermutations.size(); ++j) {
					//face0
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesFace0[j][k].dot(faceCoefficients);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithFaceCoefficients(p, faceCoefficients, 0, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsFace0Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsFace0Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
					//face1
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesFace1[j][k].dot(faceCoefficients);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithFaceCoefficients(p, faceCoefficients, 1, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsFace1Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsFace1Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
					//face2
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesFace2[j][k].dot(faceCoefficients);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithFaceCoefficients(p, faceCoefficients, 2, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsFace2Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsFace2Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
					//face3
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesFace3[j][k].dot(faceCoefficients);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithFaceCoefficients(p, faceCoefficients, 3, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsFace3Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsFace3Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
				}
				faceCoefficients.val.fill(0.0);
				faceCochain[i] = 0.0;
			}

			//for edge basis functions, the contribution from faces is excluded from the integrals (easier to take into account later)
			VectorN edgeCochain(activeSmallEdgesInEdges.size(), 0.0);
			VectorN faceValues(activeSmallEdgesInBodies.size(), 0.0);
			for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
				edgeCochain[i] = 1.0;
				solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeCochain, edgeCoefficients[i]);
				for (uint j = 0; j < edgePermutations.size(); ++j) {
					//edge0
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesEdge0[j][k].dot(edgeCoefficients[i]);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients[i], 0, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge0Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsEdge0Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
					//edge1
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesEdge1[j][k].dot(edgeCoefficients[i]);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients[i], 1, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge1Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsEdge1Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
					//edge2
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesEdge2[j][k].dot(edgeCoefficients[i]);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients[i], 2, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge2Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsEdge2Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
					//edge3
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesEdge3[j][k].dot(edgeCoefficients[i]);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients[i], 3, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge3Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsEdge3Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
					//edge4
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesEdge4[j][k].dot(edgeCoefficients[i]);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients[i], 4, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge4Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsEdge4Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
					//edge5
					for (uint k = 0; k < activeSmallEdgesInBodies.size(); ++k) {
						bodyValues[k] += bodyEdgeValuesEdge5[j][k].dot(edgeCoefficients[i]);
					}
					for (uint k = 0; k < dualFaceCount; ++k) {
						for (uint l = 0; l < dualFaceList[k].size(); ++l) {
							std::function<double(Vector3)> interpolant{
								[&](Vector3 p) -> double {
									return evaluate1FormWithEdgeCoefficients(p, edgeCoefficients[i], 5, j, body_b_Nodes, body_b_Edges, body_b_Faces,
										gradLambda0, gradLambda1, gradLambda2, gradLambda3).dot(dualFaceList[k][l].normalVector);
								}
							};
							hodgeIntegralsEdge5Edges[k][j][i] += gfd::integralAverage(interpolant, dualFaceList[k][l].p0, dualFaceList[k][l].p1, dualFaceList[k][l].p2);
						}
						hodgeIntegralsEdge5Edges[k][j][i] += bodyValues.dot(hodgeIntegralsBodyEdges[k]);
					}
					bodyValues.val.fill(0.0);
				}
				edgeCochain[i] = 0.0;
			}
		}

		if (!loadIntegralsFromFile) { //save the integrals to file
			saveHodgeIntegrals1Forms(foldername, circumcentric, scalingFactor, hodgeIntegralsBodyEdges, hodgeIntegralsFace0Edges, hodgeIntegralsFace1Edges, hodgeIntegralsFace2Edges,
				hodgeIntegralsFace3Edges, hodgeIntegralsEdge0Edges, hodgeIntegralsEdge1Edges, hodgeIntegralsEdge2Edges, hodgeIntegralsEdge3Edges, hodgeIntegralsEdge4Edges,
				hodgeIntegralsEdge5Edges, ref);
		}

		//form Hodge matrix from the precomputed integrals
		for (uint body = 0; body < refElements[ref].size(); ++body) {
			uint b = refElements[ref][body];
			//find the subsimplices and their orientations
			const Buffer<uint> nodes = getBodyNodesCustom(b);
			Buffer<uint> edges(6);
			edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
			edges[1] = mesh_old.findEdge(nodes[0], nodes[2]);
			edges[2] = mesh_old.findEdge(nodes[0], nodes[3]);
			edges[3] = mesh_old.findEdge(nodes[1], nodes[2]);
			edges[4] = mesh_old.findEdge(nodes[1], nodes[3]);
			edges[5] = mesh_old.findEdge(nodes[2], nodes[3]);
			Buffer<uint> faces(4);
			faces[0] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh_old));
			faces[1] = mesh_old.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh_old));
			faces[2] = mesh_old.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh_old));
			faces[3] = mesh_old.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh_old));
			Buffer<uint> faceNodes0 = mesh_old.getFaceNodes(faces[0]);
			Buffer<uint> faceNodes1 = mesh_old.getFaceNodes(faces[1]);
			Buffer<uint> faceNodes2 = mesh_old.getFaceNodes(faces[2]);
			Buffer<uint> faceNodes3 = mesh_old.getFaceNodes(faces[3]);
			Buffer<uint> faceEdges0(3);
			faceEdges0[0] = mesh_old.findEdge(faceNodes0[0], faceNodes0[1]);
			faceEdges0[1] = mesh_old.findEdge(faceNodes0[0], faceNodes0[2]);
			faceEdges0[2] = mesh_old.findEdge(faceNodes0[1], faceNodes0[2]);
			uint permutationIndexEdge0Face0 = getPermutation(mesh_old.getEdgeNodes(faceEdges0[0]), { faceNodes0[0], faceNodes0[1] });
			uint permutationIndexEdge1Face0 = getPermutation(mesh_old.getEdgeNodes(faceEdges0[1]), { faceNodes0[0], faceNodes0[2] });
			uint permutationIndexEdge2Face0 = getPermutation(mesh_old.getEdgeNodes(faceEdges0[2]), { faceNodes0[1], faceNodes0[2] });
			Buffer<uint> faceEdges1(3);
			faceEdges1[0] = mesh_old.findEdge(faceNodes1[0], faceNodes1[1]);
			faceEdges1[1] = mesh_old.findEdge(faceNodes1[0], faceNodes1[2]);
			faceEdges1[2] = mesh_old.findEdge(faceNodes1[1], faceNodes1[2]);
			uint permutationIndexEdge0Face1 = getPermutation(mesh_old.getEdgeNodes(faceEdges1[0]), { faceNodes1[0], faceNodes1[1] });
			uint permutationIndexEdge1Face1 = getPermutation(mesh_old.getEdgeNodes(faceEdges1[1]), { faceNodes1[0], faceNodes1[2] });
			uint permutationIndexEdge2Face1 = getPermutation(mesh_old.getEdgeNodes(faceEdges1[2]), { faceNodes1[1], faceNodes1[2] });
			Buffer<uint> faceEdges2(3);
			faceEdges2[0] = mesh_old.findEdge(faceNodes2[0], faceNodes2[1]);
			faceEdges2[1] = mesh_old.findEdge(faceNodes2[0], faceNodes2[2]);
			faceEdges2[2] = mesh_old.findEdge(faceNodes2[1], faceNodes2[2]);
			uint permutationIndexEdge0Face2 = getPermutation(mesh_old.getEdgeNodes(faceEdges2[0]), { faceNodes2[0], faceNodes2[1] });
			uint permutationIndexEdge1Face2 = getPermutation(mesh_old.getEdgeNodes(faceEdges2[1]), { faceNodes2[0], faceNodes2[2] });
			uint permutationIndexEdge2Face2 = getPermutation(mesh_old.getEdgeNodes(faceEdges2[2]), { faceNodes2[1], faceNodes2[2] });
			Buffer<uint> faceEdges3(3);
			faceEdges3[0] = mesh_old.findEdge(faceNodes3[0], faceNodes3[1]);
			faceEdges3[1] = mesh_old.findEdge(faceNodes3[0], faceNodes3[2]);
			faceEdges3[2] = mesh_old.findEdge(faceNodes3[1], faceNodes3[2]);
			uint permutationIndexEdge0Face3 = getPermutation(mesh_old.getEdgeNodes(faceEdges3[0]), { faceNodes3[0], faceNodes3[1] });
			uint permutationIndexEdge1Face3 = getPermutation(mesh_old.getEdgeNodes(faceEdges3[1]), { faceNodes3[0], faceNodes3[2] });
			uint permutationIndexEdge2Face3 = getPermutation(mesh_old.getEdgeNodes(faceEdges3[2]), { faceNodes3[1], faceNodes3[2] });
			uint permutationEdge0 = getPermutation(mesh_old.getEdgeNodes(edges[0]), { nodes[0], nodes[1] });
			uint permutationEdge1 = getPermutation(mesh_old.getEdgeNodes(edges[1]), { nodes[0], nodes[2] });
			uint permutationEdge2 = getPermutation(mesh_old.getEdgeNodes(edges[2]), { nodes[0], nodes[3] });
			uint permutationEdge3 = getPermutation(mesh_old.getEdgeNodes(edges[3]), { nodes[1], nodes[2] });
			uint permutationEdge4 = getPermutation(mesh_old.getEdgeNodes(edges[4]), { nodes[1], nodes[3] });
			uint permutationEdge5 = getPermutation(mesh_old.getEdgeNodes(edges[5]), { nodes[2], nodes[3] });
			uint permutationFace0 = getPermutation(mesh_old.getFaceNodes(faces[0]), { nodes[0], nodes[1], nodes[2] });
			uint permutationFace1 = getPermutation(mesh_old.getFaceNodes(faces[1]), { nodes[0], nodes[1], nodes[3] });
			uint permutationFace2 = getPermutation(mesh_old.getFaceNodes(faces[2]), { nodes[0], nodes[2], nodes[3] });
			uint permutationFace3 = getPermutation(mesh_old.getFaceNodes(faces[3]), { nodes[1], nodes[2], nodes[3] });

			//save node positions, outward normal vectors of the opposite faces, and the scaling factor in current body
			Vector3 p0 = mesh_old.getNodePosition3(nodes[0]);
			Vector3 p1 = mesh_old.getNodePosition3(nodes[1]);
			Vector3 p2 = mesh_old.getNodePosition3(nodes[2]);
			Vector3 p3 = mesh_old.getNodePosition3(nodes[3]);
			Vector3 p0OppositeNormal = TwoVector3(p2 - p1, p3 - p1).dual();
			if (p0OppositeNormal.dot(p1 - p0) < 0)
				p0OppositeNormal *= -1.0;
			Vector3 p1OppositeNormal = TwoVector3(p2 - p0, p3 - p0).dual();
			if (p1OppositeNormal.dot(p0 - p1) < 0)
				p1OppositeNormal *= -1.0;
			Vector3 p2OppositeNormal = TwoVector3(p1 - p0, p3 - p0).dual();
			if (p2OppositeNormal.dot(p0 - p2) < 0)
				p2OppositeNormal *= -1.0;
			Vector3 p3OppositeNormal = TwoVector3(p1 - p0, p2 - p0).dual();
			if (p3OppositeNormal.dot(p0 - p3) < 0)
				p3OppositeNormal *= -1.0;
			double scalingConstant = (p1 - p0).len() / scalingFactor;

			//find the indices of all simplices in current body and their orientations with respect to what was precomputed
			Buffer<uint> indices(dualFaceCount);
			Buffer<double> incidences(dualFaceCount);
			uint dfIndex = 0;
			for (uint e = 0; e < 6; ++e) { //small edges on edges
				for (uint i = 0; i < edgeMultiIndices.size(); ++i) {
					mi_t& mi = edgeMultiIndices[i];
					uint n0Index, n1Index;
					if (e == 0) {
						n0Index = getSmallNodeIndex({ mi[0], mi[1], 0, 0 }, 0, nodes);
						n1Index = getSmallNodeIndex({ mi[0], mi[1], 0, 0 }, 1, nodes);
					}
					else if (e == 1) {
						n0Index = getSmallNodeIndex({ mi[0], 0, mi[1], 0 }, 0, nodes);
						n1Index = getSmallNodeIndex({ mi[0], 0, mi[1], 0 }, 2, nodes);
					}
					else if (e == 2) {
						n0Index = getSmallNodeIndex({ mi[0], 0, 0, mi[1] }, 0, nodes);
						n1Index = getSmallNodeIndex({ mi[0], 0, 0, mi[1] }, 3, nodes);
					}
					else if (e == 3) {
						n0Index = getSmallNodeIndex({ 0, mi[0], mi[1], 0 }, 1, nodes);
						n1Index = getSmallNodeIndex({ 0, mi[0], mi[1], 0 }, 2, nodes);
					}
					else if (e == 4) {
						n0Index = getSmallNodeIndex({ 0, mi[0], 0, mi[1] }, 1, nodes);
						n1Index = getSmallNodeIndex({ 0, mi[0], 0, mi[1] }, 3, nodes);
					}
					else if (e == 5) {
						n0Index = getSmallNodeIndex({ 0, 0, mi[0], mi[1] }, 2, nodes);
						n1Index = getSmallNodeIndex({ 0, 0, mi[0], mi[1] }, 3, nodes);
					}
					uint edgeIndex = mesh.findEdge(n0Index, n1Index);
					indices[dfIndex] = edgeIndex;
					incidences[dfIndex] = mesh.getEdgeIncidence(edgeIndex, n1Index) * scalingConstant;
					++dfIndex;
				}
			}
			for (uint f = 0; f < 4; ++f) { //small edges on faces
				for (uint e = 0; e < 3; ++e) {
					for (uint i = 0; i < faceEdgeHoles.size(); ++i) {
						mi_t mi;
						uint n0Index, n1Index;
						if (f == 0) {
							if (e == 0) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 0, nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 1, nodes);
							}
							else if (e == 1) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 0, nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 2, nodes);
							}
							else if (e == 2) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 1, nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], mi[2], 0 }, 2, nodes);
							}
						}
						else if (f == 1) {
							if (e == 0) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 0, nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 1, nodes);
							}
							else if (e == 1) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 0, nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 3, nodes);
							}
							else if (e == 2) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
								n0Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 1, nodes);
								n1Index = getSmallNodeIndex({ mi[0], mi[1], 0, mi[2] }, 3, nodes);
							}
						}
						else if (f == 2) {
							if (e == 0) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
								n0Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 0, nodes);
								n1Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 2, nodes);
							}
							else if (e == 1) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
								n0Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 0, nodes);
								n1Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 3, nodes);
							}
							else if (e == 2) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
								n0Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 2, nodes);
								n1Index = getSmallNodeIndex({ mi[0], 0, mi[1], mi[2] }, 3, nodes);
							}
						}
						else if (f == 3) {
							if (e == 0) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 0);
								n0Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 1, nodes);
								n1Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 2, nodes);
							}
							else if (e == 1) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 1);
								n0Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 1, nodes);
								n1Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 3, nodes);
							}
							else if (e == 2) {
								mi = *getHoleEdgeMultiIndex(faceEdgeHoles[i], 2);
								n0Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 2, nodes);
								n1Index = getSmallNodeIndex({ 0, mi[0], mi[1], mi[2] }, 3, nodes);
							}
						}
						uint edgeIndex = mesh.findEdge(n0Index, n1Index);
						indices[dfIndex] = edgeIndex;
						incidences[dfIndex] = mesh.getEdgeIncidence(edgeIndex, n1Index) * scalingConstant;
						++dfIndex;
					}
				}
			}
			for (uint e = 0; e < 6; ++e) { //small edges in the interior
				for (uint i = 0; i < bodyEdgeHoles.size(); ++i) {
					mi_t mi;
					uint n0Index;
					uint n1Index;
					if (e == 0) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 0);
						n0Index = getSmallNodeIndex(mi, 0, nodes);
						n1Index = getSmallNodeIndex(mi, 1, nodes);
					}
					else if (e == 1) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 1);
						n0Index = getSmallNodeIndex(mi, 0, nodes);
						n1Index = getSmallNodeIndex(mi, 2, nodes);
					}
					else if (e == 2) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 2);
						n0Index = getSmallNodeIndex(mi, 0, nodes);
						n1Index = getSmallNodeIndex(mi, 3, nodes);
					}
					else if (e == 3) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 3);
						n0Index = getSmallNodeIndex(mi, 1, nodes);
						n1Index = getSmallNodeIndex(mi, 2, nodes);
					}
					else if (e == 4) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 4);
						n0Index = getSmallNodeIndex(mi, 1, nodes);
						n1Index = getSmallNodeIndex(mi, 3, nodes);
					}
					else if (e == 5) {
						mi = *getHoleEdgeMultiIndex(bodyEdgeHoles[i], 5);
						n0Index = getSmallNodeIndex(mi, 2, nodes);
						n1Index = getSmallNodeIndex(mi, 3, nodes);
					}
					uint edgeIndex = mesh.findEdge(n0Index, n1Index);
					indices[dfIndex] = edgeIndex;
					incidences[dfIndex] = mesh.getEdgeIncidence(edgeIndex, n1Index) * scalingConstant;
					++dfIndex;
				}
			}
			for (uint i = firstEdgeOfBody(b) + 6 * bodyEdgeHoles.size(); i < firstEdgeOfBody(b) + 6 * bodyEdgeHoles.size() + bodyFaceHoles.size(); ++i) { //extra small edges
				uint edgeIndex = i;
				indices[dfIndex] = edgeIndex;
				incidences[dfIndex] = scalingConstant;
				++dfIndex;
			}

			//get the integrals of all basis functions
			VectorN faceValues0(activeSmallEdgesInFaces.size(), 0.0);
			VectorN faceValues1(activeSmallEdgesInFaces.size(), 0.0);
			for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) { //edge basis functions
				uint edgeIndex;
				//edge0 (belongs to face0 and face1)
				if (edges[0] == faceEdges0[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face0][j].dot(edgeCoefficients[i]);
					}
				else if (edges[0] == faceEdges0[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face0][j].dot(edgeCoefficients[i]);
					}
				else if (edges[0] == faceEdges0[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face0][j].dot(edgeCoefficients[i]);
					}
				if (edges[0] == faceEdges1[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face1][j].dot(edgeCoefficients[i]);
					}
				else if (edges[0] == faceEdges1[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face1][j].dot(edgeCoefficients[i]);
					}
				else if (edges[0] == faceEdges1[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face1][j].dot(edgeCoefficients[i]);
					}
				edgeIndex = smallEdgesEdgeList(edges[0], i);
				for (uint j = 0; j < dualFaceCount; ++j) {
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsEdge0Edges[j][permutationEdge0][i];
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace0Edges[j][permutationFace0].dot(faceValues0);
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace1Edges[j][permutationFace1].dot(faceValues1);
				}
				faceValues0.val.fill(0.0);
				faceValues1.val.fill(0.0);
				//edge1 (belongs to face0 and face2)
				if (edges[1] == faceEdges0[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face0][j].dot(edgeCoefficients[i]);
					}
				else if (edges[1] == faceEdges0[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face0][j].dot(edgeCoefficients[i]);
					}
				else if (edges[1] == faceEdges0[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face0][j].dot(edgeCoefficients[i]);
					}
				if (edges[1] == faceEdges2[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face2][j].dot(edgeCoefficients[i]);
					}
				else if (edges[1] == faceEdges2[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face2][j].dot(edgeCoefficients[i]);
					}
				else if (edges[1] == faceEdges2[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face2][j].dot(edgeCoefficients[i]);
					}
				edgeIndex = smallEdgesEdgeList(edges[1], i);
				for (uint j = 0; j < dualFaceCount; ++j) {
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsEdge1Edges[j][permutationEdge1][i];
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace0Edges[j][permutationFace0].dot(faceValues0);
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace2Edges[j][permutationFace2].dot(faceValues1);
				}
				faceValues0.val.fill(0.0);
				faceValues1.val.fill(0.0);
				//edge2 (belongs to face1 and face2)
				if (edges[2] == faceEdges1[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face1][j].dot(edgeCoefficients[i]);
					}
				else if (edges[2] == faceEdges1[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face1][j].dot(edgeCoefficients[i]);
					}
				else if (edges[2] == faceEdges1[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face1][j].dot(edgeCoefficients[i]);
					}
				if (edges[2] == faceEdges2[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face2][j].dot(edgeCoefficients[i]);
					}
				else if (edges[2] == faceEdges2[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face2][j].dot(edgeCoefficients[i]);
					}
				else if (edges[2] == faceEdges2[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face2][j].dot(edgeCoefficients[i]);
					}
				edgeIndex = smallEdgesEdgeList(edges[2], i);
				for (uint j = 0; j < dualFaceCount; ++j) {
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsEdge2Edges[j][permutationEdge2][i];
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace1Edges[j][permutationFace1].dot(faceValues0);
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace2Edges[j][permutationFace2].dot(faceValues1);
				}
				faceValues0.val.fill(0.0);
				faceValues1.val.fill(0.0);
				//edge3 (belongs to face0 and face3)
				if (edges[3] == faceEdges0[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face0][j].dot(edgeCoefficients[i]);
					}
				else if (edges[3] == faceEdges0[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face0][j].dot(edgeCoefficients[i]);
					}
				else if (edges[3] == faceEdges0[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face0][j].dot(edgeCoefficients[i]);
					}
				if (edges[3] == faceEdges3[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face3][j].dot(edgeCoefficients[i]);
					}
				else if (edges[3] == faceEdges3[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face3][j].dot(edgeCoefficients[i]);
					}
				else if (edges[3] == faceEdges3[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face3][j].dot(edgeCoefficients[i]);
					}
				edgeIndex = smallEdgesEdgeList(edges[3], i);
				for (uint j = 0; j < dualFaceCount; ++j) {
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsEdge3Edges[j][permutationEdge3][i];
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace0Edges[j][permutationFace0].dot(faceValues0);
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace3Edges[j][permutationFace3].dot(faceValues1);
				}
				faceValues0.val.fill(0.0);
				faceValues1.val.fill(0.0);
				//edge4 (belongs to face1 and face3)
				if (edges[4] == faceEdges1[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face1][j].dot(edgeCoefficients[i]);
					}
				else if (edges[4] == faceEdges1[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face1][j].dot(edgeCoefficients[i]);
					}
				else if (edges[4] == faceEdges1[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face1][j].dot(edgeCoefficients[i]);
					}
				if (edges[4] == faceEdges3[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face3][j].dot(edgeCoefficients[i]);
					}
				else if (edges[4] == faceEdges3[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face3][j].dot(edgeCoefficients[i]);
					}
				else if (edges[4] == faceEdges3[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face3][j].dot(edgeCoefficients[i]);
					}
				edgeIndex = smallEdgesEdgeList(edges[4], i);
				for (uint j = 0; j < dualFaceCount; ++j) {
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsEdge4Edges[j][permutationEdge4][i];
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace1Edges[j][permutationFace1].dot(faceValues0);
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace3Edges[j][permutationFace3].dot(faceValues1);
				}
				faceValues0.val.fill(0.0);
				faceValues1.val.fill(0.0);
				//edge5 (belongs to face2 and face3)
				if (edges[5] == faceEdges2[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face2][j].dot(edgeCoefficients[i]);
					}
				else if (edges[5] == faceEdges2[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face2][j].dot(edgeCoefficients[i]);
					}
				else if (edges[5] == faceEdges2[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues0[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face2][j].dot(edgeCoefficients[i]);
					}
				if (edges[5] == faceEdges3[0])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge0[permutationIndexEdge0Face3][j].dot(edgeCoefficients[i]);
					}
				else if (edges[5] == faceEdges3[1])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge1[permutationIndexEdge1Face3][j].dot(edgeCoefficients[i]);
					}
				else if (edges[5] == faceEdges3[2])
					for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
						faceValues1[j] += faceEdgeValuesEdge2[permutationIndexEdge2Face3][j].dot(edgeCoefficients[i]);
					}
				edgeIndex = smallEdgesEdgeList(edges[5], i);
				for (uint j = 0; j < dualFaceCount; ++j) {
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsEdge5Edges[j][permutationEdge5][i];
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace2Edges[j][permutationFace2].dot(faceValues0);
					star[indices[j]][edgeIndex] += incidences[j] * hodgeIntegralsFace3Edges[j][permutationFace3].dot(faceValues1);
				}
				faceValues0.val.fill(0.0);
				faceValues1.val.fill(0.0);
			}
			for (uint i = 0; i < dualFaceCount; ++i) { //face basis functions
				for (uint j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
					star[indices[i]][smallEdgesFaceList(faces[0], j)] += incidences[i] * hodgeIntegralsFace0Edges[i][permutationFace0][j];
					star[indices[i]][smallEdgesFaceList(faces[1], j)] += incidences[i] * hodgeIntegralsFace1Edges[i][permutationFace1][j];
					star[indices[i]][smallEdgesFaceList(faces[2], j)] += incidences[i] * hodgeIntegralsFace2Edges[i][permutationFace2][j];
					star[indices[i]][smallEdgesFaceList(faces[3], j)] += incidences[i] * hodgeIntegralsFace3Edges[i][permutationFace3][j];
				}
			}
			for (uint i = 0; i < dualFaceCount; ++i) { //body basis functions
				for (uint j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
					star[indices[i]][smallEdgesBodyList(b, j)] += incidences[i] * hodgeIntegralsBodyEdges[i][j];
				}
			}
		}
	}
}

/*
The following functions return the index of the active small simplex j of the big simplex i.
Assumes that smallEdgesFaceOffsets have been stored in face 0 and smallEdgesBodyOffsets and smallFacesBodyOffsets in body 0 of the initial mesh when it was refined.
Then the indices of small simplices can be recovered from those of the big simplices and need not be explicitly saved in memory.
*/

uint SmallSimplexPartition3D::smallNodesEdgeList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + i * activeSmallNodesInEdges.size() + j;
}
uint SmallSimplexPartition3D::smallEdgesEdgeList(uint i, uint j) const {
	return i * edgeMultiIndices.size() + j;
}
uint SmallSimplexPartition3D::smallNodesFaceList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + mesh_old_ptr->getEdgeSize() * activeSmallNodesInEdges.size() + i * activeSmallNodesInFaces.size() + j;
}
uint SmallSimplexPartition3D::smallEdgesFaceList(uint i, uint j) const {
	return i * 3 * faceEdgeHoles.size() + smallEdgesFaceOffsets[j];
}
uint SmallSimplexPartition3D::smallFacesFaceList(uint i, uint j) const {
	return i * (faceMultiIndices.size() + faceEdgeHoles.size()) + j;
}
uint SmallSimplexPartition3D::smallNodesBodyList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + mesh_old_ptr->getEdgeSize() * activeSmallNodesInEdges.size() + mesh_old_ptr->getFaceSize() * activeSmallNodesInFaces.size()
		+ i * activeSmallNodesInBodies.size() + j;
}
uint SmallSimplexPartition3D::smallEdgesBodyList(uint i, uint j) const {
	return i * (6 * bodyEdgeHoles.size() + bodyFaceHoles.size()) + smallEdgesBodyOffsets[j];
}
uint SmallSimplexPartition3D::smallFacesBodyList(uint i, uint j) const {
	return i * (8 * bodyFaceHoles.size() + 4 * bodyEdgeHoles.size()) + smallFacesBodyOffsets[j];
}
uint SmallSimplexPartition3D::smallBodiesBodyList(uint i, uint j) const {
	return i * (bodyMultiIndices.size() + bodyEdgeHoles.size() + 4 * bodyFaceHoles.size()) + j;
}

/*
The following functions return the smallest index of small simplices (not necessarily active) in the interior of the given big simplex.
*/

uint SmallSimplexPartition3D::firstEdgeOfFace(uint i) const {
	return mesh_old_ptr->getEdgeSize() * edgeMultiIndices.size() + i * 3 * faceEdgeHoles.size();
}
uint SmallSimplexPartition3D::firstEdgeOfBody(uint i) const {
	return mesh_old_ptr->getEdgeSize() * edgeMultiIndices.size() + mesh_old_ptr->getFaceSize() * 3 * faceEdgeHoles.size() + i * (6 * bodyEdgeHoles.size() + bodyFaceHoles.size());
}
uint SmallSimplexPartition3D::firstFaceOfBody(uint i) const {
	return mesh_old_ptr->getFaceSize() * (faceMultiIndices.size() + faceEdgeHoles.size()) + i * (8 * bodyFaceHoles.size() + 4 * bodyEdgeHoles.size());
}

//returns the multi-index that maps the given edge to the given hole boundary
mi_t* SmallSimplexPartition3D::getHoleEdgeMultiIndex(const mi_t& holeIndex, uint edgeIndex) const {
	mi_t multiIndex{ holeIndex };
	uint i;
	if (multiIndex.size() == 3) { //hole is inverted triangle in 2D
		for (i = 0; i < multiIndex.size(); ++i) {
			if (faceEdges[edgeIndex][i] == 0)
				++multiIndex[i];
		}
		return &faceMultiIndices[indexOf(faceMultiIndices, multiIndex)];
	}
	else { //hole is inverted tetrahedron in 3D
		for (i = 0; i < multiIndex.size(); ++i) {
			if (bodyEdges[edgeIndex][i] == 0)
				++multiIndex[i];
		}
		return &bodyMultiIndices[indexOf(bodyMultiIndices, multiIndex)];
	}
}

//returns the multi-index that maps the given face to the given hole boundary
mi_t* SmallSimplexPartition3D::getHoleFaceMultiIndex(const mi_t& holeIndex, uint faceIndex) const {
	mi_t multiIndex{ holeIndex };
	uint i;
	for (i = 0; i < multiIndex.size(); ++i) {
		if (bodyFaces[faceIndex][i] == 0)
			++multiIndex[i];
	}
	return &bodyMultiIndices[indexOf(bodyMultiIndices, multiIndex)];
}

//returns the opposite edge (in a tetrahedron)
uint SmallSimplexPartition3D::getOppositeEdge(uint edge) const {
	//return 5 - edge;
	switch (edge) {
	case 0:
		return 5;
	case 1:
		return 4;
	case 2:
		return 3;
	case 3:
		return 2;
	case 4:
		return 1;
	case 5:
		return 0;
	}
}

//returns the index of the small node that is the image of multiIndex when we map the given node
uint SmallSimplexPartition3D::getSmallNodeIndex(mi_t multiIndex, uint nodeIndex, Buffer<uint> bigSimplexNodes) const {
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
	case 3:
		bigSimplexIndex = mesh_old.findFace(findEdges(bigSimplexNodes[0], bigSimplexNodes[1], bigSimplexNodes[2], mesh_old));
		nodesOrdered = mesh_old.getFaceNodes(bigSimplexIndex);
		break;
	default:
		bigSimplexIndex = mesh_old.findBody(findFaces(bigSimplexNodes[0], bigSimplexNodes[1], bigSimplexNodes[2], bigSimplexNodes[3], mesh_old));
		nodesOrdered = getBodyNodesCustom(bigSimplexIndex);
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
	case 4:
		return smallNodesBodyList(bigSimplexIndex, indexOf(activeSmallNodesInBodies, multiIndexPermuted));
	}
}

//computes the value resulting from the coefficients on the given edge in the element whose nodes, edges, and faces are given
Vector3 SmallSimplexPartition3D::evaluate1FormWithEdgeCoefficients(const Vector3& evaluationPoint, const VectorN& edgeCoefficients, uint edgeIndex, uint edgePermutation,
	const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
	const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const {
	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find barycentric coordinates
	double lambda0, lambda1, lambda2, lambda3;
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() * gradLambda0.len();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() * gradLambda1.len();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() * gradLambda2.len();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() * gradLambda3.len();

	//compute the values of the lowest order Whitney forms
	Vector3 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector3 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector3 w2 = lambda0 * gradLambda3 - lambda3 * gradLambda0;
	Vector3 w3 = lambda1 * gradLambda2 - lambda2 * gradLambda1;
	Vector3 w4 = lambda1 * gradLambda3 - lambda3 * gradLambda1;
	Vector3 w5 = lambda2 * gradLambda3 - lambda3 * gradLambda2;

	//compute the result using the coefficients given
	Buffer<double> d(6, 0.0);
	std::vector<double> barycentricCoords;

	switch (edgeIndex) {
	case 0:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[0] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[0] = -d[0];
		break;
	case 1:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[1] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[1] = -d[1];
		break;
	case 2:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda3 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[2] = -d[2];
		break;
	case 3:
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[3] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[3] = -d[3];
		break;
	case 4:
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda3 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[4] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[4] = -d[4];
		break;
	case 5:
		barycentricCoords = permute(std::vector<double>{ lambda2, lambda3 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[5] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[5] = -d[5];
		break;
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3 + d[4] * w4 + d[5] * w5;
}

//computes the value resulting from the given coefficients in the element whose nodes, edges, and faces are given
Vector3 SmallSimplexPartition3D::evaluate1FormWithCoefficients(const Vector3& evaluationPoint, const VectorN& edgeCoefficients, uint edgeIndex, uint edgePermutation, const VectorN& face0Coefficients,
	uint face0Index, uint face0Permutation, const VectorN& face1Coefficients, uint face1Index, uint face1Permutation, const VectorN& bodyCoefficients,
	const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
	const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const {

	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find barycentric coordinates
	double lambda0, lambda1, lambda2, lambda3;
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() * gradLambda0.len();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() * gradLambda1.len();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() * gradLambda2.len();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() * gradLambda3.len();

	//compute the values of the lowest order Whitney forms
	Vector3 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector3 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector3 w2 = lambda0 * gradLambda3 - lambda3 * gradLambda0;
	Vector3 w3 = lambda1 * gradLambda2 - lambda2 * gradLambda1;
	Vector3 w4 = lambda1 * gradLambda3 - lambda3 * gradLambda1;
	Vector3 w5 = lambda2 * gradLambda3 - lambda3 * gradLambda2;

	//compute the result using the coefficients given
	Buffer<double> d(6, 0.0);
	std::vector<double> barycentricCoords;
	Buffer<double> faceContributions(3, 0.0); //contributions to the edges of the faces

	switch (edgeIndex) {
	case 0:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[0] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[0] = -d[0];
		//face 0
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[face0Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face0Coefficients[i];
		}
		switch (face0Permutation) {
		case 0:
			d[0] += faceContributions[0];
			d[1] += faceContributions[1];
			d[3] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[1] += faceContributions[0];
			d[3] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[1] += faceContributions[2];
			d[3] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[1] += -faceContributions[2];
			d[3] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[1] += -faceContributions[0];
			d[3] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[1] += -faceContributions[1];
			d[3] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		//face1
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[face1Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face1Coefficients[i];
		}
		switch (face1Permutation) {
		case 0:
			d[0] += faceContributions[0];
			d[2] += faceContributions[1];
			d[4] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[2] += faceContributions[0];
			d[4] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[4] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[4] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[4] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[4] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		break;
	case 1:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[1] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[1] = -d[1];
		//face 0
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[face0Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face0Coefficients[i];
		}
		switch (face0Permutation) {
		case 0:
			d[0] += faceContributions[0];
			d[1] += faceContributions[1];
			d[3] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[1] += faceContributions[0];
			d[3] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[1] += faceContributions[2];
			d[3] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[1] += -faceContributions[2];
			d[3] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[1] += -faceContributions[0];
			d[3] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[1] += -faceContributions[1];
			d[3] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		//face2
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[face1Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face1Coefficients[i];
		}
		switch (face1Permutation) {
		case 0:
			d[1] += faceContributions[0];
			d[2] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[1] += faceContributions[1];
			d[2] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[1] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[1] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[1] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[1] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		break;
	case 2:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda3 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[2] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[2] = -d[2];
		//face1
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[face0Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face0Coefficients[i];
		}
		switch (face0Permutation) {
		case 0:
			d[0] += faceContributions[0];
			d[2] += faceContributions[1];
			d[4] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[2] += faceContributions[0];
			d[4] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[4] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[4] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[4] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[4] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		//face2
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[face1Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face1Coefficients[i];
		}
		switch (face1Permutation) {
		case 0:
			d[1] += faceContributions[0];
			d[2] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[1] += faceContributions[1];
			d[2] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[1] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[1] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[1] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[1] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		break;
	case 3:
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[3] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[3] = -d[3];
		//face 0
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[face0Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face0Coefficients[i];
		}
		switch (face0Permutation) {
		case 0:
			d[0] += faceContributions[0];
			d[1] += faceContributions[1];
			d[3] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[1] += faceContributions[0];
			d[3] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[1] += faceContributions[2];
			d[3] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[1] += -faceContributions[2];
			d[3] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[1] += -faceContributions[0];
			d[3] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[1] += -faceContributions[1];
			d[3] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		//face3
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[face1Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face1Coefficients[i];
		}
		switch (face1Permutation) {
		case 0:
			d[3] += faceContributions[0];
			d[4] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[3] += faceContributions[1];
			d[4] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[3] += -faceContributions[0];
			d[4] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[3] += -faceContributions[1];
			d[4] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[3] += faceContributions[2];
			d[4] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[3] += -faceContributions[2];
			d[4] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		break;
	case 4:
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda3 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[4] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[4] = -d[4];
		//face1
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[face0Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face0Coefficients[i];
		}
		switch (face0Permutation) {
		case 0:
			d[0] += faceContributions[0];
			d[2] += faceContributions[1];
			d[4] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[2] += faceContributions[0];
			d[4] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[4] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[4] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[4] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[4] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		//face3
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[face1Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face1Coefficients[i];
		}
		switch (face1Permutation) {
		case 0:
			d[3] += faceContributions[0];
			d[4] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[3] += faceContributions[1];
			d[4] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[3] += -faceContributions[0];
			d[4] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[3] += -faceContributions[1];
			d[4] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[3] += faceContributions[2];
			d[4] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[3] += -faceContributions[2];
			d[4] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		break;
	case 5:
		barycentricCoords = permute(std::vector<double>{ lambda2, lambda3 }, edgePermutations[edgePermutation]);
		for (uint i = 0; i < activeSmallEdgesInEdges.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInEdges[i].multiIndex;
			d[5] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1]) * edgeCoefficients[i];
		}
		if (getEdgePermutationSign(edgePermutation) < 0)
			d[5] = -d[5];
		//face2
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[face0Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face0Coefficients[i];
		}
		switch (face0Permutation) {
		case 0:
			d[1] += faceContributions[0];
			d[2] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[1] += faceContributions[1];
			d[2] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[1] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[1] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[1] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[1] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		//face3
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[face1Permutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * face1Coefficients[i];
		}
		switch (face1Permutation) {
		case 0:
			d[3] += faceContributions[0];
			d[4] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[3] += faceContributions[1];
			d[4] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[3] += -faceContributions[0];
			d[4] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[3] += -faceContributions[1];
			d[4] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[3] += faceContributions[2];
			d[4] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[3] += -faceContributions[2];
			d[4] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		faceContributions.fill(0.0);
		break;
	}

	//small simplices in the interior
	for (uint i = 0; i < activeSmallEdgesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInBodies[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInBodies[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3 + d[4] * w4 + d[5] * w5;
	
}

//computes the value resulting from the coefficients on the given face in the element whose nodes, edges, and faces are given
Vector3 SmallSimplexPartition3D::evaluate1FormWithFaceCoefficients(const Vector3& evaluationPoint, const VectorN& faceCoefficients, uint faceIndex, uint facePermutation, 
	const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
	const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const {

	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find barycentric coordinates
	double lambda0, lambda1, lambda2, lambda3;
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() * gradLambda0.len();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() * gradLambda1.len();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() * gradLambda2.len();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() * gradLambda3.len();

	//compute the values of the lowest order Whitney forms
	Vector3 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector3 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector3 w2 = lambda0 * gradLambda3 - lambda3 * gradLambda0;
	Vector3 w3 = lambda1 * gradLambda2 - lambda2 * gradLambda1;
	Vector3 w4 = lambda1 * gradLambda3 - lambda3 * gradLambda1;
	Vector3 w5 = lambda2 * gradLambda3 - lambda3 * gradLambda2;

	//compute the result using the coefficients given
	Buffer<double> d(6, 0.0);
	std::vector<double> barycentricCoords;
	Buffer<double> faceContributions(3, 0.0); //contributions to the edges of the faces

	switch (faceIndex) {
	case 0:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[facePermutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[i];
		}
		switch (facePermutation) {
		case 0:
			d[0] += faceContributions[0];
			d[1] += faceContributions[1];
			d[3] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[1] += faceContributions[0];
			d[3] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[1] += faceContributions[2];
			d[3] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[1] += -faceContributions[2];
			d[3] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[1] += -faceContributions[0];
			d[3] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[1] += -faceContributions[1];
			d[3] += -faceContributions[0];
			break;
		}
		break;
	case 1:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[facePermutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[i];
		}
		switch (facePermutation) {
		case 0:
			d[0] += faceContributions[0];
			d[2] += faceContributions[1];
			d[4] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[2] += faceContributions[0];
			d[4] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[4] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[4] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[4] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[4] += -faceContributions[0];
			break;
		}
		break;
	case 2:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[facePermutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[i];
		}
		switch (facePermutation) {
		case 0:
			d[1] += faceContributions[0];
			d[2] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[1] += faceContributions[1];
			d[2] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[1] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[1] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[1] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[1] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		break;
	case 3:
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[facePermutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[i];
		}
		switch (facePermutation) {
		case 0:
			d[3] += faceContributions[0];
			d[4] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[3] += faceContributions[1];
			d[4] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[3] += -faceContributions[0];
			d[4] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[3] += -faceContributions[1];
			d[4] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[3] += faceContributions[2];
			d[4] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[3] += -faceContributions[2];
			d[4] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		break;
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3 + d[4] * w4 + d[5] * w5;
}

//computes the value resulting from the given coefficients in the element whose nodes, edges, and faces are given
Vector3 SmallSimplexPartition3D::evaluate1FormWithCoefficients(const Vector3& evaluationPoint, const VectorN& faceCoefficients, uint faceIndex, uint facePermutation, const VectorN& bodyCoefficients,
	const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
	const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const {

	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find barycentric coordinates
	double lambda0, lambda1, lambda2, lambda3;
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() * gradLambda0.len();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() * gradLambda1.len();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() * gradLambda2.len();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() * gradLambda3.len();

	//compute the values of the lowest order Whitney forms
	Vector3 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector3 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector3 w2 = lambda0 * gradLambda3 - lambda3 * gradLambda0;
	Vector3 w3 = lambda1 * gradLambda2 - lambda2 * gradLambda1;
	Vector3 w4 = lambda1 * gradLambda3 - lambda3 * gradLambda1;
	Vector3 w5 = lambda2 * gradLambda3 - lambda3 * gradLambda2;

	//compute the result using the coefficients given
	Buffer<double> d(6, 0.0);
	std::vector<double> barycentricCoords;
	Buffer<double> faceContributions(3, 0.0); //contributions to the edges of the faces

	switch (faceIndex) {
	case 0:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda2 }, facePermutations[facePermutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[i];
		}
		switch (facePermutation) {
		case 0:
			d[0] += faceContributions[0];
			d[1] += faceContributions[1];
			d[3] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[1] += faceContributions[0];
			d[3] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[1] += faceContributions[2];
			d[3] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[1] += -faceContributions[2];
			d[3] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[1] += -faceContributions[0];
			d[3] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[1] += -faceContributions[1];
			d[3] += -faceContributions[0];
			break;
		}
		break;
	case 1:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda1, lambda3 }, facePermutations[facePermutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[i];
		}
		switch (facePermutation) {
		case 0:
			d[0] += faceContributions[0];
			d[2] += faceContributions[1];
			d[4] += faceContributions[2];
			break;
		case 1:
			d[0] += faceContributions[1];
			d[2] += faceContributions[0];
			d[4] += -faceContributions[2];
			break;
		case 2:
			d[0] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[4] += faceContributions[1];
			break;
		case 3:
			d[0] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[4] += faceContributions[0];
			break;
		case 4:
			d[0] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[4] += -faceContributions[1];
			break;
		case 5:
			d[0] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[4] += -faceContributions[0];
			break;
		}
		break;
	case 2:
		barycentricCoords = permute(std::vector<double>{ lambda0, lambda2, lambda3 }, facePermutations[facePermutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[i];
		}
		switch (facePermutation) {
		case 0:
			d[1] += faceContributions[0];
			d[2] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[1] += faceContributions[1];
			d[2] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[1] += -faceContributions[0];
			d[2] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[1] += -faceContributions[1];
			d[2] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[1] += faceContributions[2];
			d[2] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[1] += -faceContributions[2];
			d[2] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		break;
	case 3:
		barycentricCoords = permute(std::vector<double>{ lambda1, lambda2, lambda3 }, facePermutations[facePermutation]);
		for (uint i = 0; i < activeSmallEdgesInFaces.size(); ++i) {
			const mi_t& mi = *activeSmallEdgesInFaces[i].multiIndex;
			faceContributions[getEdgeIndex(*activeSmallEdgesInFaces[i].nodeIndices)] += std::pow(barycentricCoords[0], mi[0]) * std::pow(barycentricCoords[1], mi[1])
				* std::pow(barycentricCoords[2], mi[2]) * faceCoefficients[i];
		}
		switch (facePermutation) {
		case 0:
			d[3] += faceContributions[0];
			d[4] += faceContributions[1];
			d[5] += faceContributions[2];
			break;
		case 1:
			d[3] += faceContributions[1];
			d[4] += faceContributions[0];
			d[5] += -faceContributions[2];
			break;
		case 2:
			d[3] += -faceContributions[0];
			d[4] += faceContributions[2];
			d[5] += faceContributions[1];
			break;
		case 3:
			d[3] += -faceContributions[1];
			d[4] += -faceContributions[2];
			d[5] += faceContributions[0];
			break;
		case 4:
			d[3] += faceContributions[2];
			d[4] += -faceContributions[0];
			d[5] += -faceContributions[1];
			break;
		case 5:
			d[3] += -faceContributions[2];
			d[4] += -faceContributions[1];
			d[5] += -faceContributions[0];
			break;
		}
		break;
	}

	//small simplices in the interior
	for (uint i = 0; i < activeSmallEdgesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInBodies[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInBodies[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3 + d[4] * w4 + d[5] * w5;
}

//computes the value resulting from the given coefficients in the element whose nodes, edges, and faces are given
Vector3 SmallSimplexPartition3D::evaluate1FormWithBodyCoefficients(const Vector3& evaluationPoint, const VectorN& bodyCoefficients,
	const Buffer<uint>& nodes, const Buffer<uint>& edges, const Buffer<uint>& faces,
	const Vector3& gradLambda0, const Vector3& gradLambda1, const Vector3& gradLambda2, const Vector3& gradLambda3) const {

	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallSimplexPartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	Vector4 evaluationPoint4(evaluationPoint, 0.0); //for routines that require Vector4

	//find barycentric coordinates
	double lambda0, lambda1, lambda2, lambda3;
	lambda0 = mesh_old.getFaceDeviation(faces[3], evaluationPoint4).len() * gradLambda0.len();
	lambda1 = mesh_old.getFaceDeviation(faces[2], evaluationPoint4).len() * gradLambda1.len();
	lambda2 = mesh_old.getFaceDeviation(faces[1], evaluationPoint4).len() * gradLambda2.len();
	lambda3 = mesh_old.getFaceDeviation(faces[0], evaluationPoint4).len() * gradLambda3.len();

	//compute the values of the lowest order Whitney forms
	Vector3 w0 = lambda0 * gradLambda1 - lambda1 * gradLambda0;
	Vector3 w1 = lambda0 * gradLambda2 - lambda2 * gradLambda0;
	Vector3 w2 = lambda0 * gradLambda3 - lambda3 * gradLambda0;
	Vector3 w3 = lambda1 * gradLambda2 - lambda2 * gradLambda1;
	Vector3 w4 = lambda1 * gradLambda3 - lambda3 * gradLambda1;
	Vector3 w5 = lambda2 * gradLambda3 - lambda3 * gradLambda2;

	//compute the result using the coefficients given
	Buffer<double> d(6, 0.0);
	for (uint i = 0; i < activeSmallEdgesInBodies.size(); ++i) {
		const mi_t& mi = *activeSmallEdgesInBodies[i].multiIndex;
		d[getEdgeIndex(*activeSmallEdgesInBodies[i].nodeIndices)] += std::pow(lambda0, mi[0]) * std::pow(lambda1, mi[1])
			* std::pow(lambda2, mi[2]) * std::pow(lambda3, mi[3]) * bodyCoefficients[i];
	}

	return d[0] * w0 + d[1] * w1 + d[2] * w2 + d[3] * w3 + d[4] * w4 + d[5] * w5;
}

//precomputes the matrices required in evaluating higher order Whitney forms
void SmallSimplexPartition3D::formMatrices() {

	uint i, j, k;

	//0-forms
	Buffer<Buffer<double>> integrals0FormsEdges(activeSmallNodesInEdges.size());
	Buffer<Buffer<double>> integrals0FormsFaces(activeSmallNodesInFaces.size());
	Buffer<Buffer<double>> integrals0FormsBodies(activeSmallNodesInBodies.size());
	const mi_t node0EdgeMultiIndex{ {order - 1, 0} };
	const mi_t node1EdgeMultiIndex{ {0, order - 1} };
	const mi_t node0FaceMultiIndex{ {order - 1, 0, 0} };
	const mi_t node1FaceMultiIndex{ {0, order - 1, 0} };
	const mi_t node2FaceMultiIndex{ {0, 0, order - 1} };
	const mi_t node0BodyMultiIndex{ {order - 1, 0, 0, 0} };
	const mi_t node1BodyMultiIndex{ {0, order - 1, 0, 0} };
	const mi_t node2BodyMultiIndex{ {0, 0, order - 1, 0} };
	const mi_t node3BodyMultiIndex{ {0, 0, 0, order - 1} };

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

	//0-forms in bodies
	bodyValuesNode0.toVectorN(activeSmallNodesInBodies.size());
	bodyValuesNode1.toVectorN(activeSmallNodesInBodies.size());
	bodyValuesNode2.toVectorN(activeSmallNodesInBodies.size());
	bodyValuesNode3.toVectorN(activeSmallNodesInBodies.size());
	bodyNodeValuesEdge0.resize(edgePermutations.size());
	bodyNodeValuesEdge1.resize(edgePermutations.size());
	bodyNodeValuesEdge2.resize(edgePermutations.size());
	bodyNodeValuesEdge3.resize(edgePermutations.size());
	bodyNodeValuesEdge4.resize(edgePermutations.size());
	bodyNodeValuesEdge5.resize(edgePermutations.size());
	for (i = 0; i < edgePermutations.size(); ++i) {
		bodyNodeValuesEdge0[i].resize(activeSmallNodesInBodies.size());
		bodyNodeValuesEdge1[i].resize(activeSmallNodesInBodies.size());
		bodyNodeValuesEdge2[i].resize(activeSmallNodesInBodies.size());
		bodyNodeValuesEdge3[i].resize(activeSmallNodesInBodies.size());
		bodyNodeValuesEdge4[i].resize(activeSmallNodesInBodies.size());
		bodyNodeValuesEdge5[i].resize(activeSmallNodesInBodies.size());
		for (j = 0; j < activeSmallNodesInBodies.size(); ++j) {
			bodyNodeValuesEdge0[i][j].toVectorN(activeSmallNodesInEdges.size());
			bodyNodeValuesEdge1[i][j].toVectorN(activeSmallNodesInEdges.size());
			bodyNodeValuesEdge2[i][j].toVectorN(activeSmallNodesInEdges.size());
			bodyNodeValuesEdge3[i][j].toVectorN(activeSmallNodesInEdges.size());
			bodyNodeValuesEdge4[i][j].toVectorN(activeSmallNodesInEdges.size());
			bodyNodeValuesEdge5[i][j].toVectorN(activeSmallNodesInEdges.size());
		}
	}
	bodyNodeValuesFace0.resize(facePermutations.size());
	bodyNodeValuesFace1.resize(facePermutations.size());
	bodyNodeValuesFace2.resize(facePermutations.size());
	bodyNodeValuesFace3.resize(facePermutations.size());
	for (i = 0; i < facePermutations.size(); ++i) {
		bodyNodeValuesFace0[i].resize(activeSmallNodesInBodies.size());
		bodyNodeValuesFace1[i].resize(activeSmallNodesInBodies.size());
		bodyNodeValuesFace2[i].resize(activeSmallNodesInBodies.size());
		bodyNodeValuesFace3[i].resize(activeSmallNodesInBodies.size());
		for (j = 0; j < activeSmallNodesInBodies.size(); ++j) {
			bodyNodeValuesFace0[i][j].toVectorN(activeSmallNodesInFaces.size());
			bodyNodeValuesFace1[i][j].toVectorN(activeSmallNodesInFaces.size());
			bodyNodeValuesFace2[i][j].toVectorN(activeSmallNodesInFaces.size());
			bodyNodeValuesFace3[i][j].toVectorN(activeSmallNodesInFaces.size());
		}
	}
	for (i = 0; i < activeSmallNodesInBodies.size(); ++i) {
		bodyValuesNode0[i] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
			node0BodyMultiIndex, bodyNodes[0]);
		bodyValuesNode1[i] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
			node1BodyMultiIndex, bodyNodes[1]);
		bodyValuesNode2[i] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
			node2BodyMultiIndex, bodyNodes[2]);
		bodyValuesNode3[i] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
			node3BodyMultiIndex, bodyNodes[3]);
		for (j = 0; j < activeSmallNodesInEdges.size(); ++j) {
			for (k = 0; k < edgePermutations.size(); ++k) {
				bodyNodeValuesEdge0[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 0, 1), addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 0, 1));
				bodyNodeValuesEdge1[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 0, 2), addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 0, 2));
				bodyNodeValuesEdge2[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 0, 3), addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 0, 3));
				bodyNodeValuesEdge3[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 1, 2), addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 1, 2));
				bodyNodeValuesEdge4[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 1, 3), addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 1, 3));
				bodyNodeValuesEdge5[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].multiIndex, edgePermutations[k]), 2, 3), addZeros(permuteOtherWay(*activeSmallNodesInEdges[j].nodeIndices, edgePermutations[k]), 2, 3));
			}
		}
		for (j = 0; j < activeSmallNodesInFaces.size(); ++j) {
			for (k = 0; k < facePermutations.size(); ++k) {
				bodyNodeValuesFace0[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallNodesInFaces[j].multiIndex, facePermutations[k]), 3, 0), insert(permuteOtherWay(*activeSmallNodesInFaces[j].nodeIndices, facePermutations[k]), 3, 0));
				bodyNodeValuesFace1[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallNodesInFaces[j].multiIndex, facePermutations[k]), 2, 0), insert(permuteOtherWay(*activeSmallNodesInFaces[j].nodeIndices, facePermutations[k]), 2, 0));
				bodyNodeValuesFace2[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallNodesInFaces[j].multiIndex, facePermutations[k]), 1, 0), insert(permuteOtherWay(*activeSmallNodesInFaces[j].nodeIndices, facePermutations[k]), 1, 0));
				bodyNodeValuesFace3[k][i][j] = -integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallNodesInFaces[j].multiIndex, facePermutations[k]), 0, 0), insert(permuteOtherWay(*activeSmallNodesInFaces[j].nodeIndices, facePermutations[k]), 0, 0));
			}
		}
		integrals0FormsBodies[i].resize(activeSmallNodesInBodies.size());
		for (j = 0; j < activeSmallNodesInBodies.size(); ++j) {
			integrals0FormsBodies[i][j] = integrateWhitneyForm(*activeSmallNodesInBodies[i].multiIndex, *activeSmallNodesInBodies[i].nodeIndices,
				*activeSmallNodesInBodies[j].multiIndex, *activeSmallNodesInBodies[j].nodeIndices);
		}
	}

	//1-forms
	Buffer<Buffer<double>> integrals1FormsEdges(activeSmallEdgesInEdges.size());
	Buffer<Buffer<double>> integrals1FormsFaces(activeSmallEdgesInFaces.size());
	Buffer<Buffer<double>> integrals1FormsBodies(activeSmallEdgesInBodies.size());

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

	//1-forms in bodies
	bodyEdgeValuesEdge0.resize(edgePermutations.size());
	bodyEdgeValuesEdge1.resize(edgePermutations.size());
	bodyEdgeValuesEdge2.resize(edgePermutations.size());
	bodyEdgeValuesEdge3.resize(edgePermutations.size());
	bodyEdgeValuesEdge4.resize(edgePermutations.size());
	bodyEdgeValuesEdge5.resize(edgePermutations.size());
	for (i = 0; i < edgePermutations.size(); ++i) {
		bodyEdgeValuesEdge0[i].resize(activeSmallEdgesInBodies.size());
		bodyEdgeValuesEdge1[i].resize(activeSmallEdgesInBodies.size());
		bodyEdgeValuesEdge2[i].resize(activeSmallEdgesInBodies.size());
		bodyEdgeValuesEdge3[i].resize(activeSmallEdgesInBodies.size());
		bodyEdgeValuesEdge4[i].resize(activeSmallEdgesInBodies.size());
		bodyEdgeValuesEdge5[i].resize(activeSmallEdgesInBodies.size());
		for (j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
			bodyEdgeValuesEdge0[i][j].toVectorN(activeSmallEdgesInEdges.size());
			bodyEdgeValuesEdge1[i][j].toVectorN(activeSmallEdgesInEdges.size());
			bodyEdgeValuesEdge2[i][j].toVectorN(activeSmallEdgesInEdges.size());
			bodyEdgeValuesEdge3[i][j].toVectorN(activeSmallEdgesInEdges.size());
			bodyEdgeValuesEdge4[i][j].toVectorN(activeSmallEdgesInEdges.size());
			bodyEdgeValuesEdge5[i][j].toVectorN(activeSmallEdgesInEdges.size());
		}
	}
	bodyEdgeValuesFace0.resize(facePermutations.size());
	bodyEdgeValuesFace1.resize(facePermutations.size());
	bodyEdgeValuesFace2.resize(facePermutations.size());
	bodyEdgeValuesFace3.resize(facePermutations.size());
	for (i = 0; i < facePermutations.size(); ++i) {
		bodyEdgeValuesFace0[i].resize(activeSmallEdgesInBodies.size());
		bodyEdgeValuesFace1[i].resize(activeSmallEdgesInBodies.size());
		bodyEdgeValuesFace2[i].resize(activeSmallEdgesInBodies.size());
		bodyEdgeValuesFace3[i].resize(activeSmallEdgesInBodies.size());
		for (j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
			bodyEdgeValuesFace0[i][j].toVectorN(activeSmallEdgesInFaces.size());
			bodyEdgeValuesFace1[i][j].toVectorN(activeSmallEdgesInFaces.size());
			bodyEdgeValuesFace2[i][j].toVectorN(activeSmallEdgesInFaces.size());
			bodyEdgeValuesFace3[i][j].toVectorN(activeSmallEdgesInFaces.size());
		}
	}
	for (i = 0; i < activeSmallEdgesInBodies.size(); ++i) {
		for (j = 0; j < activeSmallEdgesInEdges.size(); ++j) {
			for (k = 0; k < edgePermutations.size(); ++k) {
				bodyEdgeValuesEdge0[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 0, 1), addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 0, 1));
				bodyEdgeValuesEdge1[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 0, 2), addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 0, 2));
				bodyEdgeValuesEdge2[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 0, 3), addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 0, 3));
				bodyEdgeValuesEdge3[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 1, 2), addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 1, 2));
				bodyEdgeValuesEdge4[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 1, 3), addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 1, 3));
				bodyEdgeValuesEdge5[k][i][j] = -getEdgePermutationSign(k) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].multiIndex, edgePermutations[k]), 2, 3), addZeros(permuteOtherWay(*activeSmallEdgesInEdges[j].nodeIndices, edgePermutations[k]), 2, 3));
			}
		}
		for (j = 0; j < activeSmallEdgesInFaces.size(); ++j) {
			for (k = 0; k < facePermutations.size(); ++k) {
				bodyEdgeValuesFace0[k][i][j] = -getFacePermutationSign(k, *activeSmallEdgesInFaces[j].nodeIndices) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallEdgesInFaces[j].multiIndex, facePermutations[k]), 3, 0), insert(permuteOtherWay(*activeSmallEdgesInFaces[j].nodeIndices, facePermutations[k]), 3, 0));
				bodyEdgeValuesFace1[k][i][j] = -getFacePermutationSign(k, *activeSmallEdgesInFaces[j].nodeIndices) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallEdgesInFaces[j].multiIndex, facePermutations[k]), 2, 0), insert(permuteOtherWay(*activeSmallEdgesInFaces[j].nodeIndices, facePermutations[k]), 2, 0));
				bodyEdgeValuesFace2[k][i][j] = -getFacePermutationSign(k, *activeSmallEdgesInFaces[j].nodeIndices) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallEdgesInFaces[j].multiIndex, facePermutations[k]), 1, 0), insert(permuteOtherWay(*activeSmallEdgesInFaces[j].nodeIndices, facePermutations[k]), 1, 0));
				bodyEdgeValuesFace3[k][i][j] = -getFacePermutationSign(k, *activeSmallEdgesInFaces[j].nodeIndices) * integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallEdgesInFaces[j].multiIndex, facePermutations[k]), 0, 0), insert(permuteOtherWay(*activeSmallEdgesInFaces[j].nodeIndices, facePermutations[k]), 0, 0));
			}
		}
		integrals1FormsBodies[i].resize(activeSmallEdgesInBodies.size());
		for (j = 0; j < activeSmallEdgesInBodies.size(); ++j) {
			integrals1FormsBodies[i][j] = integrateWhitneyForm(*activeSmallEdgesInBodies[i].multiIndex, *activeSmallEdgesInBodies[i].nodeIndices,
				*activeSmallEdgesInBodies[j].multiIndex, *activeSmallEdgesInBodies[j].nodeIndices);
		}
	}

	//2-forms
	Buffer<Buffer<double>> integrals2FormsFaces(activeSmallFacesInFaces.size());
	Buffer<Buffer<double>> integrals2FormsBodies(activeSmallFacesInBodies.size());

	//2-forms in faces
	for (i = 0; i < activeSmallFacesInFaces.size(); ++i) {
		integrals2FormsFaces[i].resize(activeSmallFacesInFaces.size());
		for (j = 0; j < activeSmallFacesInFaces.size(); ++j) {
			integrals2FormsFaces[i][j] = integrateWhitneyForm(*activeSmallFacesInFaces[i].multiIndex, *activeSmallFacesInFaces[i].nodeIndices,
				*activeSmallFacesInFaces[j].multiIndex, *activeSmallFacesInFaces[j].nodeIndices);
		}
	}

	//2-forms in bodies
	bodyFaceValuesFace0.resize(facePermutations.size());
	bodyFaceValuesFace1.resize(facePermutations.size());
	bodyFaceValuesFace2.resize(facePermutations.size());
	bodyFaceValuesFace3.resize(facePermutations.size());
	for (i = 0; i < facePermutations.size(); ++i) {
		bodyFaceValuesFace0[i].resize(activeSmallFacesInBodies.size());
		bodyFaceValuesFace1[i].resize(activeSmallFacesInBodies.size());
		bodyFaceValuesFace2[i].resize(activeSmallFacesInBodies.size());
		bodyFaceValuesFace3[i].resize(activeSmallFacesInBodies.size());
		for (j = 0; j < activeSmallFacesInBodies.size(); ++j) {
			bodyFaceValuesFace0[i][j].toVectorN(activeSmallFacesInFaces.size());
			bodyFaceValuesFace1[i][j].toVectorN(activeSmallFacesInFaces.size());
			bodyFaceValuesFace2[i][j].toVectorN(activeSmallFacesInFaces.size());
			bodyFaceValuesFace3[i][j].toVectorN(activeSmallFacesInFaces.size());
		}
	}
	for (i = 0; i < activeSmallFacesInBodies.size(); ++i) {
		for (j = 0; j < activeSmallFacesInFaces.size(); ++j) {
			for (k = 0; k < facePermutations.size(); ++k) {
				bodyFaceValuesFace0[k][i][j] = -getFacePermutationSign(k) * integrateWhitneyForm(*activeSmallFacesInBodies[i].multiIndex, *activeSmallFacesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallFacesInFaces[j].multiIndex, facePermutations[k]), 3, 0), insert(permuteOtherWay(*activeSmallFacesInFaces[j].nodeIndices, facePermutations[k]), 3, 0));
				bodyFaceValuesFace1[k][i][j] = -getFacePermutationSign(k) * integrateWhitneyForm(*activeSmallFacesInBodies[i].multiIndex, *activeSmallFacesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallFacesInFaces[j].multiIndex, facePermutations[k]), 2, 0), insert(permuteOtherWay(*activeSmallFacesInFaces[j].nodeIndices, facePermutations[k]), 2, 0));
				bodyFaceValuesFace2[k][i][j] = -getFacePermutationSign(k) * integrateWhitneyForm(*activeSmallFacesInBodies[i].multiIndex, *activeSmallFacesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallFacesInFaces[j].multiIndex, facePermutations[k]), 1, 0), insert(permuteOtherWay(*activeSmallFacesInFaces[j].nodeIndices, facePermutations[k]), 1, 0));
				bodyFaceValuesFace3[k][i][j] = -getFacePermutationSign(k) * integrateWhitneyForm(*activeSmallFacesInBodies[i].multiIndex, *activeSmallFacesInBodies[i].nodeIndices,
					insert(permuteOtherWay(*activeSmallFacesInFaces[j].multiIndex, facePermutations[k]), 0, 0), insert(permuteOtherWay(*activeSmallFacesInFaces[j].nodeIndices, facePermutations[k]), 0, 0));
			}
		}
		integrals2FormsBodies[i].resize(activeSmallFacesInBodies.size());
		for (j = 0; j < activeSmallFacesInBodies.size(); ++j) {
			integrals2FormsBodies[i][j] = integrateWhitneyForm(*activeSmallFacesInBodies[i].multiIndex, *activeSmallFacesInBodies[i].nodeIndices,
				*activeSmallFacesInBodies[j].multiIndex, *activeSmallFacesInBodies[j].nodeIndices);
		}
	}

	//3-forms
	Buffer<Buffer<double>> integrals3Forms(activeSmallBodiesInBodies.size());
	for (i = 0; i < activeSmallBodiesInBodies.size(); ++i) {
		integrals3Forms[i].resize(activeSmallBodiesInBodies.size());
		for (j = 0; j < activeSmallBodiesInBodies.size(); ++j) {
			integrals3Forms[i][j] = integrateWhitneyForm(*activeSmallBodiesInBodies[i].multiIndex, *activeSmallBodiesInBodies[i].nodeIndices,
				*activeSmallBodiesInBodies[j].multiIndex, *activeSmallBodiesInBodies[j].nodeIndices);
		}
	}

	//precompute LU decompositions
	double tol = 1e-15;
	decomposeLUP(integrals0FormsEdges, matrix0FormsEdges_p, tol);
	decomposeLUP(integrals0FormsFaces, matrix0FormsFaces_p, tol);
	decomposeLUP(integrals0FormsBodies, matrix0FormsBodies_p, tol);
	decomposeLUP(integrals1FormsEdges, matrix1FormsEdges_p, tol);
	decomposeLUP(integrals1FormsFaces, matrix1FormsFaces_p, tol);
	decomposeLUP(integrals1FormsBodies, matrix1FormsBodies_p, tol);
	decomposeLUP(integrals2FormsFaces, matrix2FormsFaces_p, tol);
	decomposeLUP(integrals2FormsBodies, matrix2FormsBodies_p, tol);
	decomposeLUP(integrals3Forms, matrix3Forms_p, tol);
	//convert to MatrixN
	convertMatrix(integrals0FormsEdges, matrix0FormsEdges);
	convertMatrix(integrals0FormsFaces, matrix0FormsFaces);
	convertMatrix(integrals0FormsBodies, matrix0FormsBodies);
	convertMatrix(integrals1FormsEdges, matrix1FormsEdges);
	convertMatrix(integrals1FormsFaces, matrix1FormsFaces);
	convertMatrix(integrals1FormsBodies, matrix1FormsBodies);
	convertMatrix(integrals2FormsFaces, matrix2FormsFaces);
	convertMatrix(integrals2FormsBodies, matrix2FormsBodies);
	convertMatrix(integrals3Forms, matrix3Forms);
}

//save the matrices that have been precomputed
void SmallSimplexPartition3D::saveMatrices() const {

	if (order <= 3) //not all the matrices exist, and those who do are easily formed
		return;

	Text path;

	path << "Files/3D/order" << order << "/matrix0FormsEdges.dat";
	saveMatrix(matrix0FormsEdges, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsFaces.dat";
	saveMatrix(matrix0FormsFaces, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsBodies.dat";
	saveMatrix(matrix0FormsBodies, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsEdges.dat";
	saveMatrix(matrix1FormsEdges, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsFaces.dat";
	saveMatrix(matrix1FormsFaces, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsBodies.dat";
	saveMatrix(matrix1FormsBodies, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix2FormsFaces.dat";
	saveMatrix(matrix2FormsFaces, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix2FormsBodies.dat";
	saveMatrix(matrix2FormsBodies, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix3Forms.dat";
	saveMatrix(matrix3Forms, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsEdges_p.dat";
	saveBuffer(matrix0FormsEdges_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsFaces_p.dat";
	saveBuffer(matrix0FormsFaces_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsBodies_p.dat";
	saveBuffer(matrix0FormsBodies_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsEdges_p.dat";
	saveBuffer(matrix1FormsEdges_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsFaces_p.dat";
	saveBuffer(matrix1FormsFaces_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsBodies_p.dat";
	saveBuffer(matrix1FormsBodies_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix2FormsFaces_p.dat";
	saveBuffer(matrix2FormsFaces_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix2FormsBodies_p.dat";
	saveBuffer(matrix2FormsBodies_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix3Forms_p.dat";
	saveBuffer(matrix3Forms_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/edgeValuesNode0.dat";
	saveBuffer(edgeValuesNode0.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/edgeValuesNode1.dat";
	saveBuffer(edgeValuesNode1.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceValuesNode0.dat";
	saveBuffer(faceValuesNode0.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceValuesNode1.dat";
	saveBuffer(faceValuesNode1.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceValuesNode2.dat";
	saveBuffer(faceValuesNode2.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyValuesNode0.dat";
	saveBuffer(bodyValuesNode0.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyValuesNode1.dat";
	saveBuffer(bodyValuesNode1.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyValuesNode2.dat";
	saveBuffer(bodyValuesNode2.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyValuesNode3.dat";
	saveBuffer(bodyValuesNode3.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceNodeValuesEdge0.dat";
	saveMatrix(faceNodeValuesEdge0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceNodeValuesEdge1.dat";
	saveMatrix(faceNodeValuesEdge1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceNodeValuesEdge2.dat";
	saveMatrix(faceNodeValuesEdge2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge0.dat";
	saveMatrix(bodyNodeValuesEdge0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge1.dat";
	saveMatrix(bodyNodeValuesEdge1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge2.dat";
	saveMatrix(bodyNodeValuesEdge2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge3.dat";
	saveMatrix(bodyNodeValuesEdge3, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge4.dat";
	saveMatrix(bodyNodeValuesEdge4, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge5.dat";
	saveMatrix(bodyNodeValuesEdge5, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceEdgeValuesEdge0.dat";
	saveMatrix(faceEdgeValuesEdge0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceEdgeValuesEdge1.dat";
	saveMatrix(faceEdgeValuesEdge1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceEdgeValuesEdge2.dat";
	saveMatrix(faceEdgeValuesEdge2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge0.dat";
	saveMatrix(bodyEdgeValuesEdge0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge1.dat";
	saveMatrix(bodyEdgeValuesEdge1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge2.dat";
	saveMatrix(bodyEdgeValuesEdge2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge3.dat";
	saveMatrix(bodyEdgeValuesEdge3, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge4.dat";
	saveMatrix(bodyEdgeValuesEdge4, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge5.dat";
	saveMatrix(bodyEdgeValuesEdge5, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesFace0.dat";
	saveMatrix(bodyNodeValuesFace0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesFace1.dat";
	saveMatrix(bodyNodeValuesFace1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesFace2.dat";
	saveMatrix(bodyNodeValuesFace2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesFace3.dat";
	saveMatrix(bodyNodeValuesFace3, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesFace0.dat";
	saveMatrix(bodyEdgeValuesFace0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesFace1.dat";
	saveMatrix(bodyEdgeValuesFace1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesFace2.dat";
	saveMatrix(bodyEdgeValuesFace2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesFace3.dat";
	saveMatrix(bodyEdgeValuesFace3, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyFaceValuesFace0.dat";
	saveMatrix(bodyFaceValuesFace0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyFaceValuesFace1.dat";
	saveMatrix(bodyFaceValuesFace1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyFaceValuesFace2.dat";
	saveMatrix(bodyFaceValuesFace2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyFaceValuesFace3.dat";
	saveMatrix(bodyFaceValuesFace3, path.str());
	path.clear();
}

//load the matrices that have been precomputed
void SmallSimplexPartition3D::loadMatrices() {

	if (order <= 3) { //not all the matrices exist, and those who do are easily formed
		formMatrices();
		return;
	}

	Text path;

	path << "Files/3D/order" << order << "/matrix0FormsEdges.dat";
	loadMatrix(matrix0FormsEdges, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsFaces.dat";
	loadMatrix(matrix0FormsFaces, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsBodies.dat";
	loadMatrix(matrix0FormsBodies, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsEdges.dat";
	loadMatrix(matrix1FormsEdges, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsFaces.dat";
	loadMatrix(matrix1FormsFaces, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsBodies.dat";
	loadMatrix(matrix1FormsBodies, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix2FormsFaces.dat";
	loadMatrix(matrix2FormsFaces, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix2FormsBodies.dat";
	loadMatrix(matrix2FormsBodies, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix3Forms.dat";
	loadMatrix(matrix3Forms, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsEdges_p.dat";
	loadBuffer(matrix0FormsEdges_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsFaces_p.dat";
	loadBuffer(matrix0FormsFaces_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix0FormsBodies_p.dat";
	loadBuffer(matrix0FormsBodies_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsEdges_p.dat";
	loadBuffer(matrix1FormsEdges_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsFaces_p.dat";
	loadBuffer(matrix1FormsFaces_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix1FormsBodies_p.dat";
	loadBuffer(matrix1FormsBodies_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix2FormsFaces_p.dat";
	loadBuffer(matrix2FormsFaces_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix2FormsBodies_p.dat";
	loadBuffer(matrix2FormsBodies_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/matrix3Forms_p.dat";
	loadBuffer(matrix3Forms_p, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/edgeValuesNode0.dat";
	loadBuffer(edgeValuesNode0.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/edgeValuesNode1.dat";
	loadBuffer(edgeValuesNode1.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceValuesNode0.dat";
	loadBuffer(faceValuesNode0.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceValuesNode1.dat";
	loadBuffer(faceValuesNode1.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceValuesNode2.dat";
	loadBuffer(faceValuesNode2.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyValuesNode0.dat";
	loadBuffer(bodyValuesNode0.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyValuesNode1.dat";
	loadBuffer(bodyValuesNode1.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyValuesNode2.dat";
	loadBuffer(bodyValuesNode2.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyValuesNode3.dat";
	loadBuffer(bodyValuesNode3.val, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceNodeValuesEdge0.dat";
	loadMatrix(faceNodeValuesEdge0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceNodeValuesEdge1.dat";
	loadMatrix(faceNodeValuesEdge1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceNodeValuesEdge2.dat";
	loadMatrix(faceNodeValuesEdge2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge0.dat";
	loadMatrix(bodyNodeValuesEdge0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge1.dat";
	loadMatrix(bodyNodeValuesEdge1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge2.dat";
	loadMatrix(bodyNodeValuesEdge2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge3.dat";
	loadMatrix(bodyNodeValuesEdge3, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge4.dat";
	loadMatrix(bodyNodeValuesEdge4, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesEdge5.dat";
	loadMatrix(bodyNodeValuesEdge5, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceEdgeValuesEdge0.dat";
	loadMatrix(faceEdgeValuesEdge0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceEdgeValuesEdge1.dat";
	loadMatrix(faceEdgeValuesEdge1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/faceEdgeValuesEdge2.dat";
	loadMatrix(faceEdgeValuesEdge2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge0.dat";
	loadMatrix(bodyEdgeValuesEdge0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge1.dat";
	loadMatrix(bodyEdgeValuesEdge1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge2.dat";
	loadMatrix(bodyEdgeValuesEdge2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge3.dat";
	loadMatrix(bodyEdgeValuesEdge3, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge4.dat";
	loadMatrix(bodyEdgeValuesEdge4, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesEdge5.dat";
	loadMatrix(bodyEdgeValuesEdge5, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesFace0.dat";
	loadMatrix(bodyNodeValuesFace0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesFace1.dat";
	loadMatrix(bodyNodeValuesFace1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesFace2.dat";
	loadMatrix(bodyNodeValuesFace2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyNodeValuesFace3.dat";
	loadMatrix(bodyNodeValuesFace3, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesFace0.dat";
	loadMatrix(bodyEdgeValuesFace0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesFace1.dat";
	loadMatrix(bodyEdgeValuesFace1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesFace2.dat";
	loadMatrix(bodyEdgeValuesFace2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyEdgeValuesFace3.dat";
	loadMatrix(bodyEdgeValuesFace3, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyFaceValuesFace0.dat";
	loadMatrix(bodyFaceValuesFace0, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyFaceValuesFace1.dat";
	loadMatrix(bodyFaceValuesFace1, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyFaceValuesFace2.dat";
	loadMatrix(bodyFaceValuesFace2, path.str());
	path.clear();
	path << "Files/3D/order" << order << "/bodyFaceValuesFace3.dat";
	loadMatrix(bodyFaceValuesFace3, path.str());
	path.clear();
}

/*
Save the integrals required in building the Hodge matrix for 1-forms for a specific element type (which should be indicated by foldername).
The length of the first edge of the element in which the computations were performed must be given as scalingFactor.
When this is known, the same integrals can be used in congruent elements of different sizes after scaling.
*/
void SmallSimplexPartition3D::saveHodgeIntegrals1Forms(std::string foldername, bool circumcentric, double scalingFactor, const Buffer<VectorN>& hodgeIntegralsBodyEdges,
	const Buffer<Buffer<VectorN>>& hodgeIntegralsFace0Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsFace1Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsFace2Edges,
	const Buffer<Buffer<VectorN>>& hodgeIntegralsFace3Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge0Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge1Edges,
	const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge2Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge3Edges, const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge4Edges,
	const Buffer<Buffer<VectorN>>& hodgeIntegralsEdge5Edges, int ref) const {
	std::string dual = circumcentric ? "circumcentric" : "barycentric";
	std::string refElement = (ref < 0 ? "" : "Type" + std::to_string(ref));
	Text path;
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/scalingFactor" << refElement.c_str() << ".dat";
	saveElement(scalingFactor, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsBodyEdges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsBodyEdges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsFace0Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsFace0Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsFace1Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsFace1Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsFace2Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsFace2Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsFace3Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsFace3Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge0Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge0Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge1Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge1Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge2Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge2Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge3Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge3Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge4Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge4Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername.c_str() << "/hodgeIntegralsEdge5Edges" << refElement.c_str() << ".dat";
	saveMatrix(hodgeIntegralsEdge5Edges, path.str());
	path.clear();
}

/*
Load the integrals required in building the Hodge matrix for 1-forms for a specific element type (and have been stored in a folder indicated by foldername).
Returns the length of the first edge of the element in which the integrals were computed.
When this is known, the same integrals can be used in congruent elements of different sizes after scaling.
*/
double SmallSimplexPartition3D::loadHodgeIntegrals1Forms(std::string foldername, bool circumcentric, Buffer<VectorN>& hodgeIntegralsBodyEdges, Buffer<Buffer<VectorN>>& hodgeIntegralsFace0Edges,
	Buffer<Buffer<VectorN>>& hodgeIntegralsFace1Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsFace2Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsFace3Edges,
	Buffer<Buffer<VectorN>>& hodgeIntegralsEdge0Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge1Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge2Edges,
	Buffer<Buffer<VectorN>>& hodgeIntegralsEdge3Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge4Edges, Buffer<Buffer<VectorN>>& hodgeIntegralsEdge5Edges, int ref) const {
	std::string dual = circumcentric ? "circumcentric" : "barycentric";
	std::string refElement = (ref < 0 ? "" : "Type" + std::to_string(ref));
	Text path;
	double scalingFactor;
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/scalingFactor" << refElement.c_str() << ".dat";
	loadElement(scalingFactor, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsBodyEdges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsBodyEdges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsFace0Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsFace0Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsFace1Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsFace1Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsFace2Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsFace2Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsFace3Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsFace3Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsEdge0Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge0Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsEdge1Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge1Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsEdge2Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge2Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsEdge3Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge3Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsEdge4Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge4Edges, path.str());
	path.clear();
	path << "Files/3D/order" << order << '/' << dual.c_str() << '/' << foldername << "/hodgeIntegralsEdge5Edges" << refElement.c_str() << ".dat";
	loadMatrix(hodgeIntegralsEdge5Edges, path.str());
	path.clear();
	return scalingFactor;
}
