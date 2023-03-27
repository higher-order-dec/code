#include "MeshHelperFunctions.hpp"
#include "NumericalIntegration.hpp"
#include <iostream>

namespace gfd {

	//find the edges of the triangle with the given nodes (assuming one exists)
	Buffer<uint> findEdges(const uint n0, const uint n1, const uint n2, const Mesh& mesh) {
		Buffer<uint> e(3);
		e[0] = mesh.findEdge(n0, n1);
		e[1] = mesh.findEdge(n1, n2);
		e[2] = mesh.findEdge(n2, n0);
		return e;
	}

	//find the faces of the tetrahedron with the given nodes (assuming one exists)
	Buffer<uint> findFaces(const uint n0, const uint n1, const uint n2, const uint n3, const Mesh& mesh) {
		Buffer<uint> f(4);
		f[0] = mesh.findFace(findEdges(n0, n1, n2, mesh));
		f[1] = mesh.findFace(findEdges(n0, n1, n3, mesh));
		f[2] = mesh.findFace(findEdges(n0, n2, n3, mesh));
		f[3] = mesh.findFace(findEdges(n1, n2, n3, mesh));
		return f;
	}

	//find the simplex with given nodes in simplicial mesh
	uint findSimplexWithNodes(const Buffer<uint>& simplexNodes, const Mesh& mesh) {
		switch (simplexNodes.size()) {
		case 1:
			return simplexNodes[0];
		case 2:
			return mesh.findEdge(simplexNodes[0], simplexNodes[1]);
		case 3:
			return mesh.findFace(findEdges(simplexNodes[0], simplexNodes[1], simplexNodes[2], mesh));
		case 4:
			return mesh.findBody(findFaces(simplexNodes[0], simplexNodes[1], simplexNodes[2], simplexNodes[3], mesh));
		default:
			return NONE;
		}
	}

	//find the quadrilateral with the given nodes
	uint findQuadrilateral(const uint n0, const uint n1, const uint n2, const uint n3, const Mesh& mesh) {
		Buffer<uint> e(4);
		uint i = 0;
		uint edgeIndex = mesh.findEdge(n0, n1);
		if (edgeIndex != NONE)
			e[i++] = edgeIndex;
		edgeIndex = mesh.findEdge(n0, n2);
		if (edgeIndex != NONE)
			e[i++] = edgeIndex;
		edgeIndex = mesh.findEdge(n0, n3);
		if (edgeIndex != NONE)
			e[i++] = edgeIndex;
		edgeIndex = mesh.findEdge(n1, n2);
		if (edgeIndex != NONE)
			e[i++] = edgeIndex;
		edgeIndex = mesh.findEdge(n1, n3);
		if (edgeIndex != NONE)
			e[i++] = edgeIndex;
		edgeIndex = mesh.findEdge(n2, n3);
		if (edgeIndex != NONE)
			e[i++] = edgeIndex;
		return mesh.findFace(e);
	}

	//return the nodes of face i in Cartesian 3D mesh in the order that defines the orientation
	Buffer<uint> getQuadrilateralNodes(uint i, const Mesh& mesh) {
		if (mesh.getDimension() == 2)
			return mesh.getFaceNodes(i); //in two dimensions all calls to this can be replaced with mesh.getFaceNodes(i)

		Buffer<uint> result(4, NONE);
		const Buffer<uint> n = mesh.getFaceNodes(i);
		Buffer<Vector3> pos(4);
		pos[0] = mesh.getNodePosition3(n[0]);
		pos[1] = mesh.getNodePosition3(n[1]);
		pos[2] = mesh.getNodePosition3(n[2]);
		pos[3] = mesh.getNodePosition3(n[3]);
		Vector3 faceDualVec = mesh.getFaceVector3(i).dual();
		if (isXDir(faceDualVec)) {
			double miny = pos[0].y;
			double minz = pos[0].z;
			for (uint j = 1; j < pos.size(); ++j) {
				if (pos[j].y < miny)
					miny = pos[j].y;
				if (pos[j].z < minz)
					minz = pos[j].z;
			}
			for (uint j = 0; j < pos.size(); ++j) {
				if (pos[j].z < minz + 1e-12) {
					if (pos[j].y < miny + 1e-12)
						result[0] = n[j];
					else
						result[1] = n[j];
				}
				else {
					if (pos[j].y < miny + 1e-12)
						result[3] = n[j];
					else
						result[2] = n[j];
				}
			}
		}
		else if (isYDir(faceDualVec)) {
			double minx = pos[0].x;
			double minz = pos[0].z;
			for (uint j = 1; j < pos.size(); ++j) {
				if (pos[j].x < minx)
					minx = pos[j].x;
				if (pos[j].z < minz)
					minz = pos[j].z;
			}
			for (uint j = 0; j < pos.size(); ++j) {
				if (pos[j].z < minz + 1e-12) {
					if (pos[j].x < minx + 1e-12)
						result[0] = n[j];
					else
						result[1] = n[j];
				}
				else {
					if (pos[j].x < minx + 1e-12)
						result[3] = n[j];
					else
						result[2] = n[j];
				}
			}
		}
		else {
			double minx = pos[0].x;
			double miny = pos[0].y;
			for (uint j = 1; j < pos.size(); ++j) {
				if (pos[j].x < minx)
					minx = pos[j].x;
				if (pos[j].y < miny)
					miny = pos[j].y;
			}
			for (uint j = 0; j < pos.size(); ++j) {
				if (pos[j].y < miny + 1e-12) {
					if (pos[j].x < minx + 1e-12)
						result[0] = n[j];
					else
						result[1] = n[j];
				}
				else {
					if (pos[j].x < minx + 1e-12)
						result[3] = n[j];
					else
						result[2] = n[j];
				}
			}
		}
		return result;
	}

	//return the nodes of body i in Cartesian 3D mesh in the order that Mathematica accepts
	Buffer<uint> getCubeNodes(uint i, const Mesh& mesh) {
		const Buffer<uint>& faces = mesh.getBodyFaces(i);
		return getQuadrilateralNodes(faces[0], mesh).getUnion(getQuadrilateralNodes(faces[5], mesh));
	}

	//create a mesh consisting of only one reference simplex
	void createOneElementMesh(Mesh& mesh) {
		if (mesh.getDimension() == 2) {
			uint node0 = mesh.addNode({ 0,0,0,0 });
			uint node1 = mesh.addNode({ 1,0,0,0 });
			uint node2 = mesh.addNode({ 0,1,0,0 });
			Buffer<uint> edges(3);
			edges[0] = mesh.addEdge(node0, node1);
			edges[1] = mesh.addEdge(node0, node2);
			edges[2] = mesh.addEdge(node1, node2);
			mesh.addFace(edges);
		}
		else if (mesh.getDimension() == 3) {
			uint node0 = mesh.addNode({ 0,0,0,0 });
			uint node1 = mesh.addNode({ 1,0,0,0 });
			uint node2 = mesh.addNode({ 0,1,0,0 });
			uint node3 = mesh.addNode({ 0,0,1,0 });
			uint edge01 = mesh.addEdge(node0, node1);
			uint edge02 = mesh.addEdge(node0, node2);
			uint edge03 = mesh.addEdge(node0, node3);
			uint edge12 = mesh.addEdge(node1, node2);
			uint edge13 = mesh.addEdge(node1, node3);
			uint edge23 = mesh.addEdge(node2, node3);
			Buffer<uint> edges(3);
			Buffer<uint> faces(4);
			edges[0] = edge01;
			edges[1] = edge12;
			edges[2] = edge02;
			faces[0] = mesh.addFace(edges);
			edges[1] = edge13;
			edges[2] = edge03;
			faces[1] = mesh.addFace(edges);
			edges[0] = edge02;
			edges[1] = edge23;
			faces[2] = mesh.addFace(edges);
			edges[0] = edge12;
			edges[2] = edge13;
			faces[3] = mesh.addFace(edges);
			mesh.addBody(faces);
		}
	}

	//create a 3D mesh consisting of two simplices
	void createTwoElementMesh(Mesh& mesh) {
		if (mesh.getDimension() == 3) {
			uint node0 = mesh.addNode({ 0,0,0,0 });
			uint node1 = mesh.addNode({ 1,0,0,0 });
			uint node2 = mesh.addNode({ 0,1,0,0 });
			uint node3 = mesh.addNode({ 0.5,0.5,1,0 });
			uint node4 = mesh.addNode({ 1,1,0,0 });
			uint edge01 = mesh.addEdge(node0, node1);
			uint edge02 = mesh.addEdge(node0, node2);
			uint edge03 = mesh.addEdge(node0, node3);
			uint edge12 = mesh.addEdge(node1, node2);
			uint edge13 = mesh.addEdge(node1, node3);
			uint edge23 = mesh.addEdge(node2, node3);
			Buffer<uint> edges(3);
			Buffer<uint> faces(4);
			edges[0] = edge01;
			edges[1] = edge12;
			edges[2] = edge02;
			faces[0] = mesh.addFace(edges);
			edges[1] = edge13;
			edges[2] = edge03;
			faces[1] = mesh.addFace(edges);
			edges[0] = edge02;
			edges[1] = edge23;
			faces[2] = mesh.addFace(edges);
			edges[0] = edge12;
			edges[2] = edge13;
			faces[3] = mesh.addFace(edges);
			mesh.addBody(faces);
			uint edge41 = mesh.addEdge(node4, node1);
			uint edge42 = mesh.addEdge(node4, node2);
			uint edge43 = mesh.addEdge(node4, node3);
			edges[0] = edge41;
			edges[1] = edge12;
			edges[2] = edge42;
			faces[0] = mesh.addFace(edges);
			edges[1] = edge13;
			edges[2] = edge43;
			faces[1] = mesh.addFace(edges);
			edges[0] = edge42;
			edges[1] = edge23;
			faces[2] = mesh.addFace(edges);
			edges[0] = edge12;
			edges[2] = edge13;
			faces[3] = mesh.addFace(edges);
			mesh.addBody(faces);
		}
	}

	//form the lowest order Hodge matrix for 1-forms explicitly
	void formDiagonalHodge1Forms(Buffer<Buffer<double>>& star, const Mesh& mesh) {
		star.resize(mesh.getEdgeSize());
		for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
			star[i] = Buffer<double>(mesh.getEdgeSize(), 0.0);
			star[i][i] = mesh.getEdgeHodge(i);
		}
	}

	//Check whether the point p is inside the triangle with vertices p0, p1, and p2. The outward normal vectors of the opposite edges are required as input.
	bool isInsideTriangle(const Vector2& p, const Vector2& p0, const Vector2& p1, const Vector2& p2,
		const Vector2& p0oppositeNormal, const Vector2& p1oppositeNormal, const Vector2& p2oppositeNormal) {
		return (p0oppositeNormal.dot(p1 - p) >= 0 && p1oppositeNormal.dot(p0 - p) >= 0 && p2oppositeNormal.dot(p0 - p) >= 0);
	}

	//Check whether the point p is inside the tetrahedron with vertices p0, p1, p2, and p3. The outward normal vectors of the opposite faces are required as input.
	bool isInsideTetrahedron(const Vector3& p, const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3,
		const Vector3& p0oppositeNormal, const Vector3& p1oppositeNormal, const Vector3& p2oppositeNormal, const Vector3& p3oppositeNormal) {
		return (p0oppositeNormal.dot(p1 - p) >= 0 && p1oppositeNormal.dot(p0 - p) >= 0 && p2oppositeNormal.dot(p0 - p) >= 0 && p3oppositeNormal.dot(p0 - p) >= 0);
	}

	//Print the nodes of each element in simplicial or cubical mesh in 2D or 3D.
	void printElementNodesForMathematica(const Mesh& mesh) {
		if (mesh.getDimension() == 2) {
			for (uint f = 0; f < mesh.getFaceSize(); ++f) {
				Buffer<uint> nodes = mesh.getFaceNodes(f);
				if (nodes.size() == 3) {
					Vector2 p0 = mesh.getNodePosition2(nodes[0]);
					Vector2 p1 = mesh.getNodePosition2(nodes[1]);
					Vector2 p2 = mesh.getNodePosition2(nodes[2]);
					std::cout << "{{" << p0.x << ", " << p0.y << "}, {";
					std::cout << p1.x << ", " << p1.y << "}, {";
					std::cout << p2.x << ", " << p2.y << "}},\n";
				}
				if (nodes.size() == 4) {
					nodes = getQuadrilateralNodes(f, mesh);
					Vector2 p0 = mesh.getNodePosition2(nodes[0]);
					Vector2 p1 = mesh.getNodePosition2(nodes[1]);
					Vector2 p2 = mesh.getNodePosition2(nodes[2]);
					Vector2 p3 = mesh.getNodePosition2(nodes[3]);
					std::cout << "{{" << p0.x << ", " << p0.y << "}, {";
					std::cout << p1.x << ", " << p1.y << "}, {";
					std::cout << p2.x << ", " << p2.y << "}, {";
					std::cout << p3.x << ", " << p3.y << "}},\n";
				}
			}
		}
		else {
			for (uint b = 0; b < mesh.getBodySize(); ++b) {
				Buffer<uint> nodes = mesh.getBodyNodes(b);
				if (nodes.size() == 4) {
					Vector3 p0 = mesh.getNodePosition3(nodes[0]);
					Vector3 p1 = mesh.getNodePosition3(nodes[1]);
					Vector3 p2 = mesh.getNodePosition3(nodes[2]);
					Vector3 p3 = mesh.getNodePosition3(nodes[3]);
					std::cout << "{{" << p0.x << ", " << p0.y << ", " << p0.z << "}, {";
					std::cout << p1.x << ", " << p1.y << ", " << p1.z << "}, {";
					std::cout << p2.x << ", " << p2.y << ", " << p2.z << "}, {";
					std::cout << p3.x << ", " << p3.y << ", " << p3.z << "}},\n";
				}
				else if (nodes.size() == 8) {
					nodes = getCubeNodes(b, mesh);
					Vector3 p0 = mesh.getNodePosition3(nodes[0]);
					Vector3 p1 = mesh.getNodePosition3(nodes[1]);
					Vector3 p2 = mesh.getNodePosition3(nodes[2]);
					Vector3 p3 = mesh.getNodePosition3(nodes[3]);
					Vector3 p4 = mesh.getNodePosition3(nodes[4]);
					Vector3 p5 = mesh.getNodePosition3(nodes[5]);
					Vector3 p6 = mesh.getNodePosition3(nodes[6]);
					Vector3 p7 = mesh.getNodePosition3(nodes[7]);
					std::cout << "{{" << p0.x << ", " << p0.y << ", " << p0.z << "}, {";
					std::cout << p1.x << ", " << p1.y << ", " << p1.z << "}, {";
					std::cout << p2.x << ", " << p2.y << ", " << p2.z << "}, {";
					std::cout << p3.x << ", " << p3.y << ", " << p3.z << "}, {";
					std::cout << p4.x << ", " << p4.y << ", " << p4.z << "}, {";
					std::cout << p5.x << ", " << p5.y << ", " << p5.z << "}, {";
					std::cout << p6.x << ", " << p6.y << ", " << p6.z << "}, {";
					std::cout << p7.x << ", " << p7.y << ", " << p7.z << "}},\n";
				}
			}
		}
	}

	//computes and prints the minimum, maximum, and avarege edge length in mesh
	void printEdgeLengthData(const Mesh& mesh) {
		double sum = 0.0;
		double min = 1.0e+300;
		double max = 0.0;
		for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
			double len = mesh.getEdgeVector(i).len();
			sum += len;
			if (len < min)
				min = len;
			if (len > max)
				max = len;
		}
		std::cout << "edge length minimum " << min << " maximum " << max << " average " << sum / mesh.getEdgeSize() << '\n';
	}

	//print node positions
	void printNodePositions(const Mesh& mesh) {
		if (mesh.getDimension() == 2) {
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				Vector2 pos = mesh.getNodePosition2(i);
				std::cout << '(' << pos.x << ", " << pos.y << ")\n";
			}
		}
		else {
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				Vector3 pos = mesh.getNodePosition3(i);
				std::cout << '(' << pos.x << ", " << pos.y << ", " << pos.z << ")\n";
			}
		}
	}

	//computes the length of the dual edge of edge i (assuming 2D mesh)
	double dualEdgeLength(uint i, const Mesh& mesh, bool circumcentric) {
		if (mesh.getDimension() == 2) {
			Vector2 edgePosition = circumcentric ? mesh.getEdgePosition2(i) : mesh.getEdgeAverage2(i);
			const Buffer<uint>& faces = mesh.getEdgeFaces(i);
			double sum = 0;
			for (uint i = 0; i < faces.size(); ++i) {
				Vector2 facePosition = circumcentric ? mesh.getFacePosition2(faces[i]) : mesh.getFaceAverage2(faces[i]);
				sum += (facePosition - edgePosition).len();
			}
			return sum;
		}
		else {
			return 0.0;
		}
	}

	//computes the area of the dual face of node i (if 2D mesh) or edge i (if 3D mesh)
	double dualFaceArea(uint i, const Mesh& mesh, bool circumcentric) {
		if (mesh.getDimension() == 2) {
			Vector2 p0 = mesh.getNodePosition2(i);
			Buffer<uint> faces = mesh.getNodeFaces(i);
			const Buffer<uint>& edges = mesh.getNodeEdges(i);
			double result = 0.0;
			for (uint f = 0; f < faces.size(); ++f) {
				Buffer<uint> faceEdges = edges.getIntersection(mesh.getFaceEdges(faces[f]));
				Vector2 face_bc, edge0_bc, edge1_bc;
				if (circumcentric) {
					face_bc = mesh.getFacePosition2(faces[f]);
					edge0_bc = mesh.getEdgePosition2(faceEdges[0]);
					edge1_bc = mesh.getEdgePosition2(faceEdges[1]);
				}
				else {
					face_bc = mesh.getFaceAverage2(faces[f]);
					edge0_bc = mesh.getEdgeAverage2(faceEdges[0]);
					edge1_bc = mesh.getEdgeAverage2(faceEdges[1]);
				}
				result += std::abs(TwoVector2(face_bc - p0, edge0_bc - p0).determinant()) / 2.0;
				result += std::abs(TwoVector2(face_bc - p0, edge1_bc - p0).determinant()) / 2.0;
			}
			return result;
		}
		else { //assume dimension 3
			const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
			const Buffer<uint>& faces = mesh.getEdgeFaces(i);
			const Buffer<uint> bodies = mesh.getEdgeBodies(i);
			Vector3 p0 = mesh.getNodePosition3(nodes[0]);
			Vector3 p1 = mesh.getNodePosition3(nodes[1]);
			Vector3 edgeVector = p1 - p0;
			Vector3 edge_bc = circumcentric ? mesh.getEdgePosition3(i) : mesh.getEdgeAverage3(i);
			double sum = 0.0;
			for (uint b = 0; b < bodies.size(); ++b) {
				Buffer<uint> bodyFaces = faces.getIntersection(mesh.getBodyFaces(bodies[b]));
				Vector3 face0_bc, face1_bc, body_bc;
				if (circumcentric) {
					face0_bc = mesh.getFacePosition3(bodyFaces[0]);
					face1_bc = mesh.getFacePosition3(bodyFaces[1]);
					body_bc = mesh.getBodyPosition3(bodies[b]);
				}
				else {
					face0_bc = mesh.getFaceAverage3(bodyFaces[0]);
					face1_bc = mesh.getFaceAverage3(bodyFaces[1]);
					body_bc = mesh.getBodyAverage3(bodies[b]);
				}
				sum += 0.5 * TwoVector3(face0_bc - edge_bc, body_bc - edge_bc).len();
				sum += 0.5 * TwoVector3(face1_bc - edge_bc, body_bc - edge_bc).len();
			}
			return sum;
		}
	}

	//computes the length of the boundary of the dual face of node i (assuming 2D mesh)
	double dualFaceBoundaryLength(uint i, const Mesh& mesh, bool circumcentric) {
		const Buffer<uint>& edges = mesh.getNodeEdges(i);
		double result = 0.0;
		for (uint i = 0; i < edges.size(); ++i) {
			result += dualEdgeLength(edges[i], mesh, circumcentric);
		}
		return result;
	}

	//computes the volume of the dual cell of node i (assuming 3D mesh)
	double dualCellVolume(uint i, const Mesh& mesh, bool circumcentric) {
		Vector3 p0 = mesh.getNodePosition3(i);
		Buffer<uint> bodies = mesh.getNodeBodies(i);
		Buffer<uint> faces = mesh.getNodeFaces(i);
		double result = 0.0;
		for (uint b = 0; b < bodies.size(); ++b) {
			Buffer<uint> bodyFaces = faces.getIntersection(mesh.getBodyFaces(bodies[b]));
			Vector3 body_bc, face0_bc, face1_bc, face2_bc, edge01_bc, edge02_bc, edge12_bc;
			if (circumcentric) {
				body_bc = mesh.getBodyPosition3(bodies[b]);
				face0_bc = mesh.getFacePosition3(bodyFaces[0]);
				face1_bc = mesh.getFacePosition3(bodyFaces[1]);
				face2_bc = mesh.getFacePosition3(bodyFaces[2]);
				edge01_bc = mesh.getEdgePosition3(mesh.getFaceIntersection(bodyFaces[0], bodyFaces[1]));
				edge02_bc = mesh.getEdgePosition3(mesh.getFaceIntersection(bodyFaces[0], bodyFaces[2]));
				edge12_bc = mesh.getEdgePosition3(mesh.getFaceIntersection(bodyFaces[1], bodyFaces[2]));
			}
			else {
				body_bc = mesh.getBodyAverage3(bodies[b]);
				face0_bc = mesh.getFaceAverage3(bodyFaces[0]);
				face1_bc = mesh.getFaceAverage3(bodyFaces[1]);
				face2_bc = mesh.getFaceAverage3(bodyFaces[2]);
				edge01_bc = mesh.getEdgeAverage3(mesh.getFaceIntersection(bodyFaces[0], bodyFaces[1]));
				edge02_bc = mesh.getEdgeAverage3(mesh.getFaceIntersection(bodyFaces[0], bodyFaces[2]));
				edge12_bc = mesh.getEdgeAverage3(mesh.getFaceIntersection(bodyFaces[1], bodyFaces[2]));
			}
			result += std::abs(ThreeVector3(body_bc - p0, edge01_bc - p0, face0_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge01_bc - p0, face1_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge02_bc - p0, face0_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge02_bc - p0, face2_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge12_bc - p0, face1_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge12_bc - p0, face2_bc - p0).determinant()) / 6.0;
		}
		return result;
	}

	//computes the volume of the dual cell of node i (assuming 3D Cartesian mesh)
	double dualCellVolumeCartesianMesh(uint i, const Mesh& mesh) {
		Vector3 p0 = mesh.getNodePosition3(i);
		Buffer<uint> bodies = mesh.getNodeBodies(i);
		Buffer<uint> faces = mesh.getNodeFaces(i);
		double result = 0.0;
		for (uint b = 0; b < bodies.size(); ++b) {
			Buffer<uint> bodyFaces = faces.getIntersection(mesh.getBodyFaces(bodies[b]));
			Vector3 body_bc, face0_bc, face1_bc, face2_bc, edge01_bc, edge02_bc, edge12_bc;
			edge01_bc = mesh.getEdgeAverage3(mesh.getFaceIntersection(bodyFaces[0], bodyFaces[1]));
			edge02_bc = mesh.getEdgeAverage3(mesh.getFaceIntersection(bodyFaces[0], bodyFaces[2]));
			edge12_bc = mesh.getEdgeAverage3(mesh.getFaceIntersection(bodyFaces[1], bodyFaces[2]));
			Buffer<uint> face0Nodes = mesh.getFaceNodes(bodyFaces[0]);
			face0_bc = 0.25 * (mesh.getNodePosition3(face0Nodes[0]) + mesh.getNodePosition3(face0Nodes[1]) + mesh.getNodePosition3(face0Nodes[2]) + mesh.getNodePosition3(face0Nodes[3]));
			Buffer<uint> face1Nodes = mesh.getFaceNodes(bodyFaces[1]);
			face1_bc = 0.25 * (mesh.getNodePosition3(face1Nodes[0]) + mesh.getNodePosition3(face1Nodes[1]) + mesh.getNodePosition3(face1Nodes[2]) + mesh.getNodePosition3(face1Nodes[3]));
			Buffer<uint> face2Nodes = mesh.getFaceNodes(bodyFaces[2]);
			face2_bc = 0.25 * (mesh.getNodePosition3(face2Nodes[0]) + mesh.getNodePosition3(face2Nodes[1]) + mesh.getNodePosition3(face2Nodes[2]) + mesh.getNodePosition3(face2Nodes[3]));
			body_bc = edge01_bc + (edge02_bc - p0) + (edge12_bc - p0);
			result += std::abs(ThreeVector3(body_bc - p0, edge01_bc - p0, face0_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge01_bc - p0, face1_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge02_bc - p0, face0_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge02_bc - p0, face2_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge12_bc - p0, face1_bc - p0).determinant()) / 6.0;
			result += std::abs(ThreeVector3(body_bc - p0, edge12_bc - p0, face2_bc - p0).determinant()) / 6.0;
		}
		return result;
	}

	//computes the area of the boundary of the dual cell of node i (assuming 3D mesh)
	double dualCellBoundaryArea(uint i, const Mesh& mesh, bool circumcentric) {
		const Buffer<uint>& edges = mesh.getNodeEdges(i);
		double result = 0.0;
		for (uint i = 0; i < edges.size(); ++i) {
			result += dualFaceArea(edges[i], mesh, circumcentric);
		}
		return result;
	}

	//Create BCC tetrahedral mesh of maximum edge length d such that the elements form a rhombic dodecahedron that is guaranteed to be inside the box [-r, r]^3.
	void createRhombicDodecahedronMesh(BuilderMesh& mesh, double r, double d) {
		mesh.createBccGrid(Vector3(-r, -r, -r), Vector3(r, r, r), d);
		for (uint i = mesh.getNodeSize(); i-- > 0; )
		{
			const Vector4 p = mesh.getNodePosition(i);
			if (p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.x > r - 1e-5) mesh.removeNode(i);
		}
	}

	//create triangle mesh in [0, L] x [0, T] with the requested number of steps in x- and t-directions
	void createSpacetimeMesh(Mesh& mesh, double L, double T, const uint xsteps, const uint tsteps) {
		const double dx = L / xsteps;
		const double dt = T / tsteps;
		for (uint i = 0; i < tsteps + 1; ++i) {
			for (uint j = 0; j < xsteps + 1; ++j) {
				mesh.addNode(Vector4(j * dx, i * dt, 0.0, 0.0));
			}
		}
		for (uint i = 0; i < tsteps; ++i) {
			for (uint j = 0; j < xsteps; ++j) {
				uint node00 = i * (xsteps + 1) + j;
				uint node10 = node00 + 1;
				uint node01 = (i + 1) * (xsteps + 1) + j;
				uint node11 = node01 + 1;
				Buffer<uint> edges(3);
				edges[0] = mesh.addEdge(node00, node10);
				edges[1] = mesh.addEdge(node01, node10);
				edges[2] = mesh.addEdge(node00, node01);
				mesh.addFace(edges);
				edges[0] = mesh.addEdge(node01, node11);
				edges[2] = mesh.addEdge(node10, node11);
				mesh.addFace(edges);
			}
		}
	}

	//create tetrahedral mesh in [0, Lx] x [0, Ly] x [0, T] with the requested number of steps in x-, y-, and t-directions
	void createSpacetimeMesh(Mesh& mesh, double Lx, double Ly, double T, uint xsteps, uint ysteps, uint tsteps) {
		const double dx = Lx / xsteps;
		const double dy = Ly / ysteps;
		const double dt = T / tsteps;
		for (uint i = 0; i < tsteps + 1; ++i) {
			for (uint j = 0; j < ysteps + 1; ++j) {
				for (uint k = 0; k < xsteps + 1; ++k) {
					mesh.addNode(Vector4(k * dx, j * dy, i * dt, 0.0));
				}
			}
		}
		for (uint i = 0; i < tsteps; ++i) {
			for (uint j = 0; j < ysteps; ++j) {
				for (uint k = 0; k < xsteps; ++k) {
					uint node000 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node100 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node010 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
					uint node110 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
					uint node001 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node101 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node011 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
					uint node111 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
					uint edge_000_100 = mesh.addEdge(node000, node100);
					uint edge_000_010 = mesh.addEdge(node000, node010);
					uint edge_000_001 = mesh.addEdge(node000, node001);
					uint edge_000_110 = mesh.addEdge(node000, node110);
					uint edge_001_101 = mesh.addEdge(node001, node101);
					uint edge_001_011 = mesh.addEdge(node001, node011);
					uint edge_001_100 = mesh.addEdge(node001, node100);
					uint edge_001_010 = mesh.addEdge(node001, node010);
					uint edge_001_111 = mesh.addEdge(node001, node111);
					uint edge_001_110 = mesh.addEdge(node001, node110);
					uint edge_100_101 = mesh.addEdge(node100, node101);
					uint edge_101_111 = mesh.addEdge(node101, node111);
					uint edge_010_011 = mesh.addEdge(node010, node011);
					uint edge_100_110 = mesh.addEdge(node100, node110);
					uint edge_110_111 = mesh.addEdge(node110, node111);
					uint edge_010_110 = mesh.addEdge(node010, node110);
					uint edge_011_110 = mesh.addEdge(node011, node110);
					uint edge_011_111 = mesh.addEdge(node011, node111);
					uint edge_101_110 = mesh.addEdge(node101, node110);
					Buffer<uint> edges(3);
					edges[0] = edge_000_100;
					edges[1] = edge_100_110;
					edges[2] = edge_000_110;
					uint face_000_100_110 = mesh.addFace(edges);
					edges[0] = edge_000_010;
					edges[1] = edge_000_110;
					edges[2] = edge_010_110;
					uint face_000_110_010 = mesh.addFace(edges);
					edges[0] = edge_000_001;
					edges[1] = edge_001_010;
					edges[2] = edge_000_010;
					uint face_000_001_010 = mesh.addFace(edges);
					edges[0] = edge_010_011;
					edges[1] = edge_001_011;
					edges[2] = edge_001_010;
					uint face_010_011_001 = mesh.addFace(edges);
					edges[0] = edge_000_100;
					edges[1] = edge_001_100;
					edges[2] = edge_000_001;
					uint face_000_100_001 = mesh.addFace(edges);
					edges[0] = edge_100_101;
					edges[1] = edge_001_101;
					edges[2] = edge_001_100;
					uint face_100_101_001 = mesh.addFace(edges);
					edges[0] = edge_100_101;
					edges[1] = edge_101_110;
					edges[2] = edge_100_110;
					uint face_100_110_101 = mesh.addFace(edges);
					edges[0] = edge_110_111;
					edges[1] = edge_101_111;
					edges[2] = edge_101_110;
					uint face_110_111_101 = mesh.addFace(edges);
					edges[0] = edge_010_110;
					edges[1] = edge_011_110;
					edges[2] = edge_010_011;
					uint face_010_110_011 = mesh.addFace(edges);
					edges[0] = edge_110_111;
					edges[1] = edge_011_111;
					edges[2] = edge_011_110;
					uint face_110_111_011 = mesh.addFace(edges);
					edges[0] = edge_001_101;
					edges[1] = edge_101_111;
					edges[2] = edge_001_111;
					uint face_001_101_111 = mesh.addFace(edges);
					edges[0] = edge_001_011;
					edges[1] = edge_001_111;
					edges[2] = edge_011_111;
					uint face_001_111_011 = mesh.addFace(edges);
					edges[0] = edge_001_100;
					edges[1] = edge_100_110;
					edges[2] = edge_001_110;
					uint face_001_100_110 = mesh.addFace(edges);
					edges[0] = edge_001_110;
					edges[1] = edge_011_110;
					edges[2] = edge_001_011;
					uint face_001_110_011 = mesh.addFace(edges);
					edges[0] = edge_001_101;
					edges[1] = edge_101_110;
					edges[2] = edge_001_110;
					uint face_001_101_110 = mesh.addFace(edges);
					edges[0] = edge_010_110;
					edges[1] = edge_001_110;
					edges[2] = edge_001_010;
					uint face_010_110_001 = mesh.addFace(edges);
					edges[0] = edge_000_110;
					edges[1] = edge_001_110;
					edges[2] = edge_000_001;
					uint face_000_110_001 = mesh.addFace(edges);
					edges[0] = edge_001_110;
					edges[1] = edge_110_111;
					edges[2] = edge_001_111;
					uint face_001_110_111 = mesh.addFace(edges);
					Buffer<uint> faces(4);
					faces[0] = face_000_100_110;
					faces[1] = face_000_100_001;
					faces[2] = face_000_110_001;
					faces[3] = face_001_100_110;
					mesh.addBody(faces);
					faces[0] = face_100_101_001;
					faces[1] = face_100_110_101;
					faces[2] = face_001_101_110;
					faces[3] = face_001_100_110;
					mesh.addBody(faces);
					faces[0] = face_110_111_101;
					faces[1] = face_001_101_110;
					faces[2] = face_001_101_111;
					faces[3] = face_001_110_111;
					mesh.addBody(faces);
					faces[0] = face_000_110_010;
					faces[1] = face_000_001_010;
					faces[2] = face_000_110_001;
					faces[3] = face_010_110_001;
					mesh.addBody(faces);
					faces[0] = face_010_011_001;
					faces[1] = face_010_110_011;
					faces[2] = face_001_110_011;
					faces[3] = face_010_110_001;
					mesh.addBody(faces);
					faces[0] = face_001_111_011;
					faces[1] = face_110_111_011;
					faces[2] = face_001_110_011;
					faces[3] = face_001_110_111;
					mesh.addBody(faces);
				}
			}
		}
	}

	//create tetrahedral mesh in [0, Lx] x [0, Ly] x [0, T] with the requested number of steps in x-, y-, and t-directions, with an extra node in each box
	void createSpacetimeMeshExtraNode(Mesh& mesh, double Lx, double Ly, double T, uint xsteps, uint ysteps, uint tsteps) {
		const double dx = Lx / xsteps;
		const double dy = Ly / ysteps;
		const double dt = T / tsteps;
		for (uint i = 0; i < tsteps + 1; ++i) {
			for (uint j = 0; j < ysteps + 1; ++j) {
				for (uint k = 0; k < xsteps + 1; ++k) {
					mesh.addNode(Vector4(k * dx, j * dy, i * dt, 0.0));
				}
			}
		}
		for (uint i = 0; i < tsteps; ++i) {
			for (uint j = 0; j < ysteps; ++j) {
				for (uint k = 0; k < xsteps; ++k) {
					uint node000 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node100 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node010 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
					uint node110 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
					uint node001 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node101 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node011 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
					uint node111 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
					uint nodeC = mesh.addNode(Vector4((k + 0.5) * dx, (j + 0.5) * dy, (i + 0.5) * dt, 0.0));
					uint edge_000_100 = mesh.addEdge(node000, node100);
					uint edge_000_010 = mesh.addEdge(node000, node010);
					uint edge_000_001 = mesh.addEdge(node000, node001);
					uint edge_000_110 = mesh.addEdge(node000, node110);
					uint edge_001_101 = mesh.addEdge(node001, node101);
					uint edge_001_011 = mesh.addEdge(node001, node011);
					uint edge_001_100 = mesh.addEdge(node001, node100);
					uint edge_001_010 = mesh.addEdge(node001, node010);
					uint edge_001_111 = mesh.addEdge(node001, node111);
					uint edge_100_101 = mesh.addEdge(node100, node101);
					uint edge_101_111 = mesh.addEdge(node101, node111);
					uint edge_010_011 = mesh.addEdge(node010, node011);
					uint edge_100_110 = mesh.addEdge(node100, node110);
					uint edge_110_111 = mesh.addEdge(node110, node111);
					uint edge_010_110 = mesh.addEdge(node010, node110);
					uint edge_011_110 = mesh.addEdge(node011, node110);
					uint edge_011_111 = mesh.addEdge(node011, node111);
					uint edge_101_110 = mesh.addEdge(node101, node110);
					uint edge_000_C = mesh.addEdge(node000, nodeC);
					uint edge_100_C = mesh.addEdge(node100, nodeC);
					uint edge_010_C = mesh.addEdge(node010, nodeC);
					uint edge_110_C = mesh.addEdge(node110, nodeC);
					uint edge_001_C = mesh.addEdge(node001, nodeC);
					uint edge_101_C = mesh.addEdge(node101, nodeC);
					uint edge_011_C = mesh.addEdge(node011, nodeC);
					uint edge_111_C = mesh.addEdge(node111, nodeC);
					Buffer<uint> edges(3);
					edges[0] = edge_000_100;
					edges[1] = edge_100_110;
					edges[2] = edge_000_110;
					uint face_000_100_110 = mesh.addFace(edges);
					edges[0] = edge_000_010;
					edges[1] = edge_000_110;
					edges[2] = edge_010_110;
					uint face_000_110_010 = mesh.addFace(edges);
					edges[0] = edge_000_001;
					edges[1] = edge_001_010;
					edges[2] = edge_000_010;
					uint face_000_001_010 = mesh.addFace(edges);
					edges[0] = edge_010_011;
					edges[1] = edge_001_011;
					edges[2] = edge_001_010;
					uint face_010_011_001 = mesh.addFace(edges);
					edges[0] = edge_000_100;
					edges[1] = edge_001_100;
					edges[2] = edge_000_001;
					uint face_000_100_001 = mesh.addFace(edges);
					edges[0] = edge_100_101;
					edges[1] = edge_001_101;
					edges[2] = edge_001_100;
					uint face_100_101_001 = mesh.addFace(edges);
					edges[0] = edge_100_101;
					edges[1] = edge_101_110;
					edges[2] = edge_100_110;
					uint face_100_110_101 = mesh.addFace(edges);
					edges[0] = edge_110_111;
					edges[1] = edge_101_111;
					edges[2] = edge_101_110;
					uint face_110_111_101 = mesh.addFace(edges);
					edges[0] = edge_010_110;
					edges[1] = edge_011_110;
					edges[2] = edge_010_011;
					uint face_010_110_011 = mesh.addFace(edges);
					edges[0] = edge_110_111;
					edges[1] = edge_011_111;
					edges[2] = edge_011_110;
					uint face_110_111_011 = mesh.addFace(edges);
					edges[0] = edge_001_101;
					edges[1] = edge_101_111;
					edges[2] = edge_001_111;
					uint face_001_101_111 = mesh.addFace(edges);
					edges[0] = edge_001_011;
					edges[1] = edge_001_111;
					edges[2] = edge_011_111;
					uint face_001_111_011 = mesh.addFace(edges);

					edges[0] = edge_000_100;
					edges[1] = edge_000_C;
					edges[2] = edge_100_C;
					uint face_000_100_C = mesh.addFace(edges);
					edges[0] = edge_000_010;
					edges[1] = edge_000_C;
					edges[2] = edge_010_C;
					uint face_000_010_C = mesh.addFace(edges);
					edges[0] = edge_000_001;
					edges[1] = edge_000_C;
					edges[2] = edge_001_C;
					uint face_000_001_C = mesh.addFace(edges);
					edges[0] = edge_000_110;
					edges[1] = edge_000_C;
					edges[2] = edge_110_C;
					uint face_000_110_C = mesh.addFace(edges);
					edges[0] = edge_001_101;
					edges[1] = edge_001_C;
					edges[2] = edge_101_C;
					uint face_001_101_C = mesh.addFace(edges);
					edges[0] = edge_001_011;
					edges[1] = edge_001_C;
					edges[2] = edge_011_C;
					uint face_001_011_C = mesh.addFace(edges);
					edges[0] = edge_001_100;
					edges[1] = edge_001_C;
					edges[2] = edge_100_C;
					uint face_001_100_C = mesh.addFace(edges);
					edges[0] = edge_001_010;
					edges[1] = edge_001_C;
					edges[2] = edge_010_C;
					uint face_001_010_C = mesh.addFace(edges);
					edges[0] = edge_001_111;
					edges[1] = edge_001_C;
					edges[2] = edge_111_C;
					uint face_001_111_C = mesh.addFace(edges);
					edges[0] = edge_100_101;
					edges[1] = edge_100_C;
					edges[2] = edge_101_C;
					uint face_100_101_C = mesh.addFace(edges);
					edges[0] = edge_101_111;
					edges[1] = edge_101_C;
					edges[2] = edge_111_C;
					uint face_101_111_C = mesh.addFace(edges);
					edges[0] = edge_010_011;
					edges[1] = edge_010_C;
					edges[2] = edge_011_C;
					uint face_010_011_C = mesh.addFace(edges);
					edges[0] = edge_100_110;
					edges[1] = edge_100_C;
					edges[2] = edge_110_C;
					uint face_100_110_C = mesh.addFace(edges);
					edges[0] = edge_110_111;
					edges[1] = edge_110_C;
					edges[2] = edge_111_C;
					uint face_110_111_C = mesh.addFace(edges);
					edges[0] = edge_010_110;
					edges[1] = edge_010_C;
					edges[2] = edge_110_C;
					uint face_010_110_C = mesh.addFace(edges);
					edges[0] = edge_011_110;
					edges[1] = edge_011_C;
					edges[2] = edge_110_C;
					uint face_011_110_C = mesh.addFace(edges);
					edges[0] = edge_011_111;
					edges[1] = edge_011_C;
					edges[2] = edge_111_C;
					uint face_011_111_C = mesh.addFace(edges);
					edges[0] = edge_101_110;
					edges[1] = edge_101_C;
					edges[2] = edge_110_C;
					uint face_101_110_C = mesh.addFace(edges);				
					
					Buffer<uint> faces(4);
					faces[0] = face_000_100_110;
					faces[1] = face_000_100_C;
					faces[2] = face_100_110_C;
					faces[3] = face_000_110_C;
					mesh.addBody(faces);

					faces[0] = face_000_110_010;
					faces[1] = face_000_010_C;
					faces[2] = face_010_110_C;
					faces[3] = face_000_110_C;
					mesh.addBody(faces);
					
					faces[0] = face_000_100_001;
					faces[1] = face_000_100_C;
					faces[2] = face_000_001_C;
					faces[3] = face_001_100_C;
					mesh.addBody(faces);

					faces[0] = face_100_101_001;
					faces[1] = face_100_101_C;
					faces[2] = face_001_101_C;
					faces[3] = face_001_100_C;
					mesh.addBody(faces);

					faces[0] = face_000_001_010;
					faces[1] = face_000_001_C;
					faces[2] = face_000_010_C;
					faces[3] = face_001_010_C;
					mesh.addBody(faces);

					faces[0] = face_010_011_001;
					faces[1] = face_010_011_C;
					faces[2] = face_001_011_C;
					faces[3] = face_001_010_C;
					mesh.addBody(faces);



					faces[0] = face_100_110_101;
					faces[1] = face_100_101_C;
					faces[2] = face_100_110_C;
					faces[3] = face_101_110_C;
					mesh.addBody(faces);

					faces[0] = face_110_111_101;
					faces[1] = face_110_111_C;
					faces[2] = face_101_111_C;
					faces[3] = face_101_110_C;
					mesh.addBody(faces);

					faces[0] = face_010_110_011;
					faces[1] = face_010_110_C;
					faces[2] = face_010_011_C;
					faces[3] = face_011_110_C;
					mesh.addBody(faces);

					faces[0] = face_110_111_011;
					faces[1] = face_110_111_C;
					faces[2] = face_011_111_C;
					faces[3] = face_011_110_C;
					mesh.addBody(faces);

					faces[0] = face_001_101_111;
					faces[1] = face_001_101_C;
					faces[2] = face_101_111_C;
					faces[3] = face_001_111_C;
					mesh.addBody(faces);

					faces[0] = face_001_111_011;
					faces[1] = face_001_011_C;
					faces[2] = face_011_111_C;
					faces[3] = face_001_111_C;
					mesh.addBody(faces);
				}
			}
		}
	}

	//create tetrahedral mesh in [0, Lx] x [0, Ly] x [0, T] with the requested number of steps in x-, y-, and t-directions, with extra nodes in center and faces of each box
	void createSpacetimeMeshExtraNodes(Mesh& mesh, double Lx, double Ly, double T, uint xsteps, uint ysteps, uint tsteps) {
		const double dx = Lx / xsteps;
		const double dy = Ly / ysteps;
		const double dt = T / tsteps;
		uint total_corner_nodes = (tsteps + 1) * (ysteps + 1) * (xsteps + 1);
		for (uint i = 0; i < tsteps + 1; ++i) {
			for (uint j = 0; j < ysteps + 1; ++j) {
				for (uint k = 0; k < xsteps + 1; ++k) {
					mesh.addNode(Vector4(k * dx, j * dy, i * dt, 0.0));
				}
			}
		}
		for (uint i = 0; i < tsteps + 1; ++i) {
			for (uint j = 0; j < ysteps; ++j) {
				for (uint k = 0; k < xsteps; ++k) {
					mesh.addNode(Vector4((k + 0.5) * dx, (j + 0.5) * dy, i * dt, 0.0));
				}
			}
		}
		for (uint i = 0; i < tsteps; ++i) {
			for (uint j = 0; j < ysteps + 1; ++j) {
				for (uint k = 0; k < xsteps; ++k) {
					mesh.addNode(Vector4((k + 0.5) * dx, j * dy, (i + 0.5) * dt, 0.0));
				}
			}
		}
		for (uint i = 0; i < tsteps; ++i) {
			for (uint j = 0; j < ysteps; ++j) {
				for (uint k = 0; k < xsteps + 1; ++k) {
					mesh.addNode(Vector4(k * dx, (j + 0.5) * dy, (i + 0.5) * dt, 0.0));
				}
			}
		}
		for (uint i = 0; i < tsteps; ++i) {
			for (uint j = 0; j < ysteps; ++j) {
				for (uint k = 0; k < xsteps; ++k) {
					uint node_000 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node_100 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node_010 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
					uint node_110 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
					uint node_001 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node_101 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node_011 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
					uint node_111 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
					uint node_c = mesh.addNode(Vector4((k + 0.5) * dx, (j + 0.5) * dy, (i + 0.5) * dt, 0.0));
					uint node_f_down = total_corner_nodes + i * ysteps * xsteps + j * xsteps + k;
					uint node_f_west = total_corner_nodes + (tsteps + 1) * ysteps * xsteps + tsteps * (ysteps + 1) * xsteps + i * ysteps * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node_f_south = total_corner_nodes + (tsteps + 1) * ysteps * xsteps + i * (ysteps + 1) * xsteps + j * xsteps + k;
					uint node_f_east = total_corner_nodes + (tsteps + 1) * ysteps * xsteps + tsteps * (ysteps + 1) * xsteps + i * ysteps * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node_f_north = total_corner_nodes + (tsteps + 1) * ysteps * xsteps + i * (ysteps + 1) * xsteps + (j + 1) * xsteps + k;
					uint node_f_up = total_corner_nodes + (i + 1) * ysteps * xsteps + j * xsteps + k;

					uint edge_000_100 = mesh.addEdge(node_000, node_100);
					uint edge_010_110 = mesh.addEdge(node_010, node_110);
					uint edge_000_010 = mesh.addEdge(node_000, node_010);
					uint edge_100_110 = mesh.addEdge(node_100, node_110);
					uint edge_000_f_down = mesh.addEdge(node_000, node_f_down);
					uint edge_100_f_down = mesh.addEdge(node_100, node_f_down);
					uint edge_010_f_down = mesh.addEdge(node_010, node_f_down);
					uint edge_110_f_down = mesh.addEdge(node_110, node_f_down);
					uint edge_000_c = mesh.addEdge(node_000, node_c);
					uint edge_100_c = mesh.addEdge(node_100, node_c);
					uint edge_010_c = mesh.addEdge(node_010, node_c);
					uint edge_110_c = mesh.addEdge(node_110, node_c);
					uint edge_f_down_c = mesh.addEdge(node_f_down, node_c);
					uint edge_001_011 = mesh.addEdge(node_001, node_011);
					uint edge_000_001 = mesh.addEdge(node_000, node_001);
					uint edge_010_011 = mesh.addEdge(node_010, node_011);
					uint edge_000_f_west = mesh.addEdge(node_000, node_f_west);
					uint edge_010_f_west = mesh.addEdge(node_010, node_f_west);
					uint edge_001_f_west = mesh.addEdge(node_001, node_f_west);
					uint edge_011_f_west = mesh.addEdge(node_011, node_f_west);
					uint edge_001_c = mesh.addEdge(node_001, node_c);
					uint edge_011_c = mesh.addEdge(node_011, node_c);
					uint edge_f_west_c = mesh.addEdge(node_f_west, node_c);
					uint edge_001_101 = mesh.addEdge(node_001, node_101);
					uint edge_100_101 = mesh.addEdge(node_100, node_101);
					uint edge_000_f_south = mesh.addEdge(node_000, node_f_south);
					uint edge_100_f_south = mesh.addEdge(node_100, node_f_south);
					uint edge_001_f_south = mesh.addEdge(node_001, node_f_south);
					uint edge_101_f_south = mesh.addEdge(node_101, node_f_south);
					uint edge_101_c = mesh.addEdge(node_101, node_c);
					uint edge_f_south_c = mesh.addEdge(node_f_south, node_c);
					uint edge_101_111 = mesh.addEdge(node_101, node_111);
					uint edge_110_111 = mesh.addEdge(node_110, node_111);
					uint edge_100_f_east = mesh.addEdge(node_100, node_f_east);
					uint edge_110_f_east = mesh.addEdge(node_110, node_f_east);
					uint edge_101_f_east = mesh.addEdge(node_101, node_f_east);
					uint edge_111_f_east = mesh.addEdge(node_111, node_f_east);
					uint edge_111_c = mesh.addEdge(node_111, node_c);
					uint edge_f_east_c = mesh.addEdge(node_f_east, node_c);
					uint edge_011_111 = mesh.addEdge(node_011, node_111);
					uint edge_010_f_north = mesh.addEdge(node_010, node_f_north);
					uint edge_110_f_north = mesh.addEdge(node_110, node_f_north);
					uint edge_011_f_north = mesh.addEdge(node_011, node_f_north);
					uint edge_111_f_north = mesh.addEdge(node_111, node_f_north);
					uint edge_f_north_c = mesh.addEdge(node_f_north, node_c);
					uint edge_001_f_up = mesh.addEdge(node_001, node_f_up);
					uint edge_101_f_up = mesh.addEdge(node_101, node_f_up);
					uint edge_011_f_up = mesh.addEdge(node_011, node_f_up);
					uint edge_111_f_up = mesh.addEdge(node_111, node_f_up);
					uint edge_f_up_c = mesh.addEdge(node_f_up, node_c);

					Buffer<uint> edges(3);
					Buffer<uint> faces(4);

					edges[0] = edge_000_100;
					edges[1] = edge_100_f_down;
					edges[2] = edge_000_f_down;
					uint face_edge_000_100_edge_100_f_down_edge_000_f_down = mesh.addFace(edges);
					edges[0] = edge_100_110;
					edges[1] = edge_100_f_down;
					edges[2] = edge_110_f_down;
					uint face_edge_100_110_edge_100_f_down_edge_110_f_down = mesh.addFace(edges);
					edges[0] = edge_010_110;
					edges[1] = edge_010_f_down;
					edges[2] = edge_110_f_down;
					uint face_edge_010_110_edge_010_f_down_edge_110_f_down = mesh.addFace(edges);
					edges[0] = edge_000_010;
					edges[1] = edge_000_f_down;
					edges[2] = edge_010_f_down;
					uint face_edge_000_010_edge_000_f_down_edge_010_f_down = mesh.addFace(edges);
					edges[0] = edge_000_100;
					edges[1] = edge_000_c;
					edges[2] = edge_100_c;
					uint face_edge_000_100_edge_000_c_edge_100_c = mesh.addFace(edges);
					edges[0] = edge_100_110;
					edges[1] = edge_100_c;
					edges[2] = edge_110_c;
					uint face_edge_100_110_edge_100_c_edge_110_c = mesh.addFace(edges);
					edges[0] = edge_010_110;
					edges[1] = edge_010_c;
					edges[2] = edge_110_c;
					uint face_edge_010_110_edge_010_c_edge_110_c = mesh.addFace(edges);
					edges[0] = edge_000_010;
					edges[1] = edge_000_c;
					edges[2] = edge_010_c;
					uint face_edge_000_010_edge_000_c_edge_010_c = mesh.addFace(edges);
					edges[0] = edge_000_f_down;
					edges[1] = edge_000_c;
					edges[2] = edge_f_down_c;
					uint face_edge_000_f_down_edge_000_c_edge_f_down_c = mesh.addFace(edges);
					edges[0] = edge_100_f_down;
					edges[1] = edge_100_c;
					edges[2] = edge_f_down_c;
					uint face_edge_100_f_down_edge_100_c_edge_f_down_c = mesh.addFace(edges);
					edges[0] = edge_010_f_down;
					edges[1] = edge_010_c;
					edges[2] = edge_f_down_c;
					uint face_edge_010_f_down_edge_010_c_edge_f_down_c = mesh.addFace(edges);
					edges[0] = edge_110_f_down;
					edges[1] = edge_110_c;
					edges[2] = edge_f_down_c;
					uint face_edge_110_f_down_edge_110_c_edge_f_down_c = mesh.addFace(edges);
					edges[0] = edge_000_010;
					edges[1] = edge_010_f_west;
					edges[2] = edge_000_f_west;
					uint face_edge_000_010_edge_010_f_west_edge_000_f_west = mesh.addFace(edges);
					edges[0] = edge_010_011;
					edges[1] = edge_010_f_west;
					edges[2] = edge_011_f_west;
					uint face_edge_010_011_edge_010_f_west_edge_011_f_west = mesh.addFace(edges);
					edges[0] = edge_001_011;
					edges[1] = edge_001_f_west;
					edges[2] = edge_011_f_west;
					uint face_edge_001_011_edge_001_f_west_edge_011_f_west = mesh.addFace(edges);
					edges[0] = edge_000_001;
					edges[1] = edge_000_f_west;
					edges[2] = edge_001_f_west;
					uint face_edge_000_001_edge_000_f_west_edge_001_f_west = mesh.addFace(edges);
					edges[0] = edge_010_011;
					edges[1] = edge_010_c;
					edges[2] = edge_011_c;
					uint face_edge_010_011_edge_010_c_edge_011_c = mesh.addFace(edges);
					edges[0] = edge_001_011;
					edges[1] = edge_001_c;
					edges[2] = edge_011_c;
					uint face_edge_001_011_edge_001_c_edge_011_c = mesh.addFace(edges);
					edges[0] = edge_000_001;
					edges[1] = edge_000_c;
					edges[2] = edge_001_c;
					uint face_edge_000_001_edge_000_c_edge_001_c = mesh.addFace(edges);
					edges[0] = edge_000_f_west;
					edges[1] = edge_000_c;
					edges[2] = edge_f_west_c;
					uint face_edge_000_f_west_edge_000_c_edge_f_west_c = mesh.addFace(edges);
					edges[0] = edge_010_f_west;
					edges[1] = edge_010_c;
					edges[2] = edge_f_west_c;
					uint face_edge_010_f_west_edge_010_c_edge_f_west_c = mesh.addFace(edges);
					edges[0] = edge_001_f_west;
					edges[1] = edge_001_c;
					edges[2] = edge_f_west_c;
					uint face_edge_001_f_west_edge_001_c_edge_f_west_c = mesh.addFace(edges);
					edges[0] = edge_011_f_west;
					edges[1] = edge_011_c;
					edges[2] = edge_f_west_c;
					uint face_edge_011_f_west_edge_011_c_edge_f_west_c = mesh.addFace(edges);
					edges[0] = edge_000_100;
					edges[1] = edge_100_f_south;
					edges[2] = edge_000_f_south;
					uint face_edge_000_100_edge_100_f_south_edge_000_f_south = mesh.addFace(edges);
					edges[0] = edge_100_101;
					edges[1] = edge_100_f_south;
					edges[2] = edge_101_f_south;
					uint face_edge_100_101_edge_100_f_south_edge_101_f_south = mesh.addFace(edges);
					edges[0] = edge_001_101;
					edges[1] = edge_001_f_south;
					edges[2] = edge_101_f_south;
					uint face_edge_001_101_edge_001_f_south_edge_101_f_south = mesh.addFace(edges);
					edges[0] = edge_000_001;
					edges[1] = edge_000_f_south;
					edges[2] = edge_001_f_south;
					uint face_edge_000_001_edge_000_f_south_edge_001_f_south = mesh.addFace(edges);
					edges[0] = edge_100_101;
					edges[1] = edge_100_c;
					edges[2] = edge_101_c;
					uint face_edge_100_101_edge_100_c_edge_101_c = mesh.addFace(edges);
					edges[0] = edge_001_101;
					edges[1] = edge_001_c;
					edges[2] = edge_101_c;
					uint face_edge_001_101_edge_001_c_edge_101_c = mesh.addFace(edges);
					edges[0] = edge_000_f_south;
					edges[1] = edge_000_c;
					edges[2] = edge_f_south_c;
					uint face_edge_000_f_south_edge_000_c_edge_f_south_c = mesh.addFace(edges);
					edges[0] = edge_100_f_south;
					edges[1] = edge_100_c;
					edges[2] = edge_f_south_c;
					uint face_edge_100_f_south_edge_100_c_edge_f_south_c = mesh.addFace(edges);
					edges[0] = edge_001_f_south;
					edges[1] = edge_001_c;
					edges[2] = edge_f_south_c;
					uint face_edge_001_f_south_edge_001_c_edge_f_south_c = mesh.addFace(edges);
					edges[0] = edge_101_f_south;
					edges[1] = edge_101_c;
					edges[2] = edge_f_south_c;
					uint face_edge_101_f_south_edge_101_c_edge_f_south_c = mesh.addFace(edges);
					edges[0] = edge_100_110;
					edges[1] = edge_110_f_east;
					edges[2] = edge_100_f_east;
					uint face_edge_100_110_edge_110_f_east_edge_100_f_east = mesh.addFace(edges);
					edges[0] = edge_110_111;
					edges[1] = edge_110_f_east;
					edges[2] = edge_111_f_east;
					uint face_edge_110_111_edge_110_f_east_edge_111_f_east = mesh.addFace(edges);
					edges[0] = edge_101_111;
					edges[1] = edge_101_f_east;
					edges[2] = edge_111_f_east;
					uint face_edge_101_111_edge_101_f_east_edge_111_f_east = mesh.addFace(edges);
					edges[0] = edge_100_101;
					edges[1] = edge_100_f_east;
					edges[2] = edge_101_f_east;
					uint face_edge_100_101_edge_100_f_east_edge_101_f_east = mesh.addFace(edges);
					edges[0] = edge_110_111;
					edges[1] = edge_110_c;
					edges[2] = edge_111_c;
					uint face_edge_110_111_edge_110_c_edge_111_c = mesh.addFace(edges);
					edges[0] = edge_101_111;
					edges[1] = edge_101_c;
					edges[2] = edge_111_c;
					uint face_edge_101_111_edge_101_c_edge_111_c = mesh.addFace(edges);
					edges[0] = edge_100_f_east;
					edges[1] = edge_100_c;
					edges[2] = edge_f_east_c;
					uint face_edge_100_f_east_edge_100_c_edge_f_east_c = mesh.addFace(edges);
					edges[0] = edge_110_f_east;
					edges[1] = edge_110_c;
					edges[2] = edge_f_east_c;
					uint face_edge_110_f_east_edge_110_c_edge_f_east_c = mesh.addFace(edges);
					edges[0] = edge_101_f_east;
					edges[1] = edge_101_c;
					edges[2] = edge_f_east_c;
					uint face_edge_101_f_east_edge_101_c_edge_f_east_c = mesh.addFace(edges);
					edges[0] = edge_111_f_east;
					edges[1] = edge_111_c;
					edges[2] = edge_f_east_c;
					uint face_edge_111_f_east_edge_111_c_edge_f_east_c = mesh.addFace(edges);
					edges[0] = edge_010_110;
					edges[1] = edge_110_f_north;
					edges[2] = edge_010_f_north;
					uint face_edge_010_110_edge_110_f_north_edge_010_f_north = mesh.addFace(edges);
					edges[0] = edge_110_111;
					edges[1] = edge_110_f_north;
					edges[2] = edge_111_f_north;
					uint face_edge_110_111_edge_110_f_north_edge_111_f_north = mesh.addFace(edges);
					edges[0] = edge_011_111;
					edges[1] = edge_011_f_north;
					edges[2] = edge_111_f_north;
					uint face_edge_011_111_edge_011_f_north_edge_111_f_north = mesh.addFace(edges);
					edges[0] = edge_010_011;
					edges[1] = edge_010_f_north;
					edges[2] = edge_011_f_north;
					uint face_edge_010_011_edge_010_f_north_edge_011_f_north = mesh.addFace(edges);
					edges[0] = edge_011_111;
					edges[1] = edge_011_c;
					edges[2] = edge_111_c;
					uint face_edge_011_111_edge_011_c_edge_111_c = mesh.addFace(edges);
					edges[0] = edge_010_f_north;
					edges[1] = edge_010_c;
					edges[2] = edge_f_north_c;
					uint face_edge_010_f_north_edge_010_c_edge_f_north_c = mesh.addFace(edges);
					edges[0] = edge_110_f_north;
					edges[1] = edge_110_c;
					edges[2] = edge_f_north_c;
					uint face_edge_110_f_north_edge_110_c_edge_f_north_c = mesh.addFace(edges);
					edges[0] = edge_011_f_north;
					edges[1] = edge_011_c;
					edges[2] = edge_f_north_c;
					uint face_edge_011_f_north_edge_011_c_edge_f_north_c = mesh.addFace(edges);
					edges[0] = edge_111_f_north;
					edges[1] = edge_111_c;
					edges[2] = edge_f_north_c;
					uint face_edge_111_f_north_edge_111_c_edge_f_north_c = mesh.addFace(edges);
					edges[0] = edge_001_101;
					edges[1] = edge_101_f_up;
					edges[2] = edge_001_f_up;
					uint face_edge_001_101_edge_101_f_up_edge_001_f_up = mesh.addFace(edges);
					edges[0] = edge_101_111;
					edges[1] = edge_101_f_up;
					edges[2] = edge_111_f_up;
					uint face_edge_101_111_edge_101_f_up_edge_111_f_up = mesh.addFace(edges);
					edges[0] = edge_011_111;
					edges[1] = edge_011_f_up;
					edges[2] = edge_111_f_up;
					uint face_edge_011_111_edge_011_f_up_edge_111_f_up = mesh.addFace(edges);
					edges[0] = edge_001_011;
					edges[1] = edge_001_f_up;
					edges[2] = edge_011_f_up;
					uint face_edge_001_011_edge_001_f_up_edge_011_f_up = mesh.addFace(edges);
					edges[0] = edge_001_f_up;
					edges[1] = edge_001_c;
					edges[2] = edge_f_up_c;
					uint face_edge_001_f_up_edge_001_c_edge_f_up_c = mesh.addFace(edges);
					edges[0] = edge_101_f_up;
					edges[1] = edge_101_c;
					edges[2] = edge_f_up_c;
					uint face_edge_101_f_up_edge_101_c_edge_f_up_c = mesh.addFace(edges);
					edges[0] = edge_011_f_up;
					edges[1] = edge_011_c;
					edges[2] = edge_f_up_c;
					uint face_edge_011_f_up_edge_011_c_edge_f_up_c = mesh.addFace(edges);
					edges[0] = edge_111_f_up;
					edges[1] = edge_111_c;
					edges[2] = edge_f_up_c;
					uint face_edge_111_f_up_edge_111_c_edge_f_up_c = mesh.addFace(edges);

					faces[0] = face_edge_000_100_edge_100_f_down_edge_000_f_down;
					faces[1] = face_edge_000_100_edge_000_c_edge_100_c;
					faces[2] = face_edge_000_f_down_edge_000_c_edge_f_down_c;
					faces[3] = face_edge_100_f_down_edge_100_c_edge_f_down_c;
					mesh.addBody(faces);
					faces[0] = face_edge_100_110_edge_100_f_down_edge_110_f_down;
					faces[1] = face_edge_100_110_edge_100_c_edge_110_c;
					faces[2] = face_edge_100_f_down_edge_100_c_edge_f_down_c;
					faces[3] = face_edge_110_f_down_edge_110_c_edge_f_down_c;
					mesh.addBody(faces);
					faces[0] = face_edge_010_110_edge_010_f_down_edge_110_f_down;
					faces[1] = face_edge_010_110_edge_010_c_edge_110_c;
					faces[2] = face_edge_010_f_down_edge_010_c_edge_f_down_c;
					faces[3] = face_edge_110_f_down_edge_110_c_edge_f_down_c;
					mesh.addBody(faces);
					faces[0] = face_edge_000_010_edge_000_f_down_edge_010_f_down;
					faces[1] = face_edge_000_010_edge_000_c_edge_010_c;
					faces[2] = face_edge_000_f_down_edge_000_c_edge_f_down_c;
					faces[3] = face_edge_010_f_down_edge_010_c_edge_f_down_c;
					mesh.addBody(faces);
					faces[0] = face_edge_000_010_edge_010_f_west_edge_000_f_west;
					faces[1] = face_edge_000_010_edge_000_c_edge_010_c;
					faces[2] = face_edge_000_f_west_edge_000_c_edge_f_west_c;
					faces[3] = face_edge_010_f_west_edge_010_c_edge_f_west_c;
					mesh.addBody(faces);
					faces[0] = face_edge_010_011_edge_010_f_west_edge_011_f_west;
					faces[1] = face_edge_010_011_edge_010_c_edge_011_c;
					faces[2] = face_edge_010_f_west_edge_010_c_edge_f_west_c;
					faces[3] = face_edge_011_f_west_edge_011_c_edge_f_west_c;
					mesh.addBody(faces);
					faces[0] = face_edge_001_011_edge_001_f_west_edge_011_f_west;
					faces[1] = face_edge_001_011_edge_001_c_edge_011_c;
					faces[2] = face_edge_001_f_west_edge_001_c_edge_f_west_c;
					faces[3] = face_edge_011_f_west_edge_011_c_edge_f_west_c;
					mesh.addBody(faces);
					faces[0] = face_edge_000_001_edge_000_f_west_edge_001_f_west;
					faces[1] = face_edge_000_001_edge_000_c_edge_001_c;
					faces[2] = face_edge_000_f_west_edge_000_c_edge_f_west_c;
					faces[3] = face_edge_001_f_west_edge_001_c_edge_f_west_c;
					mesh.addBody(faces);
					faces[0] = face_edge_000_100_edge_100_f_south_edge_000_f_south;
					faces[1] = face_edge_000_100_edge_000_c_edge_100_c;
					faces[2] = face_edge_000_f_south_edge_000_c_edge_f_south_c;
					faces[3] = face_edge_100_f_south_edge_100_c_edge_f_south_c;
					mesh.addBody(faces);
					faces[0] = face_edge_100_101_edge_100_f_south_edge_101_f_south;
					faces[1] = face_edge_100_101_edge_100_c_edge_101_c;
					faces[2] = face_edge_100_f_south_edge_100_c_edge_f_south_c;
					faces[3] = face_edge_101_f_south_edge_101_c_edge_f_south_c;
					mesh.addBody(faces);
					faces[0] = face_edge_001_101_edge_001_f_south_edge_101_f_south;
					faces[1] = face_edge_001_101_edge_001_c_edge_101_c;
					faces[2] = face_edge_001_f_south_edge_001_c_edge_f_south_c;
					faces[3] = face_edge_101_f_south_edge_101_c_edge_f_south_c;
					mesh.addBody(faces);
					faces[0] = face_edge_000_001_edge_000_f_south_edge_001_f_south;
					faces[1] = face_edge_000_001_edge_000_c_edge_001_c;
					faces[2] = face_edge_000_f_south_edge_000_c_edge_f_south_c;
					faces[3] = face_edge_001_f_south_edge_001_c_edge_f_south_c;
					mesh.addBody(faces);
					faces[0] = face_edge_100_110_edge_110_f_east_edge_100_f_east;
					faces[1] = face_edge_100_110_edge_100_c_edge_110_c;
					faces[2] = face_edge_100_f_east_edge_100_c_edge_f_east_c;
					faces[3] = face_edge_110_f_east_edge_110_c_edge_f_east_c;
					mesh.addBody(faces);
					faces[0] = face_edge_110_111_edge_110_f_east_edge_111_f_east;
					faces[1] = face_edge_110_111_edge_110_c_edge_111_c;
					faces[2] = face_edge_110_f_east_edge_110_c_edge_f_east_c;
					faces[3] = face_edge_111_f_east_edge_111_c_edge_f_east_c;
					mesh.addBody(faces);
					faces[0] = face_edge_101_111_edge_101_f_east_edge_111_f_east;
					faces[1] = face_edge_101_111_edge_101_c_edge_111_c;
					faces[2] = face_edge_101_f_east_edge_101_c_edge_f_east_c;
					faces[3] = face_edge_111_f_east_edge_111_c_edge_f_east_c;
					mesh.addBody(faces);
					faces[0] = face_edge_100_101_edge_100_f_east_edge_101_f_east;
					faces[1] = face_edge_100_101_edge_100_c_edge_101_c;
					faces[2] = face_edge_100_f_east_edge_100_c_edge_f_east_c;
					faces[3] = face_edge_101_f_east_edge_101_c_edge_f_east_c;
					mesh.addBody(faces);
					faces[0] = face_edge_010_110_edge_110_f_north_edge_010_f_north;
					faces[1] = face_edge_010_110_edge_010_c_edge_110_c;
					faces[2] = face_edge_010_f_north_edge_010_c_edge_f_north_c;
					faces[3] = face_edge_110_f_north_edge_110_c_edge_f_north_c;
					mesh.addBody(faces);
					faces[0] = face_edge_110_111_edge_110_f_north_edge_111_f_north;
					faces[1] = face_edge_110_111_edge_110_c_edge_111_c;
					faces[2] = face_edge_110_f_north_edge_110_c_edge_f_north_c;
					faces[3] = face_edge_111_f_north_edge_111_c_edge_f_north_c;
					mesh.addBody(faces);
					faces[0] = face_edge_011_111_edge_011_f_north_edge_111_f_north;
					faces[1] = face_edge_011_111_edge_011_c_edge_111_c;
					faces[2] = face_edge_011_f_north_edge_011_c_edge_f_north_c;
					faces[3] = face_edge_111_f_north_edge_111_c_edge_f_north_c;
					mesh.addBody(faces);
					faces[0] = face_edge_010_011_edge_010_f_north_edge_011_f_north;
					faces[1] = face_edge_010_011_edge_010_c_edge_011_c;
					faces[2] = face_edge_010_f_north_edge_010_c_edge_f_north_c;
					faces[3] = face_edge_011_f_north_edge_011_c_edge_f_north_c;
					mesh.addBody(faces);
					faces[0] = face_edge_001_101_edge_101_f_up_edge_001_f_up;
					faces[1] = face_edge_001_101_edge_001_c_edge_101_c;
					faces[2] = face_edge_001_f_up_edge_001_c_edge_f_up_c;
					faces[3] = face_edge_101_f_up_edge_101_c_edge_f_up_c;
					mesh.addBody(faces);
					faces[0] = face_edge_101_111_edge_101_f_up_edge_111_f_up;
					faces[1] = face_edge_101_111_edge_101_c_edge_111_c;
					faces[2] = face_edge_101_f_up_edge_101_c_edge_f_up_c;
					faces[3] = face_edge_111_f_up_edge_111_c_edge_f_up_c;
					mesh.addBody(faces);
					faces[0] = face_edge_011_111_edge_011_f_up_edge_111_f_up;
					faces[1] = face_edge_011_111_edge_011_c_edge_111_c;
					faces[2] = face_edge_011_f_up_edge_011_c_edge_f_up_c;
					faces[3] = face_edge_111_f_up_edge_111_c_edge_f_up_c;
					mesh.addBody(faces);
					faces[0] = face_edge_001_011_edge_001_f_up_edge_011_f_up;
					faces[1] = face_edge_001_011_edge_001_c_edge_011_c;
					faces[2] = face_edge_001_f_up_edge_001_c_edge_f_up_c;
					faces[3] = face_edge_011_f_up_edge_011_c_edge_f_up_c;
					mesh.addBody(faces);
				}
			}
		}
	}

	//create tetrahedral mesh with the requested number of parallellepipeds in x-, y-, and t-directions, each divided into six Bcc tetrahedra
	void createSpacetimeMeshBccParallellepiped(Mesh& mesh, double h, uint xsteps, uint ysteps, uint tsteps) {
		{
			const double sq2 = std::sqrt(2);
			const double dx = h * 1.5;
			const double dy = h * 2 * sq2 / 3;
			const double dt = h * sq2;
			const double diff_x_tstep = h * (-1.0 / 6.0);
			const double diff_y_tstep = h * sq2 / 3;
			const double diff_x_ystep = h * 7.0 / 6.0;
			for (uint i = 0; i < tsteps + 1; ++i) {
				for (uint j = 0; j < ysteps + 1; ++j) {
					double diff_y = i * diff_y_tstep;
					for (uint k = 0; k < xsteps + 1; ++k) {
						double diff_x = j * diff_x_ystep + i * diff_x_tstep;
						mesh.addNode(Vector4(diff_x + k * dx, diff_y + j * dy, i * dt, 0.0));
					}
				}
			}
			for (uint i = 0; i < tsteps; ++i) {
				for (uint j = 0; j < ysteps; ++j) {
					for (uint k = 0; k < xsteps; ++k) {
						uint node010 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
						uint node000 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
						uint node110 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
						uint node100 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
						uint node011 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
						uint node001 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
						uint node111 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
						uint node101 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
						uint edge_000_100 = mesh.addEdge(node000, node100);
						uint edge_000_010 = mesh.addEdge(node000, node010);
						uint edge_000_001 = mesh.addEdge(node000, node001);
						uint edge_000_110 = mesh.addEdge(node000, node110);
						uint edge_001_101 = mesh.addEdge(node001, node101);
						uint edge_001_011 = mesh.addEdge(node001, node011);
						uint edge_001_100 = mesh.addEdge(node001, node100);
						uint edge_001_010 = mesh.addEdge(node001, node010);
						uint edge_001_111 = mesh.addEdge(node001, node111);
						uint edge_001_110 = mesh.addEdge(node001, node110);
						uint edge_100_101 = mesh.addEdge(node100, node101);
						uint edge_101_111 = mesh.addEdge(node101, node111);
						uint edge_010_011 = mesh.addEdge(node010, node011);
						uint edge_100_110 = mesh.addEdge(node100, node110);
						uint edge_110_111 = mesh.addEdge(node110, node111);
						uint edge_010_110 = mesh.addEdge(node010, node110);
						uint edge_011_110 = mesh.addEdge(node011, node110);
						uint edge_011_111 = mesh.addEdge(node011, node111);
						uint edge_101_110 = mesh.addEdge(node101, node110);
						Buffer<uint> edges(3);
						edges[0] = edge_000_100;
						edges[1] = edge_100_110;
						edges[2] = edge_000_110;
						uint face_000_100_110 = mesh.addFace(edges);
						edges[0] = edge_000_010;
						edges[1] = edge_000_110;
						edges[2] = edge_010_110;
						uint face_000_110_010 = mesh.addFace(edges);
						edges[0] = edge_000_001;
						edges[1] = edge_001_010;
						edges[2] = edge_000_010;
						uint face_000_001_010 = mesh.addFace(edges);
						edges[0] = edge_010_011;
						edges[1] = edge_001_011;
						edges[2] = edge_001_010;
						uint face_010_011_001 = mesh.addFace(edges);
						edges[0] = edge_000_100;
						edges[1] = edge_001_100;
						edges[2] = edge_000_001;
						uint face_000_100_001 = mesh.addFace(edges);
						edges[0] = edge_100_101;
						edges[1] = edge_001_101;
						edges[2] = edge_001_100;
						uint face_100_101_001 = mesh.addFace(edges);
						edges[0] = edge_100_101;
						edges[1] = edge_101_110;
						edges[2] = edge_100_110;
						uint face_100_110_101 = mesh.addFace(edges);
						edges[0] = edge_110_111;
						edges[1] = edge_101_111;
						edges[2] = edge_101_110;
						uint face_110_111_101 = mesh.addFace(edges);
						edges[0] = edge_010_110;
						edges[1] = edge_011_110;
						edges[2] = edge_010_011;
						uint face_010_110_011 = mesh.addFace(edges);
						edges[0] = edge_110_111;
						edges[1] = edge_011_111;
						edges[2] = edge_011_110;
						uint face_110_111_011 = mesh.addFace(edges);
						edges[0] = edge_001_101;
						edges[1] = edge_101_111;
						edges[2] = edge_001_111;
						uint face_001_101_111 = mesh.addFace(edges);
						edges[0] = edge_001_011;
						edges[1] = edge_001_111;
						edges[2] = edge_011_111;
						uint face_001_111_011 = mesh.addFace(edges);
						edges[0] = edge_001_100;
						edges[1] = edge_100_110;
						edges[2] = edge_001_110;
						uint face_001_100_110 = mesh.addFace(edges);
						edges[0] = edge_001_110;
						edges[1] = edge_011_110;
						edges[2] = edge_001_011;
						uint face_001_110_011 = mesh.addFace(edges);
						edges[0] = edge_001_101;
						edges[1] = edge_101_110;
						edges[2] = edge_001_110;
						uint face_001_101_110 = mesh.addFace(edges);
						edges[0] = edge_010_110;
						edges[1] = edge_001_110;
						edges[2] = edge_001_010;
						uint face_010_110_001 = mesh.addFace(edges);
						edges[0] = edge_000_110;
						edges[1] = edge_001_110;
						edges[2] = edge_000_001;
						uint face_000_110_001 = mesh.addFace(edges);
						edges[0] = edge_001_110;
						edges[1] = edge_110_111;
						edges[2] = edge_001_111;
						uint face_001_110_111 = mesh.addFace(edges);
						Buffer<uint> faces(4);
						faces[0] = face_000_100_110;
						faces[1] = face_000_100_001;
						faces[2] = face_000_110_001;
						faces[3] = face_001_100_110;
						mesh.addBody(faces);
						faces[0] = face_100_101_001;
						faces[1] = face_100_110_101;
						faces[2] = face_001_101_110;
						faces[3] = face_001_100_110;
						mesh.addBody(faces);
						faces[0] = face_110_111_101;
						faces[1] = face_001_101_110;
						faces[2] = face_001_101_111;
						faces[3] = face_001_110_111;
						mesh.addBody(faces);
						faces[0] = face_000_110_010;
						faces[1] = face_000_001_010;
						faces[2] = face_000_110_001;
						faces[3] = face_010_110_001;
						mesh.addBody(faces);
						faces[0] = face_010_011_001;
						faces[1] = face_010_110_011;
						faces[2] = face_001_110_011;
						faces[3] = face_010_110_001;
						mesh.addBody(faces);
						faces[0] = face_001_111_011;
						faces[1] = face_110_111_011;
						faces[2] = face_001_110_011;
						faces[3] = face_001_110_111;
						mesh.addBody(faces);
					}
				}
			}
		}
	}

	/*
	Create Cartesian mesh with the requested number of boxes in x- and t-directions.
	*/
	void createCartesianMesh(Mesh& mesh, double Lx, double T, uint xsteps, uint tsteps) {
		const double dx = Lx / xsteps;
		const double dt = T / tsteps;
		for (uint i = 0; i < tsteps + 1; ++i) {
			for (uint k = 0; k < xsteps + 1; ++k) {
				mesh.addNode(Vector4(k * dx, i * dt, 0.0, 0.0));
			}
		}
		for (uint i = 0; i < tsteps; ++i) {
			for (uint k = 0; k < xsteps; ++k) {
				uint node00 = i * (xsteps + 1) + k;
				uint node10 = i * (xsteps + 1) + k + 1;
				uint node01 = (i + 1) * (xsteps + 1) + k;
				uint node11 = (i + 1) * (xsteps + 1) + k + 1;
				Buffer<uint> edges(4);
				edges[0] = mesh.addEdge(node00, node10);
				edges[1] = mesh.addEdge(node10, node11);
				edges[2] = mesh.addEdge(node01, node11);
				edges[3] = mesh.addEdge(node00, node01);
				mesh.addFace(edges);				
			}
		}
	}

	/*
	Create Cartesian mesh with the requested number of boxes in x-, y-, and t-directions.
	*/
	void createCartesianMesh(Mesh& mesh, double Lx, double Ly, double T, uint xsteps, uint ysteps, uint tsteps) {
		const double dx = Lx / xsteps;
		const double dy = Ly / ysteps;
		const double dt = T / tsteps;
		for (uint i = 0; i < tsteps + 1; ++i) {
			for (uint j = 0; j < ysteps + 1; ++j) {
				for (uint k = 0; k < xsteps + 1; ++k) {
					mesh.addNode(Vector4(k * dx, j * dy, i * dt, 0.0));
				}
			}
		}
		for (uint i = 0; i < tsteps; ++i) {
			for (uint j = 0; j < ysteps; ++j) {
				for (uint k = 0; k < xsteps; ++k) {
					uint node000 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node100 = i * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node010 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
					uint node110 = i * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
					uint node001 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k;
					uint node101 = (i + 1) * (ysteps + 1) * (xsteps + 1) + j * (xsteps + 1) + k + 1;
					uint node011 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k;
					uint node111 = (i + 1) * (ysteps + 1) * (xsteps + 1) + (j + 1) * (xsteps + 1) + k + 1;
					uint edge_000_100 = mesh.addEdge(node000, node100);
					uint edge_100_110 = mesh.addEdge(node100, node110);
					uint edge_010_110 = mesh.addEdge(node010, node110);
					uint edge_000_010 = mesh.addEdge(node000, node010);
					uint edge_000_001 = mesh.addEdge(node000, node001);
					uint edge_100_101 = mesh.addEdge(node100, node101);
					uint edge_110_111 = mesh.addEdge(node110, node111);
					uint edge_010_011 = mesh.addEdge(node010, node011);
					uint edge_001_101 = mesh.addEdge(node001, node101);
					uint edge_101_111 = mesh.addEdge(node101, node111);
					uint edge_011_111 = mesh.addEdge(node011, node111);
					uint edge_001_011 = mesh.addEdge(node001, node011);
					Buffer<uint> edges(4);
					edges[0] = edge_000_100;
					edges[1] = edge_100_110;
					edges[2] = edge_010_110;
					edges[3] = edge_000_010;
					uint face_down = mesh.addFace(edges);
					edges[0] = edge_000_100;
					edges[1] = edge_100_101;
					edges[2] = edge_001_101;
					edges[3] = edge_000_001;
					uint face_south = mesh.addFace(edges);
					edges[0] = edge_100_110;
					edges[1] = edge_110_111;
					edges[2] = edge_101_111;
					edges[3] = edge_100_101;
					uint face_east = mesh.addFace(edges);
					edges[0] = edge_010_110;
					edges[1] = edge_110_111;
					edges[2] = edge_011_111;
					edges[3] = edge_010_011;
					uint face_north = mesh.addFace(edges);
					edges[0] = edge_000_010;
					edges[1] = edge_010_011;
					edges[2] = edge_001_011;
					edges[3] = edge_000_001;
					uint face_west = mesh.addFace(edges);
					edges[0] = edge_001_101;
					edges[1] = edge_101_111;
					edges[2] = edge_011_111;
					edges[3] = edge_001_011;
					uint face_up = mesh.addFace(edges);
					Buffer<uint> faces(6);
					faces[0] = face_down;
					faces[1] = face_south;
					faces[2] = face_east;
					faces[3] = face_north;
					faces[4] = face_west;
					faces[5] = face_up;
					mesh.addBody(faces);
				}
			}
		}
	}

	//checks whether vec points in x-direction
	bool isXDir(const Vector3& vec) {
		return std::abs(vec.x) > 1e-12 && std::abs(vec.y) < 1e-12 && std::abs(vec.z) < 1e-12;
	}

	//checks whether vec points in y-direction
	bool isYDir(const Vector3& vec) {
		return std::abs(vec.y) > 1e-12 && std::abs(vec.x) < 1e-12 && std::abs(vec.z) < 1e-12;
	}

	//checks whether vec points in z-direction
	bool isZDir(const Vector3& vec) {
		return std::abs(vec.z) > 1e-12 && std::abs(vec.x) < 1e-12 && std::abs(vec.y) < 1e-12;
	}

	//checks whether vec points in x-, y-, or z-directions
	bool isCartesian(const Vector3& vec) {
		return (isXDir(vec) || isYDir(vec) || isZDir(vec));
	}

	//defines how vectors are printed to the console (std::cout << vec)
	std::ostream& operator<< (std::ostream& out, const Vector3& vec) {
		out << '(';
		if (std::abs(vec.x) < 1e-12)
			out << '0';
		else
			out << vec.x;
		out << ", ";
		if (std::abs(vec.y) < 1e-12)
			out << '0';
		else
			out << vec.y;
		out << ", ";
		if (std::abs(vec.z) < 1e-12)
			out << '0';
		else
			out << vec.z;
		out << ')';
		return out;
	}
}