#include "SmallCubePartition3D.hpp"
#include "WhitneyForm.hpp"
#include "NumericalIntegration.hpp"
#include "LinearAlgebraFunctions.hpp"

using namespace gfd;

SmallCubePartition3D::SmallCubePartition3D(uint order) : order{ order } {
	initialiseSmallCells();
	formMatrices();
}

void SmallCubePartition3D::initialiseSmallCells() {
	//small edges in edges
	smallEdgesInEdges.reserve(order);
	for (uint xi = 0; xi < order; ++xi) {
		mi_t mi{ {xi} };
		smallEdgesInEdges.push_back({  mi, 0 });
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
	//small edges in bodies
	smallEdgesInBodies.reserve(3 * order * (order - 1) * (order - 1));
	for (uint zi = 1; zi < order; ++zi) {
		for (uint yi = 1; yi < order; ++yi) {
			for (uint xi = 0; xi < order; ++xi) {
				smallEdgesInBodies.push_back({ {xi, yi, zi}, 0 });
			}
		}
	}
	for (uint zi = 1; zi < order; ++zi) {
		for (uint xi = 1; xi < order; ++xi) {
			for (uint yi = 0; yi < order; ++yi) {
				smallEdgesInBodies.push_back({ {xi, yi, zi}, 1 });
			}
		}
	}
	for (uint yi = 1; yi < order; ++yi) {
		for (uint xi = 1; xi < order; ++xi) {
			for (uint zi = 0; zi < order; ++zi) {
				smallEdgesInBodies.push_back({ {xi, yi, zi}, 2 });
			}
		}
	}
}

void SmallCubePartition3D::refineMesh(const BuilderMesh& mesh_old, BuilderMesh& mesh) {
	//SmallCubePartition3D is now associated with the old mesh (mesh_old) and its refinement (mesh)
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
		Vector3 p0 = mesh.getNodePosition3(nodes[0]);
		Vector3 p1 = mesh.getNodePosition3(nodes[1]);
		//small nodes on edge i
		for (uint xi = 1; xi < order; ++xi) {
			double weight_x = (double)xi / (double)order;
			mesh.addNode(Vector4((1 - weight_x) * p0 + weight_x * p1, 0.0));
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
		Vector3 p0 = mesh.getNodePosition3(nodes[0]);
		Vector3 p1 = mesh.getNodePosition3(nodes[1]);
		Vector3 p2 = mesh.getNodePosition3(nodes[2]);
		Vector3 p3 = mesh.getNodePosition3(nodes[3]);
		Vector3 xdiff = p1 - p0;
		Vector3 ydiff = p3 - p0;
		//small nodes in face i
		for (uint yi = 1; yi < order; ++yi) {
			double weight_y = (double)yi / (double)order;
			for (uint xi = 1; xi < order; ++xi) {
				double weight_x = (double)xi / (double)order;
				mesh.addNode(Vector4(p0 + weight_x * xdiff + weight_y * ydiff, 0.0));
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
				Vector3 pos0 = p0 + weight_x_start * xdiff + weight_y_start * ydiff;
				Vector3 pos1 = p0 + weight_x_end * xdiff + weight_y_start * ydiff;
				Vector3 pos2 = p0 + weight_x_end * xdiff + weight_y_end * ydiff;
				Vector3 pos3 = p0 + weight_x_start * xdiff + weight_y_end * ydiff;
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

	//refine bodies
	for (uint i = 0; i < mesh_old.getBodySize(); ++i) {
		const Buffer<uint> nodes = getCubeNodes(i, mesh_old);
		Vector3 p0 = mesh.getNodePosition3(nodes[0]);
		Vector3 p1 = mesh.getNodePosition3(nodes[1]);
		Vector3 p2 = mesh.getNodePosition3(nodes[2]);
		Vector3 p3 = mesh.getNodePosition3(nodes[3]);
		Vector3 p4 = mesh.getNodePosition3(nodes[4]);
		Vector3 p5 = mesh.getNodePosition3(nodes[5]);
		Vector3 p6 = mesh.getNodePosition3(nodes[6]);
		Vector3 p7 = mesh.getNodePosition3(nodes[7]);
		Vector3 xdiff = p1 - p0;
		Vector3 ydiff = p3 - p0;
		Vector3 zdiff = p4 - p0;
		//small nodes in body i
		for (uint zi = 1; zi < order; ++zi) {
			double weight_z = (double)zi / (double)order;
			for (uint yi = 1; yi < order; ++yi) {
				double weight_y = (double)yi / (double)order;
				for (uint xi = 1; xi < order; ++xi) {
					double weight_x = (double)xi / (double)order;
					mesh.addNode(Vector4(p0 + weight_x * xdiff + weight_y * ydiff + weight_z * zdiff, 0.0));
				}
			}
		}
		//small edges in body i
		uint face_yz0 = findQuadrilateral(nodes[0], nodes[3], nodes[4], nodes[7], mesh_old);
		uint face_yz1 = findQuadrilateral(nodes[1], nodes[2], nodes[5], nodes[6], mesh_old);
		uint face_xz0 = findQuadrilateral(nodes[0], nodes[1], nodes[4], nodes[5], mesh_old);
		uint face_xz1 = findQuadrilateral(nodes[2], nodes[3], nodes[6], nodes[7], mesh_old);
		uint face_xy0 = findQuadrilateral(nodes[0], nodes[1], nodes[2], nodes[3], mesh_old);
		uint face_xy1 = findQuadrilateral(nodes[4], nodes[5], nodes[6], nodes[7], mesh_old);
		for (uint zi = 1; zi < order; ++zi) {
			double weight_z = (double)zi / (double)order;
			for (uint yi = 1; yi < order; ++yi) {
				double weight_y = (double)yi / (double)order;
				uint prev = getSmallNodeIndex(p0 + weight_y * ydiff + weight_z * zdiff, face_yz0, 2);
				for (uint xi = 1; xi < order; ++xi) {
					double weight_x = (double)xi / (double)order;
					uint curr = getSmallNodeIndex(p0 + weight_x * xdiff + weight_y * ydiff + weight_z * zdiff, i, 3);
					mesh.addEdge(prev, curr);
					prev = curr;
				}
				mesh.addEdge(prev, getSmallNodeIndex(p1 + weight_y * ydiff + weight_z * zdiff, face_yz1, 2));
			}
		}
		for (uint zi = 1; zi < order; ++zi) {
			double weight_z = (double)zi / (double)order;
			for (uint xi = 1; xi < order; ++xi) {
				double weight_x = (double)xi / (double)order;
				uint prev = getSmallNodeIndex(p0 + weight_x * xdiff + weight_z * zdiff, face_xz0, 2);
				for (uint yi = 1; yi < order; ++yi) {
					double weight_y = (double)yi / (double)order;
					uint curr = getSmallNodeIndex(p0 + weight_x * xdiff + weight_y * ydiff + weight_z * zdiff, i, 3);
					mesh.addEdge(prev, curr);
					prev = curr;
				}
				mesh.addEdge(prev, getSmallNodeIndex(p3 + weight_x * xdiff + weight_z * zdiff, face_xz1, 2));
			}
		}
		for (uint yi = 1; yi < order; ++yi) {
			double weight_y = (double)yi / (double)order;
			for (uint xi = 1; xi < order; ++xi) {
				double weight_x = (double)xi / (double)order;
				uint prev = getSmallNodeIndex(p0 + weight_x * xdiff + weight_y * ydiff, face_xy0, 2);
				for (uint zi = 1; zi < order; ++zi) {
					double weight_z = (double)zi / (double)order;
					uint curr = getSmallNodeIndex(p0 + weight_x * xdiff + weight_y * ydiff + weight_z * zdiff, i, 3);
					mesh.addEdge(prev, curr);
					prev = curr;
				}
				mesh.addEdge(prev, getSmallNodeIndex(p4 + weight_x * xdiff + weight_y * ydiff, face_xy1, 2));
			}
		}
		//small faces in body i
		uint edge01 = mesh_old.findEdge(nodes[0], nodes[1]);
		uint edge12 = mesh_old.findEdge(nodes[1], nodes[2]);
		uint edge23 = mesh_old.findEdge(nodes[2], nodes[3]);
		uint edge03 = mesh_old.findEdge(nodes[0], nodes[3]);
		uint edge04 = mesh_old.findEdge(nodes[0], nodes[4]);
		uint edge15 = mesh_old.findEdge(nodes[1], nodes[5]);
		uint edge26 = mesh_old.findEdge(nodes[2], nodes[6]);
		uint edge37 = mesh_old.findEdge(nodes[3], nodes[7]);
		uint edge45 = mesh_old.findEdge(nodes[4], nodes[5]);
		uint edge56 = mesh_old.findEdge(nodes[5], nodes[6]);
		uint edge67 = mesh_old.findEdge(nodes[6], nodes[7]);
		uint edge47 = mesh_old.findEdge(nodes[4], nodes[7]);
		for (uint zi = 1; zi < order; ++zi) {
			double weight_z = (double)zi / (double)order;
			for (uint yi = 0; yi < order; ++yi) {
				double weight_y_start = (double)yi / (double)order;
				double weight_y_end = (double)(yi + 1) / (double)order;
				for (uint xi = 0; xi < order; ++xi) {
					double weight_x_start = (double)xi / (double)order;
					double weight_x_end = (double)(xi + 1) / (double)order;
					uint node0, node1, node2, node3;
					Vector3 pos0 = p0 + weight_x_start * xdiff + weight_y_start * ydiff + weight_z * zdiff;
					Vector3 pos1 = p0 + weight_x_end * xdiff + weight_y_start * ydiff + weight_z * zdiff;
					Vector3 pos2 = p0 + weight_x_end * xdiff + weight_y_end * ydiff + weight_z * zdiff;
					Vector3 pos3 = p0 + weight_x_start * xdiff + weight_y_end * ydiff + weight_z * zdiff;
					if (xi == 0 && yi == 0) {
						node0 = getSmallNodeIndex(pos0, edge04, 1);
					}
					else if (xi == 0) {
						node0 = getSmallNodeIndex(pos0, face_yz0, 2);
					}
					else if (yi == 0) {
						node0 = getSmallNodeIndex(pos0, face_xz0, 2);
					}
					else {
						node0 = getSmallNodeIndex(pos0, i, 3);
					}
					if (xi == order - 1 && yi == 0) {
						node1 = getSmallNodeIndex(pos1, edge15, 1);
					}
					else if (xi == order - 1) {
						node1 = getSmallNodeIndex(pos1, face_yz1, 2);
					}
					else if (yi == 0) {
						node1 = getSmallNodeIndex(pos1, face_xz0, 2);
					}
					else {
						node1 = getSmallNodeIndex(pos1, i, 3);
					}
					if (xi == order - 1 && yi == order - 1) {
						node2 = getSmallNodeIndex(pos2, edge26, 1);
					}
					else if (xi == order - 1) {
						node2 = getSmallNodeIndex(pos2, face_yz1, 2);
					}
					else if (yi == order - 1) {
						node2 = getSmallNodeIndex(pos2, face_xz1, 2);
					}
					else {
						node2 = getSmallNodeIndex(pos2, i, 3);
					}
					if (xi == 0 && yi == order - 1) {
						node3 = getSmallNodeIndex(pos3, edge37, 1);
					}
					else if (xi == 0) {
						node3 = getSmallNodeIndex(pos3, face_yz0, 2);
					}
					else if (yi == order - 1) {
						node3 = getSmallNodeIndex(pos3, face_xz1, 2);
					}
					else {
						node3 = getSmallNodeIndex(pos3, i, 3);
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
		for (uint yi = 1; yi < order; ++yi) {
			double weight_y = (double)yi / (double)order;
			for (uint zi = 0; zi < order; ++zi) {
				double weight_z_start = (double)zi / (double)order;
				double weight_z_end = (double)(zi + 1) / (double)order;
				for (uint xi = 0; xi < order; ++xi) {
					double weight_x_start = (double)xi / (double)order;
					double weight_x_end = (double)(xi + 1) / (double)order;
					uint node0, node1, node2, node3;
					Vector3 pos0 = p0 + weight_x_start * xdiff + weight_z_start * zdiff + weight_y * ydiff;
					Vector3 pos1 = p0 + weight_x_end * xdiff + weight_z_start * zdiff + weight_y * ydiff;
					Vector3 pos2 = p0 + weight_x_end * xdiff + weight_z_end * zdiff + weight_y * ydiff;
					Vector3 pos3 = p0 + weight_x_start * xdiff + weight_z_end * zdiff + weight_y * ydiff;
					if (xi == 0 && zi == 0) {
						node0 = getSmallNodeIndex(pos0, edge03, 1);
					}
					else if (xi == 0) {
						node0 = getSmallNodeIndex(pos0, face_yz0, 2);
					}
					else if (zi == 0) {
						node0 = getSmallNodeIndex(pos0, face_xy0, 2);
					}
					else {
						node0 = getSmallNodeIndex(pos0, i, 3);
					}
					if (xi == order - 1 && zi == 0) {
						node1 = getSmallNodeIndex(pos1, edge12, 1);
					}
					else if (xi == order - 1) {
						node1 = getSmallNodeIndex(pos1, face_yz1, 2);
					}
					else if (zi == 0) {
						node1 = getSmallNodeIndex(pos1, face_xy0, 2);
					}
					else {
						node1 = getSmallNodeIndex(pos1, i, 3);
					}
					if (xi == order - 1 && zi == order - 1) {
						node2 = getSmallNodeIndex(pos2, edge56, 1);
					}
					else if (xi == order - 1) {
						node2 = getSmallNodeIndex(pos2, face_yz1, 2);
					}
					else if (zi == order - 1) {
						node2 = getSmallNodeIndex(pos2, face_xy1, 2);
					}
					else {
						node2 = getSmallNodeIndex(pos2, i, 3);
					}
					if (xi == 0 && zi == order - 1) {
						node3 = getSmallNodeIndex(pos3, edge47, 1);
					}
					else if (xi == 0) {
						node3 = getSmallNodeIndex(pos3, face_yz0, 2);
					}
					else if (zi == order - 1) {
						node3 = getSmallNodeIndex(pos3, face_xy1, 2);
					}
					else {
						node3 = getSmallNodeIndex(pos3, i, 3);
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
		for (uint xi = 1; xi < order; ++xi) {
			double weight_x = (double)xi / (double)order;
			for (uint zi = 0; zi < order; ++zi) {
				double weight_z_start = (double)zi / (double)order;
				double weight_z_end = (double)(zi + 1) / (double)order;
				for (uint yi = 0; yi < order; ++yi) {
					double weight_y_start = (double)yi / (double)order;
					double weight_y_end = (double)(yi + 1) / (double)order;
					uint node0, node1, node2, node3;
					Vector3 pos0 = p0 + weight_y_start * ydiff + weight_z_start * zdiff + weight_x * xdiff;
					Vector3 pos1 = p0 + weight_y_end * ydiff + weight_z_start * zdiff + weight_x * xdiff;
					Vector3 pos2 = p0 + weight_y_end * ydiff + weight_z_end * zdiff + weight_x * xdiff;
					Vector3 pos3 = p0 + weight_y_start * ydiff + weight_z_end * zdiff + weight_x * xdiff;
					if (yi == 0 && zi == 0) {
						node0 = getSmallNodeIndex(pos0, edge01, 1);
					}
					else if (yi == 0) {
						node0 = getSmallNodeIndex(pos0, face_xz0, 2);
					}
					else if (zi == 0) {
						node0 = getSmallNodeIndex(pos0, face_xy0, 2);
					}
					else {
						node0 = getSmallNodeIndex(pos0, i, 3);
					}
					if (yi == order - 1 && zi == 0) {
						node1 = getSmallNodeIndex(pos1, edge23, 1);
					}
					else if (yi == order - 1) {
						node1 = getSmallNodeIndex(pos1, face_xz1, 2);
					}
					else if (zi == 0) {
						node1 = getSmallNodeIndex(pos1, face_xy0, 2);
					}
					else {
						node1 = getSmallNodeIndex(pos1, i, 3);
					}
					if (yi == order - 1 && zi == order - 1) {
						node2 = getSmallNodeIndex(pos2, edge67, 1);
					}
					else if (yi == order - 1) {
						node2 = getSmallNodeIndex(pos2, face_xz1, 2);
					}
					else if (zi == order - 1) {
						node2 = getSmallNodeIndex(pos2, face_xy1, 2);
					}
					else {
						node2 = getSmallNodeIndex(pos2, i, 3);
					}
					if (yi == 0 && zi == order - 1) {
						node3 = getSmallNodeIndex(pos3, edge45, 1);
					}
					else if (yi == 0) {
						node3 = getSmallNodeIndex(pos3, face_xz0, 2);
					}
					else if (zi == order - 1) {
						node3 = getSmallNodeIndex(pos3, face_xy1, 2);
					}
					else {
						node3 = getSmallNodeIndex(pos3, i, 3);
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
		//small bodies in body i
		for (uint zi = 0; zi < order; ++zi) {
			double weight_z_start = (double)zi / (double)order;
			double weight_z_end = (double)(zi + 1) / (double)order;
			for (uint yi = 0; yi < order; ++yi) {
				double weight_y_start = (double)yi / (double)order;
				double weight_y_end = (double)(yi + 1) / (double)order;
				for (uint xi = 0; xi < order; ++xi) {
					double weight_x_start = (double)xi / (double)order;
					double weight_x_end = (double)(xi + 1) / (double)order;
					uint node0, node1, node2, node3, node4, node5, node6, node7;
					Vector3 pos0 = p0 + weight_x_start * xdiff + weight_y_start * ydiff + weight_z_start * zdiff;
					Vector3 pos1 = p0 + weight_x_end * xdiff + weight_y_start * ydiff + weight_z_start * zdiff;
					Vector3 pos2 = p0 + weight_x_end * xdiff + weight_y_end * ydiff + weight_z_start * zdiff;
					Vector3 pos3 = p0 + weight_x_start * xdiff + weight_y_end * ydiff + weight_z_start * zdiff;
					Vector3 pos4 = p0 + weight_x_start * xdiff + weight_y_start * ydiff + weight_z_end * zdiff;
					Vector3 pos5 = p0 + weight_x_end * xdiff + weight_y_start * ydiff + weight_z_end * zdiff;
					Vector3 pos6 = p0 + weight_x_end * xdiff + weight_y_end * ydiff + weight_z_end * zdiff;
					Vector3 pos7 = p0 + weight_x_start * xdiff + weight_y_end * ydiff + weight_z_end * zdiff;
					if (xi == 0 && yi == 0 && zi == 0) {
						node0 = nodes[0];
					}
					else if (xi == 0 && yi == 0) {
						node0 = getSmallNodeIndex(pos0, edge04, 1);
					}
					else if (xi == 0 && zi == 0) {
						node0 = getSmallNodeIndex(pos0, edge03, 1);
					}
					else if (yi == 0 && zi == 0) {
						node0 = getSmallNodeIndex(pos0, edge01, 1);
					}
					else if (xi == 0) {
						node0 = getSmallNodeIndex(pos0, face_yz0, 2);
					}
					else if (yi == 0) {
						node0 = getSmallNodeIndex(pos0, face_xz0, 2);
					}
					else if (zi == 0) {
						node0 = getSmallNodeIndex(pos0, face_xy0, 2);
					}
					else {
						node0 = getSmallNodeIndex(pos0, i, 3);
					}
					if (xi == order - 1 && yi == 0 && zi == 0) {
						node1 = nodes[1];
					}
					else if (xi == order - 1 && yi == 0) {
						node1 = getSmallNodeIndex(pos1, edge15, 1);
					}
					else if (xi == order - 1 && zi == 0) {
						node1 = getSmallNodeIndex(pos1, edge12, 1);
					}
					else if (yi == 0 && zi == 0) {
						node1 = getSmallNodeIndex(pos1, edge01, 1);
					}
					else if (xi == order - 1) {
						node1 = getSmallNodeIndex(pos1, face_yz1, 2);
					}
					else if (yi == 0) {
						node1 = getSmallNodeIndex(pos1, face_xz0, 2);
					}
					else if (zi == 0) {
						node1 = getSmallNodeIndex(pos1, face_xy0, 2);
					}
					else {
						node1 = getSmallNodeIndex(pos1, i, 3);
					}
					if (xi == order - 1 && yi == order - 1 && zi == 0) {
						node2 = nodes[2];
					}
					else if (xi == order - 1 && yi == order - 1) {
						node2 = getSmallNodeIndex(pos2, edge26, 1);
					}
					else if (xi == order - 1 && zi == 0) {
						node2 = getSmallNodeIndex(pos2, edge12, 1);
					}
					else if (yi == order - 1 && zi == 0) {
						node2 = getSmallNodeIndex(pos2, edge23, 1);
					}
					else if (xi == order - 1) {
						node2 = getSmallNodeIndex(pos2, face_yz1, 2);
					}
					else if (yi == order - 1) {
						node2 = getSmallNodeIndex(pos2, face_xz1, 2);
					}
					else if (zi == 0) {
						node2 = getSmallNodeIndex(pos2, face_xy0, 2);
					}
					else {
						node2 = getSmallNodeIndex(pos2, i, 3);
					}
					if (xi == 0 && yi == order - 1 && zi == 0) {
						node3 = nodes[3];
					}
					else if (xi == 0 && yi == order - 1) {
						node3 = getSmallNodeIndex(pos3, edge37, 1);
					}
					else if (xi == 0 && zi == 0) {
						node3 = getSmallNodeIndex(pos3, edge03, 1);
					}
					else if (yi == order - 1 && zi == 0) {
						node3 = getSmallNodeIndex(pos3, edge23, 1);
					}
					else if (xi == 0) {
						node3 = getSmallNodeIndex(pos3, face_yz0, 2);
					}
					else if (yi == order - 1) {
						node3 = getSmallNodeIndex(pos3, face_xz1, 2);
					}
					else if (zi == 0) {
						node3 = getSmallNodeIndex(pos3, face_xy0, 2);
					}
					else {
						node3 = getSmallNodeIndex(pos3, i, 3);
					}
					if (xi == 0 && yi == 0 && zi == order - 1) {
						node4 = nodes[4];
					}
					else if (xi == 0 && yi == 0) {
						node4 = getSmallNodeIndex(pos4, edge04, 1);
					}
					else if (xi == 0 && zi == order - 1) {
						node4 = getSmallNodeIndex(pos4, edge47, 1);
					}
					else if (yi == 0 && zi == order - 1) {
						node4 = getSmallNodeIndex(pos4, edge45, 1);
					}
					else if (xi == 0) {
						node4 = getSmallNodeIndex(pos4, face_yz0, 2);
					}
					else if (yi == 0) {
						node4 = getSmallNodeIndex(pos4, face_xz0, 2);
					}
					else if (zi == order - 1) {
						node4 = getSmallNodeIndex(pos4, face_xy1, 2);
					}
					else {
						node4 = getSmallNodeIndex(pos4, i, 3);
					}
					if (xi == order - 1 && yi == 0 && zi == order - 1) {
						node5 = nodes[5];
					}
					else if (xi == order - 1 && yi == 0) {
						node5 = getSmallNodeIndex(pos5, edge15, 1);
					}
					else if (xi == order - 1 && zi == order - 1) {
						node5 = getSmallNodeIndex(pos5, edge56, 1);
					}
					else if (yi == 0 && zi == order - 1) {
						node5 = getSmallNodeIndex(pos5, edge45, 1);
					}
					else if (xi == order - 1) {
						node5 = getSmallNodeIndex(pos5, face_yz1, 2);
					}
					else if (yi == 0) {
						node5 = getSmallNodeIndex(pos5, face_xz0, 2);
					}
					else if (zi == order - 1) {
						node5 = getSmallNodeIndex(pos5, face_xy1, 2);
					}
					else {
						node5 = getSmallNodeIndex(pos5, i, 3);
					}
					if (xi == order - 1 && yi == order - 1 && zi == order - 1) {
						node6 = nodes[6];
					}
					else if (xi == order - 1 && yi == order - 1) {
						node6 = getSmallNodeIndex(pos6, edge26, 1);
					}
					else if (xi == order - 1 && zi == order - 1) {
						node6 = getSmallNodeIndex(pos6, edge56, 1);
					}
					else if (yi == order - 1 && zi == order - 1) {
						node6 = getSmallNodeIndex(pos6, edge67, 1);
					}
					else if (xi == order - 1) {
						node6 = getSmallNodeIndex(pos6, face_yz1, 2);
					}
					else if (yi == order - 1) {
						node6 = getSmallNodeIndex(pos6, face_xz1, 2);
					}
					else if (zi == order - 1) {
						node6 = getSmallNodeIndex(pos6, face_xy1, 2);
					}
					else {
						node6 = getSmallNodeIndex(pos6, i, 3);
					}
					if (xi == 0 && yi == order - 1 && zi == order - 1) {
						node7 = nodes[7];
					}
					else if (xi == 0 && yi == order - 1) {
						node7 = getSmallNodeIndex(pos7, edge37, 1);
					}
					else if (xi == 0 && zi == order - 1) {
						node7 = getSmallNodeIndex(pos7, edge47, 1);
					}
					else if (yi == order - 1 && zi == order - 1) {
						node7 = getSmallNodeIndex(pos7, edge67, 1);
					}
					else if (xi == 0) {
						node7 = getSmallNodeIndex(pos7, face_yz0, 2);
					}
					else if (yi == order - 1) {
						node7 = getSmallNodeIndex(pos7, face_xz1, 2);
					}
					else if (zi == order - 1) {
						node7 = getSmallNodeIndex(pos7, face_xy1, 2);
					}
					else {
						node7 = getSmallNodeIndex(pos7, i, 3);
					}
					Buffer<uint> faces(6);
					faces[0] = findQuadrilateral(node0, node1, node2, node3, mesh);
					faces[1] = findQuadrilateral(node0, node1, node4, node5, mesh);
					faces[2] = findQuadrilateral(node1, node2, node5, node6, mesh);
					faces[3] = findQuadrilateral(node2, node3, node6, node7, mesh);
					faces[4] = findQuadrilateral(node0, node3, node4, node7, mesh);
					faces[5] = findQuadrilateral(node4, node5, node6, node7, mesh);
					mesh.addBody(faces);
				}
			}
		}
	}
}

//solves the coefficients of the (higher order) cubical 1-form (interpolant of discreteForm)
void SmallCubePartition3D::solve1FormCoefficients(const Buffer<double>& discreteForm, Buffer<VectorN>& edgeCoefficients, Buffer<VectorN>& faceCoefficients,
	Buffer<VectorN>& bodyCoefficients) const {

	if (!mesh_old_ptr)
		return; //SmallCubePartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;

	edgeCoefficients.resize(mesh_old.getEdgeSize()); //coefficients for small edges that are on big edges
	for (uint i = 0; i < edgeCoefficients.size(); ++i) {
		VectorN edgeValues(smallEdgesInEdges.size(), 0.0);
		for (uint j = 0; j < smallEdgesInEdges.size(); ++j) {
			edgeValues[j] = discreteForm[smallEdgesEdgeList(i, j)];
		}
		solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeValues, edgeCoefficients[i]);
	}

	faceCoefficients.resize(mesh_old.getFaceSize()); //coefficients for small edges that are on big faces
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

	bodyCoefficients.resize(mesh_old.getBodySize()); //coefficients for small edges in the interior of the cube
	for (uint i = 0; i < bodyCoefficients.size(); ++i) {
		const Buffer<uint> nodes = getCubeNodes(i, mesh_old);
		Buffer<uint> edges(12);
		edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
		edges[1] = mesh_old.findEdge(nodes[1], nodes[2]);
		edges[2] = mesh_old.findEdge(nodes[2], nodes[3]);
		edges[3] = mesh_old.findEdge(nodes[0], nodes[3]);
		edges[4] = mesh_old.findEdge(nodes[0], nodes[4]);
		edges[5] = mesh_old.findEdge(nodes[1], nodes[5]);
		edges[6] = mesh_old.findEdge(nodes[2], nodes[6]);
		edges[7] = mesh_old.findEdge(nodes[3], nodes[7]);
		edges[8] = mesh_old.findEdge(nodes[4], nodes[5]);
		edges[9] = mesh_old.findEdge(nodes[5], nodes[6]);
		edges[10] = mesh_old.findEdge(nodes[6], nodes[7]);
		edges[11] = mesh_old.findEdge(nodes[4], nodes[7]);
		Buffer<uint> faces(6);
		faces[0] = findQuadrilateral(nodes[0], nodes[1], nodes[2], nodes[3], mesh_old);
		faces[1] = findQuadrilateral(nodes[0], nodes[1], nodes[4], nodes[5], mesh_old);
		faces[2] = findQuadrilateral(nodes[1], nodes[2], nodes[5], nodes[6], mesh_old);
		faces[3] = findQuadrilateral(nodes[2], nodes[3], nodes[6], nodes[7], mesh_old);
		faces[4] = findQuadrilateral(nodes[0], nodes[3], nodes[4], nodes[7], mesh_old);
		faces[5] = findQuadrilateral(nodes[4], nodes[5], nodes[6], nodes[7], mesh_old);
		VectorN bodyValues(smallEdgesInBodies.size(), 0.0);
		for (uint j = 0; j < smallEdgesInBodies.size(); ++j) {
			for (uint k = 0; k < bodyEdgeValues1Forms.size(); ++k) {
				bodyValues[j] -= bodyEdgeValues1Forms[k][j].dot(edgeCoefficients[edges[k]]);
			}
			for (uint k = 0; k < bodyFaceValues1Forms.size(); ++k) {
				bodyValues[j] -= bodyFaceValues1Forms[k][j].dot(faceCoefficients[faces[k]]);
			}
			bodyValues[j] += discreteForm[smallEdgesBodyList(i, j)];
		}
		solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyValues, bodyCoefficients[i]);
	}
}

//computes the value of the (higher order) cubical 1-form (interpolant of discreteForm) at evaluationPoint when the coefficients are known
Vector3 SmallCubePartition3D::evaluate1FormWithCoefficients(const Vector3& evaluationPoint, const Buffer<VectorN>& edgeCoefficients,
	const Buffer<VectorN>& faceCoefficients, const Buffer<VectorN>& bodyCoefficients, uint element) const {
	if (!mesh_old_ptr)
		return Vector3(0.0, 0.0, 0.0); //SmallCubePartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;

	const Buffer<uint> nodes = getCubeNodes(element, mesh_old);
	Buffer<uint> edges(12);
	edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
	edges[1] = mesh_old.findEdge(nodes[1], nodes[2]);
	edges[2] = mesh_old.findEdge(nodes[2], nodes[3]);
	edges[3] = mesh_old.findEdge(nodes[0], nodes[3]);
	edges[4] = mesh_old.findEdge(nodes[0], nodes[4]);
	edges[5] = mesh_old.findEdge(nodes[1], nodes[5]);
	edges[6] = mesh_old.findEdge(nodes[2], nodes[6]);
	edges[7] = mesh_old.findEdge(nodes[3], nodes[7]);
	edges[8] = mesh_old.findEdge(nodes[4], nodes[5]);
	edges[9] = mesh_old.findEdge(nodes[5], nodes[6]);
	edges[10] = mesh_old.findEdge(nodes[6], nodes[7]);
	edges[11] = mesh_old.findEdge(nodes[4], nodes[7]);
	Buffer<uint> faces(6);
	faces[0] = findQuadrilateral(nodes[0], nodes[1], nodes[2], nodes[3], mesh_old);
	faces[1] = findQuadrilateral(nodes[0], nodes[1], nodes[4], nodes[5], mesh_old);
	faces[2] = findQuadrilateral(nodes[1], nodes[2], nodes[5], nodes[6], mesh_old);
	faces[3] = findQuadrilateral(nodes[2], nodes[3], nodes[6], nodes[7], mesh_old);
	faces[4] = findQuadrilateral(nodes[0], nodes[3], nodes[4], nodes[7], mesh_old);
	faces[5] = findQuadrilateral(nodes[4], nodes[5], nodes[6], nodes[7], mesh_old);

	Vector3 result(0.0, 0.0, 0.0);
	Vector3 p0 = mesh_old.getNodePosition3(nodes[0]);
	Vector3 p1 = mesh_old.getNodePosition3(nodes[6]);
	double lenx = p1.x - p0.x;
	double leny = p1.y - p0.y;
	double lenz = p1.z - p0.z;
	double x = (evaluationPoint.x - p0.x) / lenx;
	double y = (evaluationPoint.y - p0.y) / leny;
	double z = (evaluationPoint.z - p0.z) / lenz;
	double x_d = (p1.x - evaluationPoint.x) / lenx;
	double y_d = (p1.y - evaluationPoint.y) / leny;
	double z_d = (p1.z - evaluationPoint.z) / lenz;
	double sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge0
		sum += edgeCoefficients[edges[0]][i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(y_d, order) * std::pow(z_d, order);
	result.x += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge1
		sum += edgeCoefficients[edges[1]][i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x, order) * std::pow(z_d, order);
	result.y += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge2
		sum += edgeCoefficients[edges[2]][i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(y, order) * std::pow(z_d, order);
	result.x += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge3
		sum += edgeCoefficients[edges[3]][i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x_d, order) * std::pow(z_d, order);
	result.y += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge4
		sum += edgeCoefficients[edges[4]][i] * std::pow(z_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(z, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x_d, order) * std::pow(y_d, order);
	result.z += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge5
		sum += edgeCoefficients[edges[5]][i] * std::pow(z_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(z, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x, order) * std::pow(y_d, order);
	result.z += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge6
		sum += edgeCoefficients[edges[6]][i] * std::pow(z_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(z, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x, order) * std::pow(y, order);
	result.z += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge7
		sum += edgeCoefficients[edges[7]][i] * std::pow(z_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(z, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x_d, order) * std::pow(y, order);
	result.z += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge8
		sum += edgeCoefficients[edges[8]][i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(y_d, order) * std::pow(z, order);
	result.x += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge9
		sum += edgeCoefficients[edges[9]][i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x, order) * std::pow(z, order);
	result.y += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge10
		sum += edgeCoefficients[edges[10]][i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(y, order) * std::pow(z, order);
	result.x += sum;
	sum = 0.0;
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) { //edge11
		sum += edgeCoefficients[edges[11]][i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
	}
	sum *= std::pow(x_d, order) * std::pow(z, order);
	result.y += sum;
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) { //face0
		const SmallEdge& se = smallEdgesInFaces[i];
		if (se.face == 0) {
			result.x += faceCoefficients[faces[0]][i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]) * std::pow(z_d, order);
		}
		else {
			result.y += faceCoefficients[faces[0]][i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - 1 - se.mi[1]) * std::pow(y, se.mi[1]) * std::pow(z_d, order);
		}

	}
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) { //face1
		const SmallEdge& se = smallEdgesInFaces[i];
		if (se.face == 0) {
			result.x += faceCoefficients[faces[1]][i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(z_d, order - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(y_d, order);
		}
		else {
			result.z += faceCoefficients[faces[1]][i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(z_d, order - 1 - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(y_d, order);
		}

	}
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) { //face2
		const SmallEdge& se = smallEdgesInFaces[i];
		if (se.face == 0) {
			result.y += faceCoefficients[faces[2]][i] * std::pow(y_d, order - 1 - se.mi[0]) * std::pow(y, se.mi[0]) * std::pow(z_d, order - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(x, order);
		}
		else {
			result.z += faceCoefficients[faces[2]][i] * std::pow(y_d, order - se.mi[0]) * std::pow(y, se.mi[0]) * std::pow(z_d, order - 1 - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(x, order);
		}

	}
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) { //face3
		const SmallEdge& se = smallEdgesInFaces[i];
		if (se.face == 0) {
			result.x += faceCoefficients[faces[3]][i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(z_d, order - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(y, order);
		}
		else {
			result.z += faceCoefficients[faces[3]][i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(z_d, order - 1 - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(y, order);
		}

	}
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) { //face4
		const SmallEdge& se = smallEdgesInFaces[i];
		if (se.face == 0) {
			result.y += faceCoefficients[faces[4]][i] * std::pow(y_d, order - 1 - se.mi[0]) * std::pow(y, se.mi[0]) * std::pow(z_d, order - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(x_d, order);
		}
		else {
			result.z += faceCoefficients[faces[4]][i] * std::pow(y_d, order - se.mi[0]) * std::pow(y, se.mi[0]) * std::pow(z_d, order - 1 - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(x_d, order);
		}

	}
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) { //face5
		const SmallEdge& se = smallEdgesInFaces[i];
		if (se.face == 0) {
			result.x += faceCoefficients[faces[5]][i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]) * std::pow(z, order);
		}
		else {
			result.y += faceCoefficients[faces[5]][i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - 1 - se.mi[1]) * std::pow(y, se.mi[1]) * std::pow(z, order);
		}

	}
	for (uint i = 0; i < smallEdgesInBodies.size(); ++i) { //body
		const SmallEdge& se = smallEdgesInBodies[i];
		if (se.face == 0) {
			result.x += bodyCoefficients[element][i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]) *
				std::pow(z_d, order - se.mi[2]) * std::pow(z, se.mi[2]);
		}
		else if (se.face == 1) {
			result.y += bodyCoefficients[element][i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - 1 - se.mi[1]) * std::pow(y, se.mi[1]) *
				std::pow(z_d, order - se.mi[2]) * std::pow(z, se.mi[2]);
		}
		else {
			result.z += bodyCoefficients[element][i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]) *
				std::pow(z_d, order - 1 - se.mi[2]) * std::pow(z, se.mi[2]);
		}
	}
	result.x /= lenx;
	result.y /= leny;
	result.z /= lenz;
	return result;
}

/*
Computes the discrete Hodge for higher order cubical 1-forms to a sparse matrix (of type Buffer<std::unordered_map<uint, double>>).
Integrals over dual faces are computed once in a single element, so all elements are assumed to have the same shape.
*/
void SmallCubePartition3D::formHodgeMatrix1Forms(Buffer<std::unordered_map<uint, double>>& star, bool minkowskiMetric) const {
	if (!mesh_old_ptr)
		return; //SmallCubePartition3D has not been associated with a mesh

	const BuilderMesh& mesh_old = *mesh_old_ptr;
	const BuilderMesh& mesh = *mesh_ptr;
	star.resize(mesh.getEdgeSize());

	//perform computations in body 0
	uint refElement = 0;
	const Buffer<uint> refElementNodes = getCubeNodes(0, mesh_old);
	Buffer<Vector3> refElementPos(8);
	for (uint i = 0; i < refElementPos.size(); ++i)
		refElementPos[i] = mesh_old.getNodePosition3(refElementNodes[i]);
	Buffer<uint> refElementEdges(12);
	refElementEdges[0] = mesh_old.findEdge(refElementNodes[0], refElementNodes[1]);
	refElementEdges[1] = mesh_old.findEdge(refElementNodes[1], refElementNodes[2]);
	refElementEdges[2] = mesh_old.findEdge(refElementNodes[2], refElementNodes[3]);
	refElementEdges[3] = mesh_old.findEdge(refElementNodes[0], refElementNodes[3]);
	refElementEdges[4] = mesh_old.findEdge(refElementNodes[0], refElementNodes[4]);
	refElementEdges[5] = mesh_old.findEdge(refElementNodes[1], refElementNodes[5]);
	refElementEdges[6] = mesh_old.findEdge(refElementNodes[2], refElementNodes[6]);
	refElementEdges[7] = mesh_old.findEdge(refElementNodes[3], refElementNodes[7]);
	refElementEdges[8] = mesh_old.findEdge(refElementNodes[4], refElementNodes[5]);
	refElementEdges[9] = mesh_old.findEdge(refElementNodes[5], refElementNodes[6]);
	refElementEdges[10] = mesh_old.findEdge(refElementNodes[6], refElementNodes[7]);
	refElementEdges[11] = mesh_old.findEdge(refElementNodes[4], refElementNodes[7]);
	Vector3 refElementEdgeX = refElementPos[1] - refElementPos[0];
	Vector3 refElementEdgeY = refElementPos[3] - refElementPos[0];
	Vector3 refElementEdgeZ = refElementPos[4] - refElementPos[0];
	double refElementLenX = refElementEdgeX.len() / order;
	double refElementLenY = refElementEdgeY.len() / order;
	double refElementLenZ = refElementEdgeZ.len() / order;
	Buffer<uint> refElementFaces(6);
	refElementFaces[0] = findQuadrilateral(refElementNodes[0], refElementNodes[1], refElementNodes[2], refElementNodes[3], mesh_old);
	refElementFaces[1] = findQuadrilateral(refElementNodes[0], refElementNodes[1], refElementNodes[4], refElementNodes[5], mesh_old);
	refElementFaces[2] = findQuadrilateral(refElementNodes[1], refElementNodes[2], refElementNodes[5], refElementNodes[6], mesh_old);
	refElementFaces[3] = findQuadrilateral(refElementNodes[2], refElementNodes[3], refElementNodes[6], refElementNodes[7], mesh_old);
	refElementFaces[4] = findQuadrilateral(refElementNodes[0], refElementNodes[3], refElementNodes[4], refElementNodes[7], mesh_old);
	refElementFaces[5] = findQuadrilateral(refElementNodes[4], refElementNodes[5], refElementNodes[6], refElementNodes[7], mesh_old);

	//form the list of dual faces
	uint dualFaceCount = 12 * order + 12 * order * (order - 1) + 3 * order * (order - 1) * (order - 1);
	struct DualFacePart
	{
		Vector3 p0;
		Vector3 p1;
		Vector3 p2;
		Vector3 p3;
		uint face;
		double area;
	};
	std::vector<DualFacePart> dualFaceList;
	dualFaceList.resize(dualFaceCount);
	uint dfIndex = 0;
	//edge 0
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p0 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[0], i));
		Vector3 p1 = p0 + refElementEdgeY / (2 * order);
		Vector3 p2 = p1 + refElementEdgeZ / (2 * order);
		Vector3 p3 = p0 + refElementEdgeZ / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ / 4 };
	}
	//edge 1
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p1 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[1], i));
		Vector3 p0 = p1 - refElementEdgeX / (2 * order);
		Vector3 p2 = p1 + refElementEdgeZ / (2 * order);
		Vector3 p3 = p0 + refElementEdgeZ / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ / 4 };
	}
	//edge 2
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p1 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[2], i));
		Vector3 p0 = p1 - refElementEdgeY / (2 * order);
		Vector3 p2 = p1 + refElementEdgeZ / (2 * order);
		Vector3 p3 = p0 + refElementEdgeZ / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ / 4 };
	}
	//edge 3
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p0 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[3], i));
		Vector3 p1 = p0 + refElementEdgeX / (2 * order);
		Vector3 p2 = p1 + refElementEdgeZ / (2 * order);
		Vector3 p3 = p0 + refElementEdgeZ / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ / 4 };
	}
	//edge 4
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p0 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[4], i));
		Vector3 p1 = p0 + refElementEdgeX / (2 * order);
		Vector3 p2 = p1 + refElementEdgeY / (2 * order);
		Vector3 p3 = p0 + refElementEdgeY / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY / 4 };
	}
	//edge 5
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p1 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[5], i));
		Vector3 p0 = p1 - refElementEdgeX / (2 * order);
		Vector3 p2 = p1 + refElementEdgeY / (2 * order);
		Vector3 p3 = p0 + refElementEdgeY / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY / 4 };
	}
	//edge 6
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p2 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[6], i));
		Vector3 p3 = p2 - refElementEdgeX / (2 * order);
		Vector3 p1 = p2 - refElementEdgeY / (2 * order);
		Vector3 p0 = p3 - refElementEdgeY / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY / 4 };
	}
	//edge 7
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p3 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[7], i));
		Vector3 p2 = p3 + refElementEdgeX / (2 * order);
		Vector3 p0 = p3 - refElementEdgeY / (2 * order);
		Vector3 p1 = p0 + refElementEdgeX / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY / 4 };
	}
	//edge 8
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p3 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[8], i));
		Vector3 p2 = p3 + refElementEdgeY / (2 * order);
		Vector3 p0 = p3 - refElementEdgeZ / (2 * order);
		Vector3 p1 = p2 - refElementEdgeZ / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ / 4 };
	}
	//edge 9
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p2 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[9], i));
		Vector3 p3 = p2 - refElementEdgeX / (2 * order);
		Vector3 p0 = p3 - refElementEdgeZ / (2 * order);
		Vector3 p1 = p2 - refElementEdgeZ / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ / 4 };
	}
	//edge 10
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p2 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[10], i));
		Vector3 p3 = p2 - refElementEdgeY / (2 * order);
		Vector3 p0 = p3 - refElementEdgeZ / (2 * order);
		Vector3 p1 = p2 - refElementEdgeZ / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ / 4 };
	}
	//edge 11
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		Vector3 p3 = mesh.getEdgeAverage3(smallEdgesEdgeList(refElementEdges[11], i));
		Vector3 p2 = p3 + refElementEdgeX / (2 * order);
		Vector3 p0 = p3 - refElementEdgeZ / (2 * order);
		Vector3 p1 = p2 - refElementEdgeZ / (2 * order);
		dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ / 4 };
	}
	//face0
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		Vector3 edgeMidpoint = mesh.getEdgeAverage3(smallEdgesFaceList(refElementFaces[0], i));
		uint face = smallEdgesInFaces[i].face;
		if (face == 0) {
			Vector3 p0 = edgeMidpoint - refElementEdgeY / (2 * order);
			Vector3 p1 = edgeMidpoint + refElementEdgeY / (2 * order);
			Vector3 p3 = p0 + refElementEdgeZ / (2 * order);
			Vector3 p2 = p1 + refElementEdgeZ / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ / 2 };
		}
		else {
			Vector3 p0 = edgeMidpoint - refElementEdgeX / (2 * order);
			Vector3 p1 = edgeMidpoint + refElementEdgeX / (2 * order);
			Vector3 p3 = p0 + refElementEdgeZ / (2 * order);
			Vector3 p2 = p1 + refElementEdgeZ / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ / 2 };
		}
	}
	//face1
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		Vector3 edgeMidpoint = mesh.getEdgeAverage3(smallEdgesFaceList(refElementFaces[1], i));
		uint face = smallEdgesInFaces[i].face;
		if (face == 0) {
			Vector3 p0 = edgeMidpoint - refElementEdgeZ / (2 * order);
			Vector3 p3 = edgeMidpoint + refElementEdgeZ / (2 * order);
			Vector3 p1 = p0 + refElementEdgeY / (2 * order);
			Vector3 p2 = p3 + refElementEdgeY / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ / 2 };
		}
		else {
			Vector3 p0 = edgeMidpoint - refElementEdgeX / (2 * order);
			Vector3 p1 = edgeMidpoint + refElementEdgeX / (2 * order);
			Vector3 p3 = p0 + refElementEdgeY / (2 * order);
			Vector3 p2 = p1 + refElementEdgeY / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY / 2 };
		}
	}
	//face2
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		Vector3 edgeMidpoint = mesh.getEdgeAverage3(smallEdgesFaceList(refElementFaces[2], i));
		uint face = smallEdgesInFaces[i].face;
		if (face == 0) {
			Vector3 p1 = edgeMidpoint - refElementEdgeZ / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeZ / (2 * order);
			Vector3 p0 = p1 - refElementEdgeX / (2 * order);
			Vector3 p3 = p2 - refElementEdgeX / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ / 2 };
		}
		else {
			Vector3 p1 = edgeMidpoint - refElementEdgeY / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeY / (2 * order);
			Vector3 p0 = p1 - refElementEdgeX / (2 * order);
			Vector3 p3 = p2 - refElementEdgeX / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY / 2 };
		}
	}
	//face3
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		Vector3 edgeMidpoint = mesh.getEdgeAverage3(smallEdgesFaceList(refElementFaces[3], i));
		uint face = smallEdgesInFaces[i].face;
		if (face == 0) {
			Vector3 p1 = edgeMidpoint - refElementEdgeZ / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeZ / (2 * order);
			Vector3 p0 = p1 - refElementEdgeY / (2 * order);
			Vector3 p3 = p2 - refElementEdgeY / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ / 2 };
		}
		else {
			Vector3 p3 = edgeMidpoint - refElementEdgeX / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeX / (2 * order);
			Vector3 p0 = p3 - refElementEdgeY / (2 * order);
			Vector3 p1 = p2 - refElementEdgeY / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY / 2 };
		}
	}
	//face4
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		Vector3 edgeMidpoint = mesh.getEdgeAverage3(smallEdgesFaceList(refElementFaces[4], i));
		uint face = smallEdgesInFaces[i].face;
		if (face == 0) {
			Vector3 p0 = edgeMidpoint - refElementEdgeZ / (2 * order);
			Vector3 p3 = edgeMidpoint + refElementEdgeZ / (2 * order);
			Vector3 p1 = p0 + refElementEdgeX / (2 * order);
			Vector3 p2 = p3 + refElementEdgeX / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ / 2 };
		}
		else {
			Vector3 p0 = edgeMidpoint - refElementEdgeY / (2 * order);
			Vector3 p3 = edgeMidpoint + refElementEdgeY / (2 * order);
			Vector3 p1 = p0 + refElementEdgeX / (2 * order);
			Vector3 p2 = p3 + refElementEdgeX / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY / 2 };
		}
	}
	//face5
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		Vector3 edgeMidpoint = mesh.getEdgeAverage3(smallEdgesFaceList(refElementFaces[5], i));
		uint face = smallEdgesInFaces[i].face;
		if (face == 0) {
			Vector3 p3 = edgeMidpoint - refElementEdgeY / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeY / (2 * order);
			Vector3 p0 = p3 - refElementEdgeZ / (2 * order);
			Vector3 p1 = p2 - refElementEdgeZ / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ / 2 };
		}
		else {
			Vector3 p3 = edgeMidpoint - refElementEdgeX / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeX / (2 * order);
			Vector3 p0 = p3 - refElementEdgeZ / (2 * order);
			Vector3 p1 = p2 - refElementEdgeZ / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ / 2 };
		}
	}
	//body
	for (uint i = 0; i < smallEdgesInBodies.size(); ++i) {
		Vector3 edgeMidpoint = mesh.getEdgeAverage3(smallEdgesBodyList(0, i));
		uint face = smallEdgesInBodies[i].face;
		if (face == 0) {
			Vector3 p0 = edgeMidpoint - refElementEdgeY / (2 * order) - refElementEdgeZ / (2 * order);
			Vector3 p1 = edgeMidpoint + refElementEdgeY / (2 * order) - refElementEdgeZ / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeY / (2 * order) + refElementEdgeZ / (2 * order);
			Vector3 p3 = edgeMidpoint - refElementEdgeY / (2 * order) + refElementEdgeZ / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 0, refElementLenY * refElementLenZ };
		}
		else if (face == 1) {
			Vector3 p0 = edgeMidpoint - refElementEdgeX / (2 * order) - refElementEdgeZ / (2 * order);
			Vector3 p1 = edgeMidpoint + refElementEdgeX / (2 * order) - refElementEdgeZ / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeX / (2 * order) + refElementEdgeZ / (2 * order);
			Vector3 p3 = edgeMidpoint - refElementEdgeX / (2 * order) + refElementEdgeZ / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 1, refElementLenX * refElementLenZ };
		}
		else {
			Vector3 p0 = edgeMidpoint - refElementEdgeX / (2 * order) - refElementEdgeY / (2 * order);
			Vector3 p1 = edgeMidpoint + refElementEdgeX / (2 * order) - refElementEdgeY / (2 * order);
			Vector3 p2 = edgeMidpoint + refElementEdgeX / (2 * order) + refElementEdgeY / (2 * order);
			Vector3 p3 = edgeMidpoint - refElementEdgeX / (2 * order) + refElementEdgeY / (2 * order);
			dualFaceList[dfIndex++] = { p0, p1, p2, p3, 2, refElementLenX * refElementLenY };
		}
	}

	if (minkowskiMetric) {
		for (uint i = 0; i < dfIndex; ++i) {
			if (dualFaceList[i].face == 2)
				dualFaceList[i].area *= -1.0;
		}
	}

	//compute the integrals of basis functions
	Buffer<VectorN> hodgeIntegralsBodyEdges(dualFaceCount);
	Buffer<Buffer<VectorN>> hodgeIntegralsFaceEdges(6);
	Buffer<Buffer<VectorN>> hodgeIntegralsEdgeEdges(12);
	for (uint i = 0; i < dualFaceCount; ++i) {
		hodgeIntegralsBodyEdges[i].toVectorN(smallEdgesInBodies.size());
	}
	for (uint i = 0; i < 6; ++i) {
		hodgeIntegralsFaceEdges[i].resize(dualFaceCount);
		for (uint j = 0; j < dualFaceCount; ++j) {
			hodgeIntegralsFaceEdges[i][j].toVectorN(smallEdgesInFaces.size());
		}
	}
	for (uint i = 0; i < 12; ++i) {
		hodgeIntegralsEdgeEdges[i].resize(dualFaceCount);
		for (uint j = 0; j < dualFaceCount; ++j) {
			hodgeIntegralsEdgeEdges[i][j].toVectorN(smallEdgesInEdges.size());
		}
	}

	//first integrate body basis functions
	VectorN bodyCochain(smallEdgesInBodies.size(), 0.0);
	VectorN bodyCoefficients(smallEdgesInBodies.size(), 0.0);
	for (uint i = 0; i < smallEdgesInBodies.size(); ++i) {
		bodyCochain[i] = 1.0;
		solveLUP(matrix1FormsBodies, matrix1FormsBodies_p, bodyCochain, bodyCoefficients);
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (smallEdgesInBodies[i].face != face)
				continue;
			std::function<double(Vector3)> interpolant;
			if (face == 0)
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithBodyCoefficients(p, refElementPos, bodyCoefficients).x; };
			else if (face == 1)
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithBodyCoefficients(p, refElementPos, bodyCoefficients).y; };
			else
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithBodyCoefficients(p, refElementPos, bodyCoefficients).z; };
			hodgeIntegralsBodyEdges[j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));

		}
		bodyCoefficients.val.fill(0.0);
		bodyCochain[i] = 0.0;
	}

	//next integrate face basis functions
	VectorN faceCochain(smallEdgesInFaces.size(), 0.0);
	VectorN faceCoefficients(smallEdgesInFaces.size(), 0.0);
	VectorN bodyValues(smallEdgesInBodies.size(), 0.0);
	for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
		faceCochain[i] = 1.0;
		solveLUP(matrix1FormsFaces, matrix1FormsFaces_p, faceCochain, faceCoefficients);
		//face0
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyFaceValues1Forms[0][k].dot(faceCoefficients);
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if ((smallEdgesInFaces[i].face == 0 && face != 0) || (smallEdgesInFaces[i].face == 1 && face != 1))
				continue;
			std::function<double(Vector3)> interpolant;
			if (face == 0)
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 0).x; };
			else
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 0).y; };
			hodgeIntegralsFaceEdges[0][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsFaceEdges[0][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
		}
		bodyValues.val.fill(0.0);
		//face1
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyFaceValues1Forms[1][k].dot(faceCoefficients);
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if ((smallEdgesInFaces[i].face == 0 && face != 0) || (smallEdgesInFaces[i].face == 1 && face != 2))
				continue;
			std::function<double(Vector3)> interpolant;
			if (face == 0)
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 1).x; };
			else
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 1).z; };
			hodgeIntegralsFaceEdges[1][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsFaceEdges[1][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
		}
		bodyValues.val.fill(0.0);
		//face2
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyFaceValues1Forms[2][k].dot(faceCoefficients);
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if ((smallEdgesInFaces[i].face == 0 && face != 1) || (smallEdgesInFaces[i].face == 1 && face != 2))
				continue;
			std::function<double(Vector3)> interpolant;
			if (face == 1)
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 2).y; };
			else
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 2).z; };
			hodgeIntegralsFaceEdges[2][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsFaceEdges[2][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
		}
		bodyValues.val.fill(0.0);
		//face3
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyFaceValues1Forms[3][k].dot(faceCoefficients);
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if ((smallEdgesInFaces[i].face == 0 && face != 0) || (smallEdgesInFaces[i].face == 1 && face != 2))
				continue;
			std::function<double(Vector3)> interpolant;
			if (face == 0)
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 3).x; };
			else
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 3).z; };
			hodgeIntegralsFaceEdges[3][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsFaceEdges[3][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
		}
		bodyValues.val.fill(0.0);
		//face4
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyFaceValues1Forms[4][k].dot(faceCoefficients);
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if ((smallEdgesInFaces[i].face == 0 && face != 1) || (smallEdgesInFaces[i].face == 1 && face != 2))
				continue;
			std::function<double(Vector3)> interpolant;
			if (face == 1)
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 4).y; };
			else
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 4).z; };
			hodgeIntegralsFaceEdges[4][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsFaceEdges[4][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
		}
		bodyValues.val.fill(0.0);
		//face5
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyFaceValues1Forms[5][k].dot(faceCoefficients);
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if ((smallEdgesInFaces[i].face == 0 && face != 0) || (smallEdgesInFaces[i].face == 1 && face != 1))
				continue;
			std::function<double(Vector3)> interpolant;
			if (face == 0)
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 5).x; };
			else
				interpolant = [&](Vector3 p) -> double { return evaluate1FormWithFaceCoefficients(p, refElementPos, faceCoefficients, 5).y; };
			hodgeIntegralsFaceEdges[5][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsFaceEdges[5][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
		}
		bodyValues.val.fill(0.0);
		faceCoefficients.val.fill(0.0);
		faceCochain[i] = 0.0;
	}

	//finally integrate edge basis functions
	VectorN edgeCochain(smallEdgesInEdges.size(), 0.0);
	VectorN edgeCoefficients(smallEdgesInEdges.size(), 0.0);
	VectorN face0Values(smallEdgesInFaces.size(), 0.0);
	VectorN face1Values(smallEdgesInFaces.size(), 0.0);
	for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
		edgeCochain[i] = 1.0;
		solveLUP(matrix1FormsEdges, matrix1FormsEdges_p, edgeCochain, edgeCoefficients);
		//edge0
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[0][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[0][k].dot(edgeCoefficients); //0 in face 0
			face1Values[k] += faceEdgeValues1Forms[0][k].dot(edgeCoefficients); //0 in face 1
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 0)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 0).x; };
			hodgeIntegralsEdgeEdges[0][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[0][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[0][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[0][j]);
			hodgeIntegralsEdgeEdges[0][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[1][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge1
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[1][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[1][k].dot(edgeCoefficients); //1 in face 0
			face1Values[k] += faceEdgeValues1Forms[0][k].dot(edgeCoefficients); //0 in face 2
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 1)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 1).y; };
			hodgeIntegralsEdgeEdges[1][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[1][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[1][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[0][j]);
			hodgeIntegralsEdgeEdges[1][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[2][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge2
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[2][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[2][k].dot(edgeCoefficients); //2 in face 0
			face1Values[k] += faceEdgeValues1Forms[0][k].dot(edgeCoefficients); //0 in face 3
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 0)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 2).x; };
			hodgeIntegralsEdgeEdges[2][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[2][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[2][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[0][j]);
			hodgeIntegralsEdgeEdges[2][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[3][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge3
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[3][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[3][k].dot(edgeCoefficients); //3 in face 0
			face1Values[k] += faceEdgeValues1Forms[0][k].dot(edgeCoefficients); //0 in face 4
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 1)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 3).y; };
			hodgeIntegralsEdgeEdges[3][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[3][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[3][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[0][j]);
			hodgeIntegralsEdgeEdges[3][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[4][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge4
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[4][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[3][k].dot(edgeCoefficients); //3 in face 1
			face1Values[k] += faceEdgeValues1Forms[3][k].dot(edgeCoefficients); //3 in face 4
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 2)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 4).z; };
			hodgeIntegralsEdgeEdges[4][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[4][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[4][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[1][j]);
			hodgeIntegralsEdgeEdges[4][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[4][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge5
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[5][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[1][k].dot(edgeCoefficients); //1 in face 1
			face1Values[k] += faceEdgeValues1Forms[3][k].dot(edgeCoefficients); //3 in face 2
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 2)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 5).z; };
			hodgeIntegralsEdgeEdges[5][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[5][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[5][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[1][j]);
			hodgeIntegralsEdgeEdges[5][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[2][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge6
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[6][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[1][k].dot(edgeCoefficients); //1 in face 2
			face1Values[k] += faceEdgeValues1Forms[1][k].dot(edgeCoefficients); //1 in face 3
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 2)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 6).z; };
			hodgeIntegralsEdgeEdges[6][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[6][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[6][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[2][j]);
			hodgeIntegralsEdgeEdges[6][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[3][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge7
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[7][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[3][k].dot(edgeCoefficients); //3 in face 3
			face1Values[k] += faceEdgeValues1Forms[1][k].dot(edgeCoefficients); //1 in face 4
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 2)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 7).z; };
			hodgeIntegralsEdgeEdges[7][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[7][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[7][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[3][j]);
			hodgeIntegralsEdgeEdges[7][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[4][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge8
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[8][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[2][k].dot(edgeCoefficients); //2 in face 1
			face1Values[k] += faceEdgeValues1Forms[0][k].dot(edgeCoefficients); //0 in face 5
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 0)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 8).x; };
			hodgeIntegralsEdgeEdges[8][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[8][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[8][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[1][j]);
			hodgeIntegralsEdgeEdges[8][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[5][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge9
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[9][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[2][k].dot(edgeCoefficients); //2 in face 2
			face1Values[k] += faceEdgeValues1Forms[1][k].dot(edgeCoefficients); //1 in face 5
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 1)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 9).y; };
			hodgeIntegralsEdgeEdges[9][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[9][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[9][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[2][j]);
			hodgeIntegralsEdgeEdges[9][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[5][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge10
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[10][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[2][k].dot(edgeCoefficients); //2 in face 3
			face1Values[k] += faceEdgeValues1Forms[2][k].dot(edgeCoefficients); //2 in face 5
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 0)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 10).x; };
			hodgeIntegralsEdgeEdges[10][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[10][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[10][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[3][j]);
			hodgeIntegralsEdgeEdges[10][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[5][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		//edge11
		for (uint k = 0; k < smallEdgesInBodies.size(); ++k) {
			bodyValues[k] += bodyEdgeValues1Forms[11][k].dot(edgeCoefficients);
		}
		for (uint k = 0; k < smallEdgesInFaces.size(); ++k) {
			face0Values[k] += faceEdgeValues1Forms[2][k].dot(edgeCoefficients); //2 in face 4
			face1Values[k] += faceEdgeValues1Forms[3][k].dot(edgeCoefficients); //3 in face 5
		}
		for (uint j = 0; j < dualFaceCount; ++j) {
			uint face = dualFaceList[j].face;
			if (face != 1)
				continue;
			std::function<double(Vector3)> interpolant = [&](Vector3 p) -> double { return evaluate1FormWithEdgeCoefficients(p, refElementPos, edgeCoefficients, 11).y; };
			hodgeIntegralsEdgeEdges[11][j][i] += 0.5 * dualFaceList[j].area * (integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p1, dualFaceList[j].p2) +
				integralAverage(interpolant, dualFaceList[j].p0, dualFaceList[j].p2, dualFaceList[j].p3));
			hodgeIntegralsEdgeEdges[11][j][i] -= bodyValues.dot(hodgeIntegralsBodyEdges[j]);
			hodgeIntegralsEdgeEdges[11][j][i] -= face0Values.dot(hodgeIntegralsFaceEdges[4][j]);
			hodgeIntegralsEdgeEdges[11][j][i] -= face1Values.dot(hodgeIntegralsFaceEdges[5][j]);
		}
		bodyValues.val.fill(0.0);
		face0Values.val.fill(0.0);
		face1Values.val.fill(0.0);
		edgeCoefficients.val.fill(0.0);
		edgeCochain[i] = 0.0;
	}

	//form Hodge matrix from the precomputed integrals
	for (uint b = 0; b < mesh_old.getBodySize(); ++b) {
		const Buffer<uint> nodes = getCubeNodes(b, mesh_old);
		Buffer<Vector3> pos(8);
		for (uint i = 0; i < pos.size(); ++i)
			pos[i] = mesh_old.getNodePosition3(nodes[i]);
		Buffer<uint> edges(12);
		edges[0] = mesh_old.findEdge(nodes[0], nodes[1]);
		edges[1] = mesh_old.findEdge(nodes[1], nodes[2]);
		edges[2] = mesh_old.findEdge(nodes[2], nodes[3]);
		edges[3] = mesh_old.findEdge(nodes[0], nodes[3]);
		edges[4] = mesh_old.findEdge(nodes[0], nodes[4]);
		edges[5] = mesh_old.findEdge(nodes[1], nodes[5]);
		edges[6] = mesh_old.findEdge(nodes[2], nodes[6]);
		edges[7] = mesh_old.findEdge(nodes[3], nodes[7]);
		edges[8] = mesh_old.findEdge(nodes[4], nodes[5]);
		edges[9] = mesh_old.findEdge(nodes[5], nodes[6]);
		edges[10] = mesh_old.findEdge(nodes[6], nodes[7]);
		edges[11] = mesh_old.findEdge(nodes[4], nodes[7]);
		Vector3 edgeX = pos[1] - pos[0];
		Vector3 edgeY = pos[3] - pos[0];
		Vector3 edgeZ = pos[4] - pos[0];
		double lenX = edgeX.len() / order;
		double lenY = edgeY.len() / order;
		double lenZ = edgeZ.len() / order;
		Buffer<uint> faces(6);
		faces[0] = findQuadrilateral(nodes[0], nodes[1], nodes[2], nodes[3], mesh_old);
		faces[1] = findQuadrilateral(nodes[0], nodes[1], nodes[4], nodes[5], mesh_old);
		faces[2] = findQuadrilateral(nodes[1], nodes[2], nodes[5], nodes[6], mesh_old);
		faces[3] = findQuadrilateral(nodes[2], nodes[3], nodes[6], nodes[7], mesh_old);
		faces[4] = findQuadrilateral(nodes[0], nodes[3], nodes[4], nodes[7], mesh_old);
		faces[5] = findQuadrilateral(nodes[4], nodes[5], nodes[6], nodes[7], mesh_old);

		//find the indices
		Buffer<uint> indices(dualFaceCount);
		dfIndex = 0;
		//edge 0
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[0], i);
		}
		//edge 1
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[1], i);
		}
		//edge 2
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[2], i);
		}
		//edge 3
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[3], i);
		}
		//edge 4
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[4], i);
		}
		//edge 5
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[5], i);
		}
		//edge 6
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[6], i);
		}
		//edge 7
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[7], i);
		}
		//edge 8
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[8], i);
		}
		//edge 9
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[9], i);
		}
		//edge 10
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[10], i);
		}
		//edge 11
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			indices[dfIndex++] = smallEdgesEdgeList(edges[11], i);
		}
		//face0
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			indices[dfIndex++] = smallEdgesFaceList(faces[0], i);
		}
		//face1
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			indices[dfIndex++] = smallEdgesFaceList(faces[1], i);
		}
		//face2
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			indices[dfIndex++] = smallEdgesFaceList(faces[2], i);
		}
		//face3
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			indices[dfIndex++] = smallEdgesFaceList(faces[3], i);
		}
		//face4
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			indices[dfIndex++] = smallEdgesFaceList(faces[4], i);
		}
		//face5
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			indices[dfIndex++] = smallEdgesFaceList(faces[5], i);
		}
		//body
		for (uint i = 0; i < smallEdgesInBodies.size(); ++i) {
			indices[dfIndex++] = smallEdgesBodyList(b, i);
		}

		for (uint i = 0; i < dualFaceCount; ++i) { //edge basis functions
			for (uint j = 0; j < smallEdgesInEdges.size(); ++j) {
				star[indices[i]][smallEdgesEdgeList(edges[0], j)] += hodgeIntegralsEdgeEdges[0][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[1], j)] += hodgeIntegralsEdgeEdges[1][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[2], j)] += hodgeIntegralsEdgeEdges[2][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[3], j)] += hodgeIntegralsEdgeEdges[3][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[4], j)] += hodgeIntegralsEdgeEdges[4][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[5], j)] += hodgeIntegralsEdgeEdges[5][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[6], j)] += hodgeIntegralsEdgeEdges[6][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[7], j)] += hodgeIntegralsEdgeEdges[7][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[8], j)] += hodgeIntegralsEdgeEdges[8][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[9], j)] += hodgeIntegralsEdgeEdges[9][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[10], j)] += hodgeIntegralsEdgeEdges[10][i][j];
				star[indices[i]][smallEdgesEdgeList(edges[11], j)] += hodgeIntegralsEdgeEdges[11][i][j];
			}
		}
		for (uint i = 0; i < dualFaceCount; ++i) { //face basis functions
			for (uint j = 0; j < smallEdgesInFaces.size(); ++j) {
				star[indices[i]][smallEdgesFaceList(faces[0], j)] += hodgeIntegralsFaceEdges[0][i][j];
				star[indices[i]][smallEdgesFaceList(faces[1], j)] += hodgeIntegralsFaceEdges[1][i][j];
				star[indices[i]][smallEdgesFaceList(faces[2], j)] += hodgeIntegralsFaceEdges[2][i][j];
				star[indices[i]][smallEdgesFaceList(faces[3], j)] += hodgeIntegralsFaceEdges[3][i][j];
				star[indices[i]][smallEdgesFaceList(faces[4], j)] += hodgeIntegralsFaceEdges[4][i][j];
				star[indices[i]][smallEdgesFaceList(faces[5], j)] += hodgeIntegralsFaceEdges[5][i][j];
			}
		}
		for (uint i = 0; i < dualFaceCount; ++i) { //body basis functions
			for (uint j = 0; j < smallEdgesInBodies.size(); ++j) {
				star[indices[i]][smallEdgesBodyList(b, j)] += hodgeIntegralsBodyEdges[i][j];
			}
		}
	}
}

//The following functions return the index of the small cell j of the big cell i.

uint SmallCubePartition3D::smallNodesEdgeList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + i * (order - 1) + j;
}
uint SmallCubePartition3D::smallEdgesEdgeList(uint i, uint j) const {
	return i * order + j;
}
uint SmallCubePartition3D::smallNodesFaceList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + mesh_old_ptr->getEdgeSize() * (order - 1) + i * (order - 1) * (order - 1) + j;
}
uint SmallCubePartition3D::smallEdgesFaceList(uint i, uint j) const {
	return mesh_old_ptr->getEdgeSize() * order + i * 2 * order * (order - 1) + j;
}
uint SmallCubePartition3D::smallFacesFaceList(uint i, uint j) const {
	return i * order * order + j;
}
uint SmallCubePartition3D::smallNodesBodyList(uint i, uint j) const {
	return mesh_old_ptr->getNodeSize() + mesh_old_ptr->getEdgeSize() * (order - 1) + mesh_old_ptr->getFaceSize() * (order - 1) * (order - 1) 
		+ i * (order - 1) * (order - 1) * (order - 1) + j;
}
uint SmallCubePartition3D::smallEdgesBodyList(uint i, uint j) const {
	return mesh_old_ptr->getEdgeSize() * order + mesh_old_ptr->getFaceSize() * 2 * order * (order - 1) + i * 3 * order * (order - 1) * (order - 1) + j;
}
uint SmallCubePartition3D::smallFacesBodyList(uint i, uint j) const {
	return mesh_old_ptr->getFaceSize() * order * order + i * 3 * order * order * (order - 1) + j;
}
uint SmallCubePartition3D::smallBodiesBodyList(uint i, uint j) const {
	return i * order * order * order + j;
}

//returns the index of the small node located at pos in the element of given dimension
uint SmallCubePartition3D::getSmallNodeIndex(const Vector3& pos, uint element, uint dim) {
	const BuilderMesh& mesh = *mesh_ptr;
	if (dim == 1) {
		for (uint i = 0; i < order - 1; ++i) {
			uint node = smallNodesEdgeList(element, i);
			if ((mesh.getNodePosition3(node) - pos).lensq() < 1e-24)
				return node;
		}
	}
	else if (dim == 2) {
		for (uint i = 0; i < (order - 1) * (order - 1); ++i) {
			uint node = smallNodesFaceList(element, i);
			if ((mesh.getNodePosition3(node) - pos).lensq() < 1e-24)
				return node;
		}
	}
	else if (dim == 3) {
		for (uint i = 0; i < (order - 1) * (order - 1) * (order - 1); ++i) {
			uint node = smallNodesBodyList(element, i);
			if ((mesh.getNodePosition3(node) - pos).lensq() < 1e-24)
				return node;
		}
	}
	return NONE;
}

//computes the value resulting from the coefficients on the given edge in the element whose nodes are given in pos
Vector3 SmallCubePartition3D::evaluate1FormWithEdgeCoefficients(const Vector3& evaluationPoint, const Buffer<Vector3>& pos, const VectorN& edgeCoefficients, uint edge) const {
	Vector3 result(0.0, 0.0, 0.0);
	double lenx = pos[6].x - pos[0].x;
	double leny = pos[6].y - pos[0].y;
	double lenz = pos[6].z - pos[0].z;
	double x = (evaluationPoint.x - pos[0].x) / lenx;
	double y = (evaluationPoint.y - pos[0].y) / leny;
	double z = (evaluationPoint.z - pos[0].z) / lenz;
	double x_d = (pos[6].x - evaluationPoint.x) / lenx;
	double y_d = (pos[6].y - evaluationPoint.y) / leny;
	double z_d = (pos[6].z - evaluationPoint.z) / lenz;
	double sum = 0.0;
	if (edge == 0) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(y_d, order) * std::pow(z_d, order);
		result.x += sum;
		sum = 0.0;
	}
	else if (edge == 1) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x, order) * std::pow(z_d, order);
		result.y += sum;
		sum = 0.0;
	}
	else if (edge == 2) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(y, order) * std::pow(z_d, order);
		result.x += sum;
		sum = 0.0;
	}
	else if (edge == 3) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x_d, order) * std::pow(z_d, order);
		result.y += sum;
		sum = 0.0;
	}
	else if (edge == 4) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(z_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(z, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x_d, order) * std::pow(y_d, order);
		result.z += sum;
		sum = 0.0;
	}
	else if (edge == 5) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(z_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(z, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x, order) * std::pow(y_d, order);
		result.z += sum;
		sum = 0.0;
	}
	else if (edge == 6) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(z_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(z, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x, order) * std::pow(y, order);
		result.z += sum;
		sum = 0.0;
	}
	else if (edge == 7) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(z_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(z, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x_d, order) * std::pow(y, order);
		result.z += sum;
		sum = 0.0;
	}
	else if (edge == 8) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(y_d, order) * std::pow(z, order);
		result.x += sum;
		sum = 0.0;
	}
	else if (edge == 9) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x, order) * std::pow(z, order);
		result.y += sum;
		sum = 0.0;
	}
	else if (edge == 10) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(x_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(x, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(y, order) * std::pow(z, order);
		result.x += sum;
		sum = 0.0;
	}
	else if (edge == 11) {
		for (uint i = 0; i < smallEdgesInEdges.size(); ++i) {
			sum += edgeCoefficients[i] * std::pow(y_d, order - 1 - smallEdgesInEdges[i].mi[0]) * std::pow(y, smallEdgesInEdges[i].mi[0]);
		}
		sum *= std::pow(x_d, order) * std::pow(z, order);
	}
	result.y += sum;
	result.x /= lenx;
	result.y /= leny;
	result.z /= lenz;
	return result;
}

//computes the value resulting from the coefficients on the given face in the element whose nodes are given in pos
Vector3 SmallCubePartition3D::evaluate1FormWithFaceCoefficients(const Vector3& evaluationPoint, const Buffer<Vector3>& pos, const VectorN& faceCoefficients, uint face) const {
	Vector3 result(0.0, 0.0, 0.0);
	double lenx = pos[6].x - pos[0].x;
	double leny = pos[6].y - pos[0].y;
	double lenz = pos[6].z - pos[0].z;
	double x = (evaluationPoint.x - pos[0].x) / lenx;
	double y = (evaluationPoint.y - pos[0].y) / leny;
	double z = (evaluationPoint.z - pos[0].z) / lenz;
	double x_d = (pos[6].x - evaluationPoint.x) / lenx;
	double y_d = (pos[6].y - evaluationPoint.y) / leny;
	double z_d = (pos[6].z - evaluationPoint.z) / lenz;
	if (face == 0) {
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			const SmallEdge& se = smallEdgesInFaces[i];
			if (se.face == 0) {
				result.x += faceCoefficients[i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]) * std::pow(z_d, order);
			}
			else {
				result.y += faceCoefficients[i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - 1 - se.mi[1]) * std::pow(y, se.mi[1]) * std::pow(z_d, order);
			}

		}
	}
	else if (face == 1) {
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			const SmallEdge& se = smallEdgesInFaces[i];
			if (se.face == 0) {
				result.x += faceCoefficients[i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(z_d, order - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(y_d, order);
			}
			else {
				result.z += faceCoefficients[i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(z_d, order - 1 - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(y_d, order);
			}

		}
	}
	else if (face == 2) {
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			const SmallEdge& se = smallEdgesInFaces[i];
			if (se.face == 0) {
				result.y += faceCoefficients[i] * std::pow(y_d, order - 1 - se.mi[0]) * std::pow(y, se.mi[0]) * std::pow(z_d, order - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(x, order);
			}
			else {
				result.z += faceCoefficients[i] * std::pow(y_d, order - se.mi[0]) * std::pow(y, se.mi[0]) * std::pow(z_d, order - 1 - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(x, order);
			}

		}
	}
	else if (face == 3) {
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			const SmallEdge& se = smallEdgesInFaces[i];
			if (se.face == 0) {
				result.x += faceCoefficients[i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(z_d, order - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(y, order);
			}
			else {
				result.z += faceCoefficients[i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(z_d, order - 1 - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(y, order);
			}

		}
	}
	else if (face == 4) {
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			const SmallEdge& se = smallEdgesInFaces[i];
			if (se.face == 0) {
				result.y += faceCoefficients[i] * std::pow(y_d, order - 1 - se.mi[0]) * std::pow(y, se.mi[0]) * std::pow(z_d, order - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(x_d, order);
			}
			else {
				result.z += faceCoefficients[i] * std::pow(y_d, order - se.mi[0]) * std::pow(y, se.mi[0]) * std::pow(z_d, order - 1 - se.mi[1]) * std::pow(z, se.mi[1]) * std::pow(x_d, order);
			}

		}
	}
	else if (face == 5) {
		for (uint i = 0; i < smallEdgesInFaces.size(); ++i) {
			const SmallEdge& se = smallEdgesInFaces[i];
			if (se.face == 0) {
				result.x += faceCoefficients[i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]) * std::pow(z, order);
			}
			else {
				result.y += faceCoefficients[i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - 1 - se.mi[1]) * std::pow(y, se.mi[1]) * std::pow(z, order);
			}

		}
	}
	result.x /= lenx;
	result.y /= leny;
	result.z /= lenz;
	return result;
}

//computes the value resulting from the coefficients in the element whose nodes are given in pos
Vector3 SmallCubePartition3D::evaluate1FormWithBodyCoefficients(const Vector3& evaluationPoint, const Buffer<Vector3>& pos, const VectorN& bodyCoefficients) const {
	Vector3 result(0.0, 0.0, 0.0);
	double lenx = pos[6].x - pos[0].x;
	double leny = pos[6].y - pos[0].y;
	double lenz = pos[6].z - pos[0].z;
	double x = (evaluationPoint.x - pos[0].x) / lenx;
	double y = (evaluationPoint.y - pos[0].y) / leny;
	double z = (evaluationPoint.z - pos[0].z) / lenz;
	double x_d = (pos[6].x - evaluationPoint.x) / lenx;
	double y_d = (pos[6].y - evaluationPoint.y) / leny;
	double z_d = (pos[6].z - evaluationPoint.z) / lenz;
	for (uint i = 0; i < smallEdgesInBodies.size(); ++i) {
		const SmallEdge& se = smallEdgesInBodies[i];
		if (se.face == 0) {
			result.x += bodyCoefficients[i] * std::pow(x_d, order - 1 - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]) *
				std::pow(z_d, order - se.mi[2]) * std::pow(z, se.mi[2]);
		}
		else if (se.face == 1) {
			result.y += bodyCoefficients[i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - 1 - se.mi[1]) * std::pow(y, se.mi[1]) *
				std::pow(z_d, order - se.mi[2]) * std::pow(z, se.mi[2]);
		}
		else {
			result.z += bodyCoefficients[i] * std::pow(x_d, order - se.mi[0]) * std::pow(x, se.mi[0]) * std::pow(y_d, order - se.mi[1]) * std::pow(y, se.mi[1]) *
				std::pow(z_d, order - 1 - se.mi[2]) * std::pow(z, se.mi[2]);
		}
	}
	result.x /= lenx;
	result.y /= leny;
	result.z /= lenz;
	return result;
}

//precomputes the matrices required in evaluating higher order cubical forms
void SmallCubePartition3D::formMatrices() {
	integrals1Forms.resize(smallEdgesInEdges.size());
	Buffer<Buffer<double>> integrals1FormsEdges(smallEdgesInEdges.size());
	Buffer<Buffer<double>> integrals1FormsFaces(smallEdgesInFaces.size());
	Buffer<Buffer<double>> integrals1FormsBodies(smallEdgesInBodies.size());

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
				if ((se_i.face == 0 && (k == 1 || k == 3)) || (se_i.face == 1 && (k == 0 || k == 2))) {
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
			if (se_i.face != se_j.face)
				integrals1FormsFaces[i][j] = 0.0;
			else {
				if (se_i.face == 0) { //dx-component
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

	//1-forms in bodies
	bodyEdgeValues1Forms.resize(12);
	for (uint k = 0; k < 12; ++k) {
		bodyEdgeValues1Forms[k].resize(smallEdgesInBodies.size());
		for (uint i = 0; i < smallEdgesInBodies.size(); ++i) {
			bodyEdgeValues1Forms[k][i].toVectorN(smallEdgesInEdges.size());
			const SmallEdge& se_i = smallEdgesInBodies[i]; //small edge i is the domain of integration
			for (uint j = 0; j < smallEdgesInEdges.size(); ++j) {
				const SmallEdge& se_j = smallEdgesInEdges[j]; //the form corresponding to small edge j (in edge k) is the integrand
				if ((se_i.face != 0 && (k == 0 || k == 2 || k == 8 || k == 10)) || (se_i.face != 1 && (k == 1 || k == 3 || k == 9 || k == 11)) ||
					(se_i.face != 2 && (k == 4 || k == 5 || k == 6 || k == 7))) {
					bodyEdgeValues1Forms[k][i][j] = 0.0;
				}
				else {
					if (se_i.face == 0) {
						double y = (double)(se_i.mi[1]) / (double)(order);
						double z = (double)(se_i.mi[2]) / (double)(order);
						if (k == 0)
							bodyEdgeValues1Forms[k][i][j] = std::pow(1 - y, order) * std::pow(1 - z, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
						else if (k == 2)
							bodyEdgeValues1Forms[k][i][j] = std::pow(y, order) * std::pow(1 - z, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
						else if (k == 8)
							bodyEdgeValues1Forms[k][i][j] = std::pow(1 - y, order) * std::pow(z, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
						else
							bodyEdgeValues1Forms[k][i][j] = std::pow(y, order) * std::pow(z, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
					}
					else if (se_i.face == 1) {
						double x = (double)(se_i.mi[0]) / (double)(order);
						double z = (double)(se_i.mi[2]) / (double)(order);
						if (k == 1)
							bodyEdgeValues1Forms[k][i][j] = std::pow(x, order) * std::pow(1 - z, order) * integrals1Forms[se_i.mi[1]][se_j.mi[0]];
						else if (k == 3)
							bodyEdgeValues1Forms[k][i][j] = std::pow(1 - x, order) * std::pow(1 - z, order) * integrals1Forms[se_i.mi[1]][se_j.mi[0]];
						else if (k == 9)
							bodyEdgeValues1Forms[k][i][j] = std::pow(x, order) * std::pow(z, order) * integrals1Forms[se_i.mi[1]][se_j.mi[0]];
						else
							bodyEdgeValues1Forms[k][i][j] = std::pow(1 - x, order) * std::pow(z, order) * integrals1Forms[se_i.mi[1]][se_j.mi[0]];
					}
					else {
						double x = (double)(se_i.mi[0]) / (double)(order);
						double y = (double)(se_i.mi[1]) / (double)(order);
						if (k == 4)
							bodyEdgeValues1Forms[k][i][j] = std::pow(1 - x, order) * std::pow(1 - y, order) * integrals1Forms[se_i.mi[2]][se_j.mi[0]];
						else if (k == 5)
							bodyEdgeValues1Forms[k][i][j] = std::pow(x, order) * std::pow(1 - y, order) * integrals1Forms[se_i.mi[2]][se_j.mi[0]];
						else if (k == 6)
							bodyEdgeValues1Forms[k][i][j] = std::pow(x, order) * std::pow(y, order) * integrals1Forms[se_i.mi[2]][se_j.mi[0]];
						else
							bodyEdgeValues1Forms[k][i][j] = std::pow(1 - x, order) * std::pow(y, order) * integrals1Forms[se_i.mi[2]][se_j.mi[0]];
					}
				}
			}
		}
	}
	bodyFaceValues1Forms.resize(6);
	for (uint k = 0; k < 6; ++k) {
		bodyFaceValues1Forms[k].resize(smallEdgesInBodies.size());
		for (uint i = 0; i < smallEdgesInBodies.size(); ++i) {
			bodyFaceValues1Forms[k][i].toVectorN(smallEdgesInFaces.size());
			const SmallEdge& se_i = smallEdgesInBodies[i]; //small edge i is the domain of integration
			for (uint j = 0; j < smallEdgesInFaces.size(); ++j) {
				const SmallEdge& se_j = smallEdgesInFaces[j]; //the form corresponding to small edge j (in face k) is the integrand
				if ((se_i.face != 0 && ((k == 0 && se_j.face == 0) || (k == 1 && se_j.face == 0) || (k == 3 && se_j.face == 0) || (k == 5 && se_j.face == 0))) ||
					(se_i.face != 1 && ((k == 0 && se_j.face == 1) || (k == 2 && se_j.face == 0) || (k == 4 && se_j.face == 0) || (k == 5 && se_j.face == 1))) ||
					(se_i.face != 2 && ((k == 1 && se_j.face == 1) || (k == 2 && se_j.face == 1) || (k == 3 && se_j.face == 1) || (k == 4 && se_j.face == 1)))) {
					bodyFaceValues1Forms[k][i][j] = 0.0;
				}
				else {
					if (se_i.face == 0) {
						double y = (double)(se_i.mi[1]) / (double)(order);
						double z = (double)(se_i.mi[2]) / (double)(order);
						if (k == 0) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - y, order - se_j.mi[1]) * std::pow(y, se_j.mi[1]) * std::pow(1 - z, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
						}
						else if (k == 1) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - z, order - se_j.mi[1]) * std::pow(z, se_j.mi[1]) * std::pow(1 - y, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
						}
						else if (k == 3) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - z, order - se_j.mi[1]) * std::pow(z, se_j.mi[1]) * std::pow(y, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
						}
						else {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - y, order - se_j.mi[1]) * std::pow(y, se_j.mi[1]) * std::pow(z, order) * integrals1Forms[se_i.mi[0]][se_j.mi[0]];
						}
					}
					else if (se_i.face == 1) {
						double x = (double)(se_i.mi[0]) / (double)(order);
						double z = (double)(se_i.mi[2]) / (double)(order);
						if (k == 0) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - x, order - se_j.mi[0]) * std::pow(x, se_j.mi[0]) * std::pow(1 - z, order) * integrals1Forms[se_i.mi[1]][se_j.mi[1]];
						}
						else if (k == 2) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - z, order - se_j.mi[1]) * std::pow(z, se_j.mi[1]) * std::pow(x, order) * integrals1Forms[se_i.mi[1]][se_j.mi[0]];
						}
						else if (k == 4) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - z, order - se_j.mi[1]) * std::pow(z, se_j.mi[1]) * std::pow(1 - x, order) * integrals1Forms[se_i.mi[1]][se_j.mi[0]];
						}
						else {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - x, order - se_j.mi[0]) * std::pow(x, se_j.mi[0]) * std::pow(z, order) * integrals1Forms[se_i.mi[1]][se_j.mi[1]];
						}
					}
					else {
						double x = (double)(se_i.mi[0]) / (double)(order);
						double y = (double)(se_i.mi[1]) / (double)(order);
						if (k == 1) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - x, order - se_j.mi[0]) * std::pow(x, se_j.mi[0]) * std::pow(1 - y, order) * integrals1Forms[se_i.mi[2]][se_j.mi[1]];
						}
						else if (k == 2) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - y, order - se_j.mi[0]) * std::pow(y, se_j.mi[0]) * std::pow(x, order) * integrals1Forms[se_i.mi[2]][se_j.mi[1]];
						}
						else if (k == 3) {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - x, order - se_j.mi[0]) * std::pow(x, se_j.mi[0]) * std::pow(y, order) * integrals1Forms[se_i.mi[2]][se_j.mi[1]];
						}
						else {
							bodyFaceValues1Forms[k][i][j] = std::pow(1 - y, order - se_j.mi[0]) * std::pow(y, se_j.mi[0]) * std::pow(1 - x, order) * integrals1Forms[se_i.mi[2]][se_j.mi[1]];
						}
					}
				}
			}
		}
	}
	for (uint i = 0; i < smallEdgesInBodies.size(); ++i) {
		integrals1FormsBodies[i].resize(smallEdgesInBodies.size());
		const SmallEdge& se_i = smallEdgesInBodies[i]; //small edge i is the domain of integration
		for (uint j = 0; j < smallEdgesInBodies.size(); ++j) {
			const SmallEdge& se_j = smallEdgesInBodies[j]; //the form corresponding to small edge j is the integrand
			if (se_i.face != se_j.face)
				integrals1FormsBodies[i][j] = 0.0;
			else {
				if (se_i.face == 0) { //dx-component
					double y = (double)(se_i.mi[1]) / (double)(order);
					double z = (double)(se_i.mi[2]) / (double)(order);
					integrals1FormsBodies[i][j] = std::pow(1 - y, order - se_j.mi[1]) * std::pow(y, se_j.mi[1]) * std::pow(1 - z, order - se_j.mi[2]) * std::pow(z, se_j.mi[2])
						* integrals1Forms[se_i.mi[0]][se_j.mi[0]];
				}
				else if (se_i.face == 1) { //dy-component
					double x = (double)(se_i.mi[0]) / (double)(order);
					double z = (double)(se_i.mi[2]) / (double)(order);
					integrals1FormsBodies[i][j] = std::pow(1 - x, order - se_j.mi[0]) * std::pow(x, se_j.mi[0]) * std::pow(1 - z, order - se_j.mi[2]) * std::pow(z, se_j.mi[2])
						* integrals1Forms[se_i.mi[1]][se_j.mi[1]];
				}
				else { //dz-component
					double x = (double)(se_i.mi[0]) / (double)(order);
					double y = (double)(se_i.mi[1]) / (double)(order);
					integrals1FormsBodies[i][j] = std::pow(1 - x, order - se_j.mi[0]) * std::pow(x, se_j.mi[0]) * std::pow(1 - y, order - se_j.mi[1]) * std::pow(y, se_j.mi[1])
						* integrals1Forms[se_i.mi[2]][se_j.mi[2]];
				}
			}
		}
	}

	//precompute LU decompositions
	double tol = 1e-15;
	decomposeLUP(integrals1FormsEdges, matrix1FormsEdges_p, tol);
	decomposeLUP(integrals1FormsFaces, matrix1FormsFaces_p, tol);
	decomposeLUP(integrals1FormsBodies, matrix1FormsBodies_p, tol);
	//convert to MatrixN
	convertMatrix(integrals1FormsEdges, matrix1FormsEdges);
	convertMatrix(integrals1FormsFaces, matrix1FormsFaces);
	convertMatrix(integrals1FormsBodies, matrix1FormsBodies);
}
