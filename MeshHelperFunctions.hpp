/*
MeshHelperFunctions contains some useful functions related to Mesh class.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _MESHHELPERFUNCTIONS_HPP_INCLUDED_
#define _MESHHELPERFUNCTIONS_HPP_INCLUDED_

#include "GFD/Types/Types.hpp"
#include "GFD/Types/Vector.hpp"
#include "GFD/Mesh/BuilderMesh.hpp"
#include <array>
#include <functional>

namespace gfd
{
	Buffer<uint> findEdges(const uint n0, const uint n1, const uint n2, const Mesh& mesh);
	Buffer<uint> findFaces(const uint n0, const uint n1, const uint n2, const uint n3, const Mesh& mesh);
	uint findSimplexWithNodes(const Buffer<uint>& simplexNodes, const Mesh& mesh);
	uint findQuadrilateral(const uint n0, const uint n1, const uint n2, const uint n3, const Mesh& mesh);
	Buffer<uint> getQuadrilateralNodes(uint i, const Mesh& mesh);
	Buffer<uint> getCubeNodes(uint i, const Mesh& mesh);
	void createOneElementMesh(Mesh& mesh);
	void createTwoElementMesh(Mesh& mesh);
	void formDiagonalHodge1Forms(Buffer<Buffer<double>>& star, const Mesh& mesh);
	bool isInsideTriangle(const Vector2& p, const Vector2& p0, const Vector2& p1, const Vector2& p2,
		const Vector2& p0oppositeNormal, const Vector2& p1oppositeNormal, const Vector2& p2oppositeNormal);
	bool isInsideTetrahedron(const Vector3& p, const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3,
		const Vector3& p0oppositeNormal, const Vector3& p1oppositeNormal, const Vector3& p2oppositeNormal, const Vector3& p3oppositeNormal);
	void printElementNodesForMathematica(const Mesh& mesh);
	void printEdgeLengthData(const Mesh& mesh);
	void printNodePositions(const Mesh& mesh);
	double dualEdgeLength(uint i, const Mesh& mesh, bool circumcentric = false);
	double dualFaceArea(uint i, const Mesh& mesh, bool circumcentric = false);
	double dualFaceBoundaryLength(uint i, const Mesh& mesh, bool circumcentric = false);
	double dualCellVolume(uint i, const Mesh& mesh, bool circumcentric = false);
	double dualCellVolumeCartesianMesh(uint i, const Mesh& mesh);
	double dualCellBoundaryArea(uint i, const Mesh& mesh, bool circumcentric = false);
	void createRhombicDodecahedronMesh(BuilderMesh& mesh, double r, double d);
	void createSpacetimeMesh(Mesh& mesh, double L, double T, uint xsteps, uint tsteps);
	void createSpacetimeMesh(Mesh& mesh, double Lx, double Ly, double T, uint xsteps, uint ysteps, uint tsteps);
	void createSpacetimeMeshExtraNode(Mesh& mesh, double Lx, double Ly, double T, uint xsteps, uint ysteps, uint tsteps);
	void createSpacetimeMeshExtraNodes(Mesh& mesh, double Lx, double Ly, double T, uint xsteps, uint ysteps, uint tsteps);
	void createSpacetimeMeshBccParallellepiped(Mesh& mesh, double h, uint xsteps, uint ysteps, uint tsteps);
	void createCartesianMesh(Mesh& mesh, double Lx, double T, uint xsteps, uint tsteps);
	void createCartesianMesh(Mesh& mesh, double Lx, double Ly, double T, uint xsteps, uint ysteps, uint tsteps);
	bool isXDir(const Vector3& vec);
	bool isYDir(const Vector3& vec);
	bool isZDir(const Vector3& vec);
	bool isCartesian(const Vector3& vec);
	std::ostream& operator<< (std::ostream& out, const Vector3& vec);

	const std::array<Vector2, 3> refSimplexNodes2D{ { {0.0,0.0}, {1.0,0.0}, {0.0,1.0} } }; // (0,0), (1,0), and (0,1)
	const std::array<Vector3, 4> refSimplexNodes3D{ { {0.0,0.0,0.0}, {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} } }; // (0,0,0), (1,0,0), (0,1,0), and (0,0,1)
}

#endif //_MESHHELPERFUNCTIONS_HPP_INCLUDED_