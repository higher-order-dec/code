/**
 * DelaunayMesh includes insert and erase routines of nodes for convex weighted Delaunay mesh.
 * The metric should be positive definite and constant through the domain.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 * Modified by Jonni Lohi, University of Jyväskylä, 2023: added the member function "findElement".
 */

#ifndef _DELAUNAYMESH_HPP_INCLUDED_
#define _DELAUNAYMESH_HPP_INCLUDED_

#include "Mesh.hpp"
#include "../Types/Types.hpp"
#include "../Types/Buffer.hpp"
#include "../Types/Vector.hpp"

namespace gfd
{

class DelaunayMesh : public Mesh
{
public:
	DelaunayMesh(const uint dim = 4) : Mesh(dim) {}
	virtual ~DelaunayMesh() { clear(); }

	bool includesCell(const Vector4 &p, uint &curr) const; // searches the nearest cell. returns true if the cell (of id curr) includes p. otherwise returns false and curr is the id of a boundary cell.
	uint searchNode(const Vector4 &p, uint curr = 0) const; // returns the id of the node corresponding to the Voronoi-cell including p
	bool findElement(const Vector4& p, uint& elementIndex) const; //Finds the index of the element containing p in 2D or 3D mesh. Similar to includesCell but more reliable.

	// insert and erase weighted Delaunay nodes
	uint insertNode(const Vector4 &p, const double w, uint near, const bool forced);
	//uint insertNode(const Vector4 &p, double w = 0.0, const ushort s = 0, uint near = 0, const bool forced = false);
	bool eraseNode(const uint n);

	// insert Delaunay mesh
	bool insertMesh(const DelaunayMesh &mesh);

	// metric
	double getRadiusSq(const Vector4 &p, const uint node) const;
	double getRadiusSq(const Vector4 &p, const Buffer<uint> &n) const;


protected:

	struct CellSet
	{
		// tables of detached cells
		uint ns;
		Buffer<uint> n;
		uint es;
		Buffer<uint> e;
		uint fs;
		Buffer<uint> f;
		uint bs;
		Buffer<uint> b;
		uint qs;
		Buffer<uint> q;
		CellSet() { ns = es = fs = bs = qs = 0; }
	};

	// detach and attach cells (Use detach to remove without changing existing ids. Finalize with removeDetached.)
	uint attachNode(const Vector4 &p, CellSet &detached);
	uint attachEdge(const uint n0, const uint n1, CellSet &detached);
	uint attachFace(const Buffer<uint> &e, CellSet &detached);
	uint attachBody(const Buffer<uint> &f, CellSet &detached);
	uint attachQuad(const Buffer<uint> &b, CellSet &detached);
	void detachNode(const uint n, CellSet &detached);
	void detachEdge(const uint e, CellSet &detached);
	void detachFace(const uint f, CellSet &detached);
	void detachBody(const uint b, CellSet &detached);
	void detachQuad(const uint q, CellSet &detached);
	void detachEdgeRecursive(const uint e, CellSet &detached);
	void detachFaceRecursive(const uint f, CellSet &detached);
	void detachBodyRecursive(const uint b, CellSet &detached);
	void detachQuadRecursive(const uint q, CellSet &detached);
	void removeDetached(CellSet &detached);

	// building scripts
	void increaseDimension();
	uint mergeFace(const uint f, CellSet &detached);
	uint mergeBody(const uint b, CellSet &detached);
	uint mergeQuad(const uint q, CellSet &detached);

	// for orthogonal cells
	Vector4 getNodeAverage(const Buffer<uint> &n) const;
	Vector4 getEdgeDirection(const uint e, const uint n) const { return getEdgeAverage(e) - getNodePosition(n); }
	Vector4 getFaceDirection(const uint f, const uint e) const { return getNodeAverage(getFaceNodes(f)) - getNodePosition(getEdgeAnyNode(e)); }
	Vector4 getBodyDirection(const uint b, const uint f) const { return getNodeAverage(getBodyNodes(b)) - getNodePosition(getFaceAnyNode(f)); }
	Vector4 getQuadDirection(const uint q, const uint b) const { return getNodeAverage(getQuadNodes(q)) - getNodePosition(getBodyAnyNode(b)); }

	// metric
	double getLengthSq(const Vector4 &p, const uint node) const;
	bool isInsideSphere(const Vector4 &p, const double sq, const uint node) const;
	bool isInsideSphere(const Vector4 &p, const double sq, const Buffer<uint> &n) const;
	bool isOutsideSphere(const Vector4 &p, const double sq, const Buffer<uint> &n) const;
	uint eraseNodesInside(const Vector4 &p, const double sq, uint near);

};

}

#endif //_DELAUNAYMESH_HPP_INCLUDED_
