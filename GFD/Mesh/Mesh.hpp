/**
 * Mesh is a class for partitioning a 1--4 dimensional domain with polyhedral cells.
 * We use the following naming: Node = 0-cell, Edge = 1-cell, Face = 2-cell, Body = 3-cell, and Quad = 4-cell.
 * Cells are constructed recursively and are linked with their boundary cells and parent cells.
 * Each cell can be assigned with a flag (unsigned int)
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _MESH_HPP_INCLUDED_
#define _MESH_HPP_INCLUDED_

#include "../Types/Types.hpp"
#include "../Types/Buffer.hpp"
#include "../Types/Text.hpp"
#include "../Types/Vector.hpp"
#include "../Types/UintSet.hpp"
#include <fstream>

using namespace std;

namespace gfd
{

struct Node
{
	Buffer<uint> e;
};

struct Edge
{
	Buffer<uint> n;
	Buffer<uint> f;
};

struct Face
{
	Buffer<uint> e;
	Buffer<uint> b;
};

struct Body
{
	Buffer<uint> f;
	Buffer<uint> q;
};

struct Quad
{
	Buffer<uint> b;
};

class Mesh
{
public:
	Mesh(const uint dim = 4);
	virtual ~Mesh() { clear(); }
	virtual void clear();
	void swap(Mesh &mesh);

	// load and save mesh
	//bool loadMesh(const std::string &path);
	//bool saveMesh(const std::string &path) const;
	bool loadJRMesh(const std::string &path);
	bool saveJRMesh(const std::string &path) const;

	// statistics
	void writeStatistics(Text &text, const UintSet &flag = UINTSETALL) const;

	// vector space dimension
	uint getDimension() const { return m_dim; }

	// number of cells
	uint getNodeSize() const { return m_nsize; }
	uint getEdgeSize() const { return m_esize; }
	uint getFaceSize() const { return m_fsize; }
	uint getBodySize() const { return m_bsize; }
	uint getQuadSize() const { return m_qsize; }

	// intersections
	uint getEdgeIntersection(const uint e0, const uint e1) const { return m_e[e0].n.getFirstIntersection(m_e[e1].n, NONE); }
	uint getFaceIntersection(const uint f0, const uint f1) const { return m_f[f0].e.getFirstIntersection(m_f[f1].e, NONE); }
	uint getBodyIntersection(const uint b0, const uint b1) const { return m_b[b0].f.getFirstIntersection(m_b[b1].f, NONE); }
	uint getQuadIntersection(const uint q0, const uint q1) const { return m_q[q0].b.getFirstIntersection(m_q[q1].b, NONE); }

	// neighbors
	const Buffer<uint> &getNodeEdges(const uint n) const { return m_n[n].e; }
	Buffer<uint> getNodeFaces(const uint n) const;
	Buffer<uint> getNodeBodies(const uint n) const;
	Buffer<uint> getNodeQuads(const uint n) const;
	const Buffer<uint> &getEdgeNodes(const uint e) const { return m_e[e].n; }
	const Buffer<uint> &getEdgeFaces(const uint e) const { return m_e[e].f; }
	Buffer<uint> getEdgeBodies(const uint e) const;
	Buffer<uint> getEdgeQuads(const uint e) const;
	Buffer<uint> getFaceNodes(const uint f) const;
	const Buffer<uint> &getFaceEdges(const uint f) const { return m_f[f].e; }
	const Buffer<uint> &getFaceBodies(const uint f) const { return m_f[f].b; }
	Buffer<uint> getFaceQuads(const uint f) const;
	Buffer<uint> getBodyNodes(const uint b) const;
	Buffer<uint> getBodyEdges(const uint b) const;
	const Buffer<uint> &getBodyFaces(const uint b) const { return m_b[b].f; }
	const Buffer<uint> &getBodyQuads(const uint b) const { return m_b[b].q; }
	Buffer<uint> getQuadNodes(const uint q) const;
	Buffer<uint> getQuadEdges(const uint q) const;
	Buffer<uint> getQuadFaces(const uint q) const;
	const Buffer<uint> &getQuadBodies(const uint q) const { return m_q[q].b; }

	// get single neighbor
	uint getNodeAnyEdge(const uint n) const { return m_n[n].e[0]; }
	uint getNodeAnyFace(const uint n) const { return m_e[m_n[n].e[0]].f[0]; }
	uint getNodeAnyBody(const uint n) const { return m_f[m_e[m_n[n].e[0]].f[0]].b[0]; }
	uint getNodeAnyQuad(const uint n) const { return m_b[m_f[m_e[m_n[n].e[0]].f[0]].b[0]].q[0]; }
	uint getEdgeAnyNode(const uint e) const { return m_e[e].n[0]; }
	uint getEdgeAnyFace(const uint e) const { return m_e[e].f[0]; }
	uint getEdgeAnyBody(const uint e) const { return m_f[m_e[e].f[0]].b[0]; }
	uint getEdgeAnyQuad(const uint e) const { return m_b[m_f[m_e[e].f[0]].b[0]].q[0]; }
	uint getFaceAnyNode(const uint f) const { return m_e[m_f[f].e[0]].n[0]; }
	uint getFaceAnyEdge(const uint f) const { return m_f[f].e[0]; }
	uint getFaceAnyBody(const uint f) const { return m_f[f].b[0]; }
	uint getFaceAnyQuad(const uint f) const { return m_b[m_f[f].b[0]].q[0]; }
	uint getBodyAnyNode(const uint b) const { return m_e[m_f[m_b[b].f[0]].e[0]].n[0]; }
	uint getBodyAnyEdge(const uint b) const { return m_f[m_b[b].f[0]].e[0]; }
	uint getBodyAnyFace(const uint b) const { return m_b[b].f[0]; }
	uint getBodyAnyQuad(const uint b) const { return m_b[b].q[0]; }
	uint getQuadAnyNode(const uint q) const { return m_e[m_f[m_b[m_q[q].b[0]].f[0]].e[0]].n[0]; }
	uint getQuadAnyEdge(const uint q) const { return m_f[m_b[m_q[q].b[0]].f[0]].e[0]; }
	uint getQuadAnyFace(const uint q) const { return m_b[m_q[q].b[0]].f[0]; }
	uint getQuadAnyBody(const uint q) const { return m_q[q].b[0]; }

	// more neighbors
	uint getEdgeOtherNode(const uint e, const uint n) const;

	// incidence
	sign getEdgeIncidence(const uint e, const uint n) const;
	sign getFaceIncidence(const uint f, const uint e) const;
	sign getBodyIncidence(const uint b, const uint f) const;
	sign getQuadIncidence(const uint q, const uint b) const;

	// circumcenter positions
	double getNodePosition1(const uint n) const { return m_p[n * m_dim]; }
	Vector2 getNodePosition2(const uint n) const { const uint i = n * m_dim; return Vector2(m_p[i], m_p[i+1]); }
	Vector3 getNodePosition3(const uint n) const { const uint i = n * m_dim; return Vector3(m_p[i], m_p[i+1], m_p[i+2]); }
	Vector4 getNodePosition4(const uint n) const { const uint i = n * m_dim; return Vector4(m_p[i], m_p[i+1], m_p[i+2], m_p[i+3]); }
	Vector4 getNodePosition(const uint n) const;
	double getEdgePosition1(const uint e) const;
	Vector2 getEdgePosition2(const uint e) const;
	Vector3 getEdgePosition3(const uint e) const;
	Vector4 getEdgePosition4(const uint e) const;
	Vector4 getEdgePosition(const uint e) const;
	Vector2 getFacePosition2(const uint f) const;
	Vector3 getFacePosition3(const uint f) const;
	Vector4 getFacePosition4(const uint f) const;
	Vector4 getFacePosition(const uint f) const;
	Vector3 getBodyPosition3(const uint b) const;
	Vector4 getBodyPosition4(const uint b) const;
	Vector4 getBodyPosition(const uint b) const;
	Vector4 getQuadPosition(const uint q) const;

	// average positions
	double getEdgeAverage1(const uint e) const;
	Vector2 getEdgeAverage2(const uint e) const;
	Vector3 getEdgeAverage3(const uint e) const;
	Vector4 getEdgeAverage4(const uint e) const;
	Vector4 getEdgeAverage(const uint e) const;
	Vector2 getFaceAverage2(const uint f) const;
	Vector3 getFaceAverage3(const uint f) const;
	Vector4 getFaceAverage4(const uint f) const;
	Vector4 getFaceAverage(const uint f) const;
	Vector3 getBodyAverage3(const uint b) const;
	Vector4 getBodyAverage4(const uint b) const;
	Vector4 getBodyAverage(const uint b) const;
	Vector4 getQuadAverage(const uint q) const;

	// primal volume vectors
	double getEdgeVector1(const uint e) const;
	Vector2 getEdgeVector2(const uint e) const;
	Vector3 getEdgeVector3(const uint e) const;
	Vector4 getEdgeVector4(const uint e) const;
	Vector4 getEdgeVector(const uint e) const;
	TwoVector2 getFaceVector2(const uint f) const;
	TwoVector3 getFaceVector3(const uint f) const;
	TwoVector4 getFaceVector4(const uint f) const;
	TwoVector4 getFaceVector(const uint f) const;
	ThreeVector3 getBodyVector3(const uint b) const;
	ThreeVector4 getBodyVector4(const uint b) const;
	ThreeVector4 getBodyVector(const uint b) const;
	FourVector4 getQuadVector(const uint q) const;

	// circumcenter determinants
	double getNodeWeight(const uint n) const { if(n < m_w.size()) return m_w[n]; return 0.0; }
	SymMatrix4 getMetric(const uint flag = 0) const;
	uint getMetricSize() const { return 2 * m_m.size() / (m_dim * (m_dim + 1)); }
	Vector4 getTransformed(const Vector4 &r, const uint flag = 0) const;

	// dual average positions
	double getNodeDualAverage1(const uint n) const;
	Vector2 getNodeDualAverage2(const uint n) const;
	Vector3 getNodeDualAverage3(const uint n) const;
	Vector4 getNodeDualAverage4(const uint n) const;
	Vector4 getNodeDualAverage(const uint n) const;
	Vector2 getEdgeDualAverage2(const uint e) const;
	Vector3 getEdgeDualAverage3(const uint e) const;
	Vector4 getEdgeDualAverage4(const uint e) const;
	Vector4 getEdgeDualAverage(const uint e) const;
	Vector3 getFaceDualAverage3(const uint f) const;
	Vector4 getFaceDualAverage4(const uint f) const;
	Vector4 getFaceDualAverage(const uint f) const;
	Vector4 getBodyDualAverage4(const uint b) const;
	Vector4 getBodyDualAverage(const uint b) const;

	// dual volume vectors
	double getNodeDualVector1(const uint n) const;
	TwoVector2 getNodeDualVector2(const uint n) const;
	ThreeVector3 getNodeDualVector3(const uint n) const;
	FourVector4 getNodeDualVector4(const uint n) const;
	FourVector4 getNodeDualVector(const uint n) const;
	double getEdgeDualVector1(const uint e) const;
	Vector2 getEdgeDualVector2(const uint e) const;
	TwoVector3 getEdgeDualVector3(const uint e) const;
	ThreeVector4 getEdgeDualVector4(const uint e) const;
	ThreeVector4 getEdgeDualVector(const uint e) const;
	double getFaceDualVector2(const uint f) const;
	Vector3 getFaceDualVector3(const uint f) const;
	TwoVector4 getFaceDualVector4(const uint f) const;
	TwoVector4 getFaceDualVector(const uint f) const;
	double getBodyDualVector3(const uint b) const;
	Vector4 getBodyDualVector4(const uint b) const;
	Vector4 getBodyDualVector(const uint b) const;
	double getQuadDualVector(const uint q) const;

	// diagonal Hodge terms
	double getNodeHodge(const uint n) const;
	double getNodeHodge(const uint n, const double &metric) const;
	double getEdgeHodge(const uint e) const;
	double getEdgeHodge(const uint e, const SymMatrix4 &metric) const;
	double getFaceHodge(const uint f) const;
	double getFaceHodge(const uint f, const SymTwoMatrix4 &metric) const;
	double getBodyHodge(const uint b) const;
	double getBodyHodge(const uint b, const SymThreeMatrix4 &metric) const;
	double getQuadHodge(const uint q) const;
	double getQuadHodge(const uint q, const SymFourMatrix4 &metric) const;

	// get deviation vector from the cell plane to position p (independent of metric)
	Vector4 getEdgeDeviation(const uint e, const Vector4 &p) const;
	Vector4 getFaceDeviation(const uint f, const Vector4 &p) const;
	Vector4 getBodyDeviation(const uint b, const Vector4 &p) const;

	// get projection of vector d that is orthogonal to the cell (depend on the metric)
	//Vector4 getEdgeOrthogonal(const uint e, const Vector4 &d) const;
	//Vector4 getFaceOrthogonal(const uint f, const Vector4 &d) const;
	//Vector4 getBodyOrthogonal(const uint b, const Vector4 &d) const;

	// find elements
	uint findNode(const Vector4 &p, const double zerolensq = 1e-13, uint curr = 0, const bool assured = true) const;
	uint findEdge(const uint n0, const uint n1) const;
	uint findFace(const Buffer<uint> &e) const;
	uint findBody(const Buffer<uint> &f) const;
	uint findQuad(const Buffer<uint> &b) const;

	// neighbors for advanced use only
	void setNodeEdges(const uint n, const Buffer<uint> &e) { m_n[n].e = e; }
	void setEdgeNodes(const uint e, const Buffer<uint> &n) { m_e[e].n = n; }
	void setEdgeFaces(const uint e, const Buffer<uint> &f) { m_e[e].f = f; }
	void setFaceEdges(const uint f, const Buffer<uint> &e) { m_f[f].e = e; }
	void setFaceBodies(const uint f, const Buffer<uint> &b) { m_f[f].b = b; }
	void setBodyFaces(const uint b, const Buffer<uint> &f) { m_b[b].f = f; }
	void setBodyQuads(const uint b, const Buffer<uint> &q) { m_b[b].q = q; }
	void setQuadBodies(const uint q, const Buffer<uint> &b) { m_q[q].b = b; }

	// flags
	uint getNodeFlag(const uint i) const { if(i < m_nflag.size()) return m_nflag[i]; return 0; }
	uint getEdgeFlag(const uint i) const { if(i < m_eflag.size()) return m_eflag[i]; return 0; }
	uint getFaceFlag(const uint i) const { if(i < m_fflag.size()) return m_fflag[i]; return 0; }
	uint getBodyFlag(const uint i) const { if(i < m_bflag.size()) return m_bflag[i]; return 0; }
	uint getQuadFlag(const uint i) const { if(i < m_qflag.size()) return m_qflag[i]; return 0; }

	// add and remove cells
	uint addNode(const Vector4 &p);
	uint addEdge(const uint n0, const uint n1);
	uint addFace(const Buffer<uint> &e);
	uint addBody(const Buffer<uint> &f);
	uint addQuad(const Buffer<uint> &b);
	void removeNode(const uint n);
	void removeEdge(const uint e);
	void removeFace(const uint f);
	void removeBody(const uint b);
	void removeQuad(const uint q);

	// transformation and relocation
	void setNodePosition(const uint n, const Vector4 &p);
	void transform(const Matrix4 &mat);
	void move(const Vector4 &vec);
	void setNodeWeight(const uint n, const double w);
	void setMetric(const SymMatrix4 &m, const uint flag = 0);

	// modify flags
	void setNodeFlag(const uint n, const uint flag);
	void setEdgeFlag(const uint e, const uint flag);
	void setFaceFlag(const uint f, const uint flag);
	void setBodyFlag(const uint b, const uint flag);
	void setQuadFlag(const uint q, const uint flag);

	// resizing the element buffers (use this only for optimization, if you know the element sizes in advance)
	void resizeNodeBuffer(const uint size);
	void resizeEdgeBuffer(const uint size);
	void resizeFaceBuffer(const uint size);
	void resizeBodyBuffer(const uint size);
	void resizeQuadBuffer(const uint size);

	Buffer<Vector4> &getEdgeSimplices(const uint e, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getFaceSimplices(const uint f, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getBodySimplices(const uint b, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getQuadSimplices(const uint q, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getNodeEdgeSimplices(const uint n, Buffer<uint> &e, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getNodeFaceSimplices(const uint n, Buffer<uint> &f, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getNodeBodySimplices(const uint n, Buffer<uint> &b, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getNodeQuadSimplices(const uint n, Buffer<uint> &q, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getEdgeNodeSimplices(const uint e, Buffer<uint> &n, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getEdgeFaceSimplices(const uint e, Buffer<uint> &f, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getEdgeBodySimplices(const uint e, Buffer<uint> &b, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getEdgeQuadSimplices(const uint e, Buffer<uint> &q, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getFaceNodeSimplices(const uint f, Buffer<uint> &n, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getFaceEdgeSimplices(const uint f, Buffer<uint> &e, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getFaceBodySimplices(const uint f, Buffer<uint> &b, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getFaceQuadSimplices(const uint f, Buffer<uint> &q, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getBodyNodeSimplices(const uint b, Buffer<uint> &n, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getBodyEdgeSimplices(const uint b, Buffer<uint> &e, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getBodyFaceSimplices(const uint b, Buffer<uint> &f, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getBodyQuadSimplices(const uint b, Buffer<uint> &q, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getQuadNodeSimplices(const uint q, Buffer<uint> &n, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getQuadEdgeSimplices(const uint q, Buffer<uint> &e, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getQuadFaceSimplices(const uint q, Buffer<uint> &f, Buffer<Vector4> &p) const;
	Buffer<Vector4> &getQuadBodySimplices(const uint q, Buffer<uint> &b, Buffer<Vector4> &p) const;

	Vector4 getNodeEdgeVectors(const uint n, Buffer<uint> &e, Buffer<Vector4> &v) const;
	Vector4 getNodeFaceVectors(const uint n, Buffer<uint> &f, Buffer<TwoVector4> &v) const;
	Vector4 getNodeBodyVectors(const uint n, Buffer<uint> &b, Buffer<ThreeVector4> &v) const;
	Vector4 getNodeQuadVectors(const uint n, Buffer<uint> &q, Buffer<FourVector4> &v) const;
	Vector4 getEdgeNodeVectors(const uint e, Buffer<uint> &n, Buffer<Vector4> &v) const;
	Vector4 getEdgeFaceVectors(const uint e, Buffer<uint> &f, Buffer<Vector4> &v) const;
	Vector4 getEdgeBodyVectors(const uint e, Buffer<uint> &b, Buffer<TwoVector4> &v) const;
	Vector4 getEdgeQuadVectors(const uint e, Buffer<uint> &q, Buffer<ThreeVector4> &v) const;
	Vector4 getFaceNodeVectors(const uint f, Buffer<uint> &n, Buffer<TwoVector4> &v) const;
	Vector4 getFaceEdgeVectors(const uint f, Buffer<uint> &e, Buffer<Vector4> &v) const;
	Vector4 getFaceBodyVectors(const uint f, Buffer<uint> &b, Buffer<Vector4> &v) const;
	Vector4 getFaceQuadVectors(const uint f, Buffer<uint> &q, Buffer<TwoVector4> &v) const;
	Vector4 getBodyNodeVectors(const uint b, Buffer<uint> &n, Buffer<ThreeVector4> &v) const;
	Vector4 getBodyEdgeVectors(const uint b, Buffer<uint> &e, Buffer<TwoVector4> &v) const;
	Vector4 getBodyFaceVectors(const uint b, Buffer<uint> &f, Buffer<Vector4> &v) const;
	Vector4 getBodyQuadVectors(const uint b, Buffer<uint> &q, Buffer<Vector4> &v) const;
	Vector4 getQuadNodeVectors(const uint q, Buffer<uint> &n, Buffer<FourVector4> &v) const;
	Vector4 getQuadEdgeVectors(const uint q, Buffer<uint> &e, Buffer<ThreeVector4> &v) const;
	Vector4 getQuadFaceVectors(const uint q, Buffer<uint> &f, Buffer<TwoVector4> &v) const;
	Vector4 getQuadBodyVectors(const uint q, Buffer<uint> &b, Buffer<Vector4> &v) const;

protected:

	uint m_dim; // vector space dimension
	Buffer<double> m_p; // vector coordinates for each node

	// mesh elements
	uint m_nsize;
	Buffer<Node> m_n; // nodes
	uint m_esize;
	Buffer<Edge> m_e; // edges
	uint m_fsize;
	Buffer<Face> m_f; // faces
	uint m_bsize;
	Buffer<Body> m_b; // bodies
	uint m_qsize;
	Buffer<Quad> m_q; // quads

	// flags
	Buffer<uint> m_nflag; // node flags (optional)
	Buffer<uint> m_eflag; // edge flags (optional)
	Buffer<uint> m_fflag; // face flags (optional)
	Buffer<uint> m_bflag; // body flags (optional)
	Buffer<uint> m_qflag; // quad flags (optional)

	// circumcenter computation (squared distance of v is v.dot(m * v) + w)
	Buffer<double> m_m; // symmetric matrices to determine dot product (optional)
	Buffer<double> m_w; // node weights (optional)

	virtual const string getJRMeshType() const { return "JRM1"; }
	virtual bool loadJRMeshMore(std::ifstream &fs) { return !fs.fail(); }
	virtual bool saveJRMeshMore(std::ofstream &fs) const { return !fs.fail(); }
	Vector2 getCellPosition2(const Buffer<uint> &n, const uint flag, const SymMatrix2 &a0) const;
	Vector3 getCellPosition3(const Buffer<uint> &n, const uint flag, const SymMatrix3 &a0) const;
	Vector4 getCellPosition4(const Buffer<uint> &n, const uint flag, const SymMatrix4 &a0) const;
	template<typename T> void setTerm(const uint i, const T &term, const T &zero, const uint bufs, Buffer<T> &buf) const;
	double getMetric1(const uint flag) const;
	SymMatrix2 getMetric2(const uint flag) const;
	SymMatrix3 getMetric3(const uint flag) const;
	SymMatrix4 getMetric4(const uint flag) const;
	double getTransformed1(const double &r, const uint flag) const;
	Vector2 getTransformed2(const Vector2 &r, const uint flag) const;
	Vector3 getTransformed3(const Vector3 &r, const uint flag) const;
	Vector4 getTransformed4(const Vector4 &r, const uint flag) const;

public:
	void orderFaceEdges(const uint f);
	void orderBodyFaces(const uint b);
	void orderQuadBodies(const uint q);
};

}

#endif //_MESH_HPP_INCLUDED_
