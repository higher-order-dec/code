/**
 * BuilderMesh is a Mesh with several building routines.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _BUILDERMESH_HPP_INCLUDED_
#define _BUILDERMESH_HPP_INCLUDED_

#include "DelaunayMesh.hpp"
#include "../Types/Types.hpp"
#include "../Types/Buffer.hpp"
#include "../Types/Vector.hpp"
#include "../Types/UintSet.hpp"

namespace gfd
{

class BuilderMesh : public DelaunayMesh
{
public:
	BuilderMesh(const uint dim = 4);
	virtual ~BuilderMesh() { clear(); }
	void clear();

	void createCopy(const Mesh &mesh);

	// create regular grids
	void createGrid(const Vector4 &minp, const Vector4 &maxp, const double h);
	void createGrid(const Vector4 &minp, const Vector4 &maxp, const Vector4 &h);

	// create 2d grids
	void createTriangleGrid(const Vector2 &minp, const Vector2 &maxp, const double h, const bool rect = false);
	void createHexagonGrid(const Vector2 &minp, const Vector2 &maxp, const double h);
	void createSnubSquareGrid(const Vector2 &minp, const Vector2 &maxp, const double h);
	void createTetrilleGrid(const Vector2 &minp, const Vector2 &maxp, const double h);
	void createCircleBoundaryGrid(const double r, const double h);

	// create 3d grids
	void createFccGrid(const Vector3 &minp, const Vector3 &maxp, const double h);
	void createBccGrid(const Vector3 &minp, const Vector3 &maxp, const double h);
	void createTruncatedOctahedraGrid(const Vector3 &minp, const Vector3 &maxp, const double h); // dual of BCC
	void createA15Grid(const Vector3 &minp, const Vector3 &maxp, const double h);
	void createC15Grid(const Vector3 &minp, const Vector3 &maxp, const double h);
	void createZGrid(const Vector3 &minp, const Vector3 &maxp, const double h);
	void createSphereBoundaryGrid(const double r, const double h, const uint optimizationSteps = 0, const bool integerDivision = false);
	void createSphereBoundaryPentaGrid(const double r, const double h);
	void createFaceSplit(const Mesh &mesh, const uint div);

	// create 4d grids
	void createQccGrid(const Vector4 &minp, const Vector4 &maxp, const Vector4 &h);
	void createAsyncQccGrid(const Vector4 &minp, const Vector4 &maxp, const Vector4 &h);
	void createAsyncQccGrid2(const Vector4 &minp, const Vector4 &maxp, const Vector4 &h);

	// create modifications
	bool createRepeated(const Mesh &mesh, const Vector4 &pos, const Vector4 &step, const uint steps) { return createRepeated(mesh, pos, step, step, steps); }
	bool createRepeated(const Mesh &mesh, const Vector4 &pos, const Vector4 &dir, const Vector4 &step, const uint steps);
	void createIntersection(const Mesh &mesh, const Vector4 &v, const double dot); // create intersection mesh at plane x, where v.dot(x) = dot
	void createDualMesh(const Mesh &mesh);

	// flags
	void clearNodeFlags() { m_nflag.clear(); }
	void clearEdgeFlags() { m_eflag.clear(); }
	void clearFaceFlags() { m_fflag.clear(); }
	void clearBodyFlags() { m_bflag.clear(); }
	void clearQuadFlags() { m_qflag.clear(); }
	void clearFlags();
	void fillNodeFlags(const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillEdgeFlags(const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillFaceFlags(const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillBodyFlags(const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillQuadFlags(const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillFlags(const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillBoundaryFlags(const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillNodeRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillEdgeRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillFaceRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillBodyRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillQuadRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag = UINTSETALL);
	void fillRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag = UINTSETALL);
	void expandFlags(const uint flag, const UintSet &oldflag = UINTSETALL, uint layers = uint(-1));
	void removeByFlags(const UintSet &flag = UINTSETALL);

	// combine with other mesh
	void combine(const Mesh &mesh);

	// stretch mesh
	void stretchLinear(const Vector4 &v, const uint steps, const UintSet &flag = UINTSETALL, const uint flagMiddle = 0, const uint flagEnd = 0);

	// repeat mesh blocks
	bool repeatMiddle(const Vector4 &pos, const Vector4 &step, const uint steps) { return repeatMiddle(pos, step, step, steps); }
	bool repeatMiddle(const Vector4 &pos, const Vector4 &dir, const Vector4 &step, const uint steps);

	// hodge optimization
	void improveNodeByHodge(const uint n, const bool position = true, const bool weight = true);
	void optimizeNodes(const UintSet &flag = UINTSETALL, const uint iterations = 100, const bool position = true, const bool weight = true);
	bool optimizeNodesIteration(const Buffer<uint> &n, const bool position = true, const bool weight = true);

protected:

	// stretch mesh
	void stretch(const Buffer<uint> &n, const uint steps, const UintSet &flag, const uint flagMiddle, const uint flagEnd);
	uint getIndex(uint i, uint num, const Buffer<uint> &ind, const Buffer<uint> &slot, const Buffer<uint> &link, const Buffer<uint> &s);
	double getCellLengthSq(const uint i) const;

};

}

#endif //_BUILDERMESH_HPP_INCLUDED_
