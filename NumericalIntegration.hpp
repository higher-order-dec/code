/*
NumericalIntegration.hpp provides functions for numerical integration on simplices.
Currently uses quadrature formulas that are
-29th order accurate for lines with 15 evaluation points, taken from https://pomax.github.io/bezierinfo/legendre-gauss.html
-11th order accurate for triangles with 28 evaluation points, taken from https://global-sci.org/intro/article_detail/jcm/8561.html
-11th order accurate for tetrahedra with 94 evaluation points, taken from https://doi.org/10.1002/NME.6313
Integration over more general cells can be done by dividing them into simplices.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _NUMERICALINTEGRATION_HPP_INCLUDED_
#define _NUMERICALINTEGRATION_HPP_INCLUDED_

#include "GFD/Types/Types.hpp"
#include "GFD/Types/Vector.hpp"
#include "GFD/Types/Buffer.hpp"
#include "GFD/Mesh/Mesh.hpp"
#include <functional>

namespace gfd {
	double integralAverage(const std::function<double(double)>& fn, double x0, double x1);
	double integralAverage(const std::function<double(Vector2)>& fn, const Vector2& node0, const Vector2& node1);
	double integralAverage(const std::function<double(Vector3)>& fn, const Vector3& node0, const Vector3& node1);
	double integralAverage(const std::function<double(Vector2)>& fn, const Vector2& node0, const Vector2& node1, const Vector2& node2);
	double integralAverage(const std::function<double(Vector3)>& fn, const Vector3& node0, const Vector3& node1, const Vector3& node2);
	double integralAverage(const std::function<double(Vector3)>& fn, const Vector3& node0, const Vector3& node1, const Vector3& node2, const Vector3& node3);
	void discretise0Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<double(Vector2)>& fn);
	void discretise0Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<double(Vector3)>& fn);
	void discretise0Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<double(Vector4)>& fn);
	void discretise1Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<Vector2(Vector2)>& fn);
	void discretise1Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<Vector3(Vector3)>& fn);
	void discretise1Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<Vector4(Vector4)>& fn);
	void discretise1FormDualMesh(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<Vector2(Vector2)>& fn);
	void discretise1FormWithPotential(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<double(Vector3)>& fn);
	void discretise2Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<double(Vector2)>& fn);
	void discretise2Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<Vector3(Vector3)>& fn);
	void discretise3Form(const Mesh& mesh, Buffer<double>& discreteForm, const std::function<double(Vector3)>& fn);
	double integrateDual1Cell(const Mesh& mesh, uint edge, const std::function<Vector2(Vector2)>& fn, bool circumcentric = false);
	double integrateDual2Cell(const Mesh& mesh, uint node, const std::function<double(Vector2)>& fn, bool circumcentric = false);
	double integrateDual2CellInParts(const Mesh& mesh, uint node, const std::function<double(Vector2)>& fn, double maxSize, bool circumcentric = false);
	double integrateDual2Cell(const Mesh& mesh, uint edge, const std::function<Vector3(Vector3)>& fn, bool circumcentric = false);
	double integrateDual2CellCartesianMesh(const Mesh& mesh, uint edge, const std::function<Vector3(Vector3)>& fn);
	double integrateDual3Cell(const Mesh& mesh, uint node, const std::function<double(Vector3)>& fn, bool circumcentric = false);
	double integrateDual3CellInParts(const Mesh& mesh, uint node, const std::function<double(Vector3)>& fn, double maxSize, bool circumcentric = false);
	double integrateDual3CellCartesianMesh(const Mesh& mesh, uint node, const std::function<double(Vector3)>& fn);
	double integrateNodeDualBoundary(const Mesh& mesh, uint node, const std::function<Vector2(Vector2)>& fn, bool circumcentric = false, bool minkowskiMetric = false);
	double integrateNodeDualBoundary(const Mesh& mesh, uint node, const std::function<Vector3(Vector3)>& fn, bool circumcentric = false, bool minkowskiMetric = false);
	double integrateNodeDualBoundaryCartesianMesh(const Mesh& mesh, uint node, const std::function<Vector3(Vector3)>& fn, bool minkowskiMetric);
	void integrateEdges(const Mesh& mesh, Buffer<double>& means, const std::function<double(Vector4)>& fn);
	void integrateTriangles(const Mesh& mesh, Buffer<double>& means, const std::function<double(Vector4)>& fn);
	double computeL2NormSquared(const Mesh& mesh, const std::function<double(Vector2)>& fn);
	double computeL2NormSquared(const Mesh& mesh, const std::function<double(Vector3)>& fn);
	double computeL2NormSquared(const Mesh& mesh, const std::function<Vector2(Vector2)>& fn);
	double computeL2NormSquared(const Mesh& mesh, const std::function<Vector3(Vector3)>& fn);
	double computeL2NormSquaredCartesianMesh(const Mesh& mesh, const std::function<Vector3(Vector3)>& fn);
	double triangleIntegral(const std::function<double(Vector2)>& fn, const Vector2& node0, const Vector2& node1, const Vector2& node2, const double maxSize);
	double tetrahedronIntegral(const std::function<double(Vector3)>& fn, const Vector3& node0, const Vector3& node1, const Vector3& node2, const Vector3& node3, const double maxSize);
	double parallelogramIntegral(const std::function<double(Vector2)>& fn, const Vector2& p0, const Vector2& p1, const Vector2& p2, const Vector2& p3, const double maxSize);
	double parallelepipedIntegral(const std::function<double(Vector3)>& fn, const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3,
		const Vector3& p4, const Vector3& p5, const Vector3& p6, const Vector3& p7, const double maxSize);
}
#endif //_NUMERICALINTEGRATION_HPP_INCLUDED_