/*
WhitneyForm.hpp provides functions for computations with Whitney forms.
Computations for affine-invariant quantities can be performed in a reference simplex.
Reference simplex has nodes (0,0), (1,0), (0,1) in 2D and (0,0,0), (1,0,0), (0,1,0), (0,0,1) in 3D (defined in MeshHelperFunctions.hpp).
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _WHITNEYFORM_HPP_INCLUDED_
#define _WHITNEYFORM_HPP_INCLUDED_

#include "MultiIndex.hpp"
#include "MeshHelperFunctions.hpp"
#include "GFD/Types/Types.hpp"
#include "GFD/Types/Vector.hpp"
#include "GFD/Types/Buffer.hpp"
#include "GFD/Mesh/DelaunayMesh.hpp"
#include <functional>

namespace gfd {
	void findBarycentricCoordinates(const Vector3& node0, const Vector3& node1, const Vector3& node2, const Vector3& node3, const Vector3& evaluationPoint,
		double& lambda0, double& lambda1, double& lambda2, double& lambda3);
	double evaluateBarycentricFunction(const Vector3& node0, const Vector3& node1, const Vector3& node2, const Vector3& node3, const Vector3& evaluationPoint, uint node);
	Vector2 evaluateWhitney1FormRefSimplex(double lambda0, double lambda1, double lambda2, uint whitneyForm);
	Vector3 evaluateWhitney1FormRefSimplex(double lambda0, double lambda1, double lambda2, double lambda3, uint whitneyForm);	
	Vector2 evaluateWhitney1FormRefSimplex(const Vector2& evaluationPoint, const mi_t& edge);
	Vector3 evaluateWhitney1FormRefSimplex(const Vector3& evaluationPoint, const mi_t& edge);
	Vector3 evaluateWhitney2FormRefSimplex(const Vector3& evaluationPoint, const mi_t& face);
	double evaluate0Form(const DelaunayMesh& mesh, const Buffer<double>& discreteForm, const Vector3& evaluationPoint);
	Vector3 evaluate1Form(const DelaunayMesh& mesh, const Buffer<double>& discreteForm, const Vector3& evaluationPoint);
	Vector3 evaluate2Form(const DelaunayMesh& mesh, const Buffer<double>& discreteForm, const Vector3& evaluationPoint);
	double evaluate3Form(const DelaunayMesh& mesh, const Buffer<double>& discreteForm, const Vector3& evaluationPoint);
	double evaluateBarycentricProduct(const DelaunayMesh& mesh, const mi_t& exponent, const Vector3& evaluationPoint);
	double evaluateBarycentricProductRefSimplex(const mi_t& exponent, const Vector2& evaluationPoint);
	double evaluateBarycentricProductRefSimplex(const mi_t& exponent, const Vector3& evaluationPoint);
	double integrateBarycentricProduct(const mi_t& exponent);
	double integrateBarycentricProduct(const mi_t& exponent, const mi_t& subsimplex);
}

#endif //_WHITNEYFORM_HPP_INCLUDED_