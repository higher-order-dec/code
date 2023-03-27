/*
ErrorAnalysisFunctions contains functions for testing and error analysis.
Author: Jonni Lohi, University of Jyväskylä, 2023.
*/

#ifndef _ERRORANALYSISFUNCTIONS_HPP_INCLUDED_
#define _ERRORANALYSISFUNCTIONS_HPP_INCLUDED_

#include "SmallSimplexPartition2D.hpp"
#include "SmallSimplexPartition3D.hpp"
#include "SmallCubePartition2D.hpp"
#include "SmallCubePartition3D.hpp"
#include "GFD/Mesh/BuilderMesh.hpp"
#include "GFD/Types/Types.hpp"

namespace gfd
{
	//examples
	
	void polynomialTestWhitney0Forms2D();
	void polynomialTestWhitney0Forms3D();
	void polynomialTestWhitney1Forms2D();
	void polynomialTestWhitney1Forms3D();
	void polynomialTestWhitney2Forms2D();
	void polynomialTestWhitney2Forms3D();
	void polynomialTestWhitney3Forms3D();
	void testWhitney1Forms3D();
	void testCubical1Forms2D();
	void testCubical1Forms3D();
	void hodgeTestWhitney1Forms2D();
	void hodgeTestWhitney1Forms3D();
	void laplacianTestWhitney2D();
	void laplacianTestWhitney3D();
	void minkowskiLaplacianTestWhitney2D();
	void minkowskiLaplacianTestWhitney3D();
	void polynomialPoissonExampleWhitney2D(const char testFunc, bool circumcentric = false);
	void poissonExampleWhitney2D(bool circumcentric = false);
	void polynomialPoissonExampleWhitney3D(const char testFunc, bool circumcentric = false);
	void poissonExampleWhitney3D(bool circumcentric = false);
	void polynomialWaveExampleWhitney2D(const char testFunc);
	void waveExampleWhitney2D();
	void polynomialWaveExampleWhitney3D(const char testFunc);
	void waveExampleWhitney3D();
	void waveExampleWhitney2D_TimeStepping(uint order, const double T, const uint T_REPEATS);
	void waveExampleWhitney2D_TimeSteppingFullSolution(uint order, const double T, const uint T_REPEATS);
	void waveExampleWhitney3D_TimeStepping(uint order, const double T, const uint T_REPEATS);
	void waveExampleWhitney3D_TimeSteppingFullSolution(uint order, const double T, const uint T_REPEATS);
	void waveExampleCubical2D_TimeStepping(uint order, const double T, const uint T_REPEATS);
	void waveExampleCubical3D_TimeStepping(uint order, const double T, const uint T_REPEATS);
	void waveExampleCubical3D_TimeSteppingFullSolution(uint order, const double T, const uint T_REPEATS);
	void waveExampleSmallSteps3D(uint order, const double T, const uint T_REPEATS);

	//analysis of stability and system matrix structure

	void poissonExampleWhitney2D_Stability(bool loadFromFile = false, bool circumcentric = false);
	void poissonExampleWhitney3D_Stability(bool loadFromFile = false, bool circumcentric = false);
	void waveExampleWhitney2D_Stability(uint order, const double T, const uint T_REPEATS);
	void waveExampleWhitney3D_Stability(uint order, const double T, const uint T_REPEATS);
	void waveExampleCubical2D_Stability(uint order, const double T, const uint T_REPEATS);
	void waveExampleCubical3D_Stability(uint order, const double T, const uint T_REPEATS);
	void waveExampleSmallSteps3D_Stability(uint order, const double T, const uint T_REPEATS);
	void waveExampleWhitney2D_PermutationTest(uint order, const uint meshNumber);
	void waveExampleWhitney3D_PermutationTest(uint order, const uint meshNumber);
	void waveExampleCubical3D_PermutationTest(uint order, const uint meshNumber);
	void waveExampleCubical2D_MatrixBlocks(uint order, const double T, const uint meshNumber);
	void waveExampleCubical3D_MatrixBlocks(uint order, const double T, const uint meshNumber);
	void waveExampleSmallSteps3D_MatrixBlocks(uint order, const double T, const uint meshNumber);

	//auxiliary functions

	double computeDifferenceL2Norm0Form(const Mesh& mesh_old, const SmallSimplexPartition2D& ssp, const Buffer<double>& discreteForm, const std::function<double(Vector2)>& fn,
		uint minElements = 0);
	double computeDifferenceL2Norm0Form(const Mesh& mesh_old, const SmallSimplexPartition3D& ssp, const Buffer<double>& discreteForm, const std::function<double(Vector3)>& fn,
		uint minElements = 0);
	double computeDifferenceL2Norm1Form(const Mesh& mesh_old, const SmallSimplexPartition2D& ssp, const Buffer<double>& discreteForm, const std::function<Vector2(Vector2)>& fn,
		uint minElements = 0);
	double computeDifferenceL2Norm1Form(const Mesh& mesh_old, const SmallSimplexPartition3D& ssp, const Buffer<double>& discreteForm, const std::function<Vector3(Vector3)>& fn,
		uint minElements = 0);
	double computeDifferenceL2Norm2Form(const Mesh& mesh_old, const SmallSimplexPartition2D& ssp, const Buffer<double>& discreteForm, const std::function<double(Vector2)>& fn,
		uint minElements = 0);
	double computeDifferenceL2Norm2Form(const Mesh& mesh_old, const SmallSimplexPartition3D& ssp, const Buffer<double>& discreteForm, const std::function<Vector3(Vector3)>& fn,
		uint minElements = 0);
	double computeDifferenceL2Norm3Form(const Mesh& mesh_old, const SmallSimplexPartition3D& ssp, const Buffer<double>& discreteForm, const std::function<double(Vector3)>& fn,
		uint minElements = 0);
	double computeDifferenceL2Norm1Form(const Mesh& mesh_old, const SmallCubePartition2D& scp, const Buffer<double>& discreteForm, const std::function<Vector2(Vector2)>& fn, double maxLen);
	double computeDifferenceL2Norm1Form(const Mesh& mesh_old, const SmallCubePartition3D& scp, const Buffer<double>& discreteForm, const std::function<Vector3(Vector3)>& fn, double maxLen);
	void solveSystem(const Buffer<std::unordered_map<uint, double>>& mat, Buffer<double>& x, const Buffer<double>& b, const Mesh& mesh, int rowFlag = -1, int colFlag = -1);
	void getDECSolutionPoissonEq(const BuilderMesh& mesh, const SmallSimplexPartition2D& ssp, const std::function<double(Vector2)>& sourceTerm, const std::function<double(Vector2)>& dirichletBC,
		const std::function<Vector2(Vector2)>& neumannBC, Buffer<double>& sol, bool circumcentric, bool oneElement = false);
	void getDECSolutionPoissonEq(const BuilderMesh& mesh, const SmallSimplexPartition3D& ssp, const std::function<double(Vector3)>& sourceTerm, const std::function<double(Vector3)>& dirichletBC,
		const std::function<Vector3(Vector3)>& neumannBC, Buffer<double>& sol, bool circumcentric, bool oneElement = false);
	void getDECSolutionWaveEq(const BuilderMesh& mesh_old, const BuilderMesh& mesh, const SmallSimplexPartition2D& ssp, const std::function<double(Vector2)>& sourceTerm,
		const std::function<double(Vector2)>& dirichletBC, const std::function<Vector2(Vector2)>& neumannBC, Buffer<double>& sol);
	void getDECSolutionWaveEq(const BuilderMesh& mesh_old, const BuilderMesh& mesh, const SmallSimplexPartition3D& ssp, const std::function<double(Vector3)>& sourceTerm,
		const std::function<double(Vector3)>& dirichletBC, const std::function<Vector3(Vector3)>& neumannBC, Buffer<double>& sol);
}

#endif //_ERRORANALYSISFUNCTIONS_HPP_INCLUDED_