#include "WhitneyForm.hpp"

namespace gfd {

	//gradients of barycentric functions in the reference simplex
	const std::array<Vector2, 3> gradientsInRefSimplex2D{ { {-1.0,-1.0}, {1.0,0.0}, {0.0,1.0} } };
	const std::array<Vector3, 4> gradientsInRefSimplex3D{ { {-1.0,-1.0,-1.0}, {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} } };

	//finds the barycentric coordinates of evaluationPoint with respect to the given nodes
	void findBarycentricCoordinates(const Vector3& node0, const Vector3& node1, const Vector3& node2, const Vector3& node3, const Vector3& evaluationPoint,
		double& lambda0, double& lambda1, double& lambda2, double& lambda3) {
		double denominator = ThreeVector3(node1 - node0, node2 - node1, node3 - node2).determinant();
		lambda0 = std::abs(ThreeVector3(node1 - evaluationPoint, node2 - node1, node3 - node2).determinant() / denominator);
		lambda1 = std::abs(ThreeVector3(node0 - evaluationPoint, node2 - node0, node3 - node2).determinant() / denominator);
		lambda2 = std::abs(ThreeVector3(node0 - evaluationPoint, node1 - node0, node3 - node1).determinant() / denominator);
		lambda3 = std::abs(ThreeVector3(node0 - evaluationPoint, node1 - node0, node2 - node1).determinant() / denominator);
		//lambda3 = 1.0 - lambda0 - lambda1 - lambda2; (may be less precise)
	}

	//Evaluates the barycentric function corresponding to node at evaluationPoint. Barycentric functions are with respect to the given nodes.
	double evaluateBarycentricFunction(const Vector3& node0, const Vector3& node1, const Vector3& node2, const Vector3& node3, const Vector3& evaluationPoint, uint node) {
		switch (node) {
		case 0:
			return ThreeVector3(node1 - evaluationPoint, node2 - evaluationPoint, node3 - evaluationPoint).determinant() / ThreeVector3(node1 - node0, node2 - node0, node3 - node0).determinant();
		case 1:
			return ThreeVector3(node0 - evaluationPoint, node2 - evaluationPoint, node3 - evaluationPoint).determinant() / ThreeVector3(node0 - node1, node2 - node1, node3 - node1).determinant();
		case 2:
			return ThreeVector3(node0 - evaluationPoint, node1 - evaluationPoint, node3 - evaluationPoint).determinant() / ThreeVector3(node0 - node2, node1 - node2, node3 - node2).determinant();
		case 3:
			return ThreeVector3(node0 - evaluationPoint, node1 - evaluationPoint, node2 - evaluationPoint).determinant() / ThreeVector3(node0 - node3, node1 - node3, node2 - node3).determinant();
		}
	}

	//computes the value of the (lowest order) Whitney 1-form corresponding to the given edge at the point with given barycentric coordinates in 2D reference simplex
	Vector2 evaluateWhitney1FormRefSimplex(double lambda0, double lambda1, double lambda2, uint whitneyForm) {
		switch (whitneyForm) {
		case 0:
			return lambda0 * gradientsInRefSimplex2D[1] - lambda1 * gradientsInRefSimplex2D[0];
		case 1:
			return lambda0 * gradientsInRefSimplex2D[2] - lambda2 * gradientsInRefSimplex2D[0];
		case 2:
			return lambda1 * gradientsInRefSimplex2D[2] - lambda2 * gradientsInRefSimplex2D[1];
		default:
			return Vector2(0.0, 0.0);
		}
	}

	//computes the value of the (lowest order) Whitney 1-form corresponding to the given edge at the point with given barycentric coordinates in 3D reference simplex
	Vector3 evaluateWhitney1FormRefSimplex(double lambda0, double lambda1, double lambda2, double lambda3, uint whitneyForm) {
		switch (whitneyForm) {
		case 0:
			return lambda0 * gradientsInRefSimplex3D[1] - lambda1 * gradientsInRefSimplex3D[0];
		case 1:
			return lambda0 * gradientsInRefSimplex3D[2] - lambda2 * gradientsInRefSimplex3D[0];
		case 2:
			return lambda0 * gradientsInRefSimplex3D[3] - lambda3 * gradientsInRefSimplex3D[0];
		case 3:
			return lambda1 * gradientsInRefSimplex3D[2] - lambda2 * gradientsInRefSimplex3D[1];
		case 4:
			return lambda1 * gradientsInRefSimplex3D[3] - lambda3 * gradientsInRefSimplex3D[1];
		case 5:
			return lambda2 * gradientsInRefSimplex3D[3] - lambda3 * gradientsInRefSimplex3D[2];
		default:
			return Vector3(0, 0, 0);
		}
	}

	//computes the value of the (lowest order) Whitney 1-form corresponding to the given edge at evaluationPoint in 2D reference simplex
	Vector2 evaluateWhitney1FormRefSimplex(const Vector2& evaluationPoint, const mi_t& edge) {
		for (uint i = 0; i < faceEdges.size(); ++i)
			if (equals(edge, faceEdges[i]))
				return evaluateWhitney1FormRefSimplex(1.0 - evaluationPoint.x - evaluationPoint.y, evaluationPoint.x, evaluationPoint.y, i);
		return Vector2(0.0, 0.0);
	}

	//computes the value of the (lowest order) Whitney 1-form corresponding to the given edge at evaluationPoint in 3D reference simplex
	Vector3 evaluateWhitney1FormRefSimplex(const Vector3& evaluationPoint, const mi_t& edge) {
		for (uint i = 0; i < bodyEdges.size(); ++i)
			if (equals(edge, bodyEdges[i]))
				return evaluateWhitney1FormRefSimplex(1.0 - evaluationPoint.x - evaluationPoint.y - evaluationPoint.z, evaluationPoint.x, evaluationPoint.y, evaluationPoint.z, i);
		return Vector3(0.0, 0.0, 0.0);
	}

	//computes the value of the (lowest order) Whitney 2-form corresponding to the given face at evaluationPoint in 3D reference simplex
	Vector3 evaluateWhitney2FormRefSimplex(const Vector3& evaluationPoint, const mi_t& face) {
		uint oppositeNode = indexOf(face, 0);
		return (oppositeNode % 2 == 1 ? 1.0 : -1.0) * 2 * (refSimplexNodes3D[oppositeNode] - evaluationPoint);
	}

	//computes the value of the (lowest order) Whitney 0-form (interpolant of discreteForm) at evaluationPoint
	double evaluate0Form(const DelaunayMesh& mesh, const Buffer<double>& discreteForm, const Vector3& evaluationPoint) {

		//find the element containing evaluationPoint
		uint element = 0;
		if (!mesh.findElement(Vector4(evaluationPoint, 0.0), element))
			return 0.0;

		//find barycentric coordinates of evaluationPoint
		double lambda0, lambda1, lambda2, lambda3;
		const Buffer<uint> nodes = mesh.getBodyNodes(element);
		findBarycentricCoordinates(mesh.getNodePosition3(nodes[0]), mesh.getNodePosition3(nodes[1]), mesh.getNodePosition3(nodes[2]), mesh.getNodePosition3(nodes[3]),
			evaluationPoint, lambda0, lambda1, lambda2, lambda3);

		//compute the result
		return discreteForm[nodes[0]] * lambda0 + discreteForm[nodes[1]] * lambda1 + discreteForm[nodes[2]] * lambda2 + discreteForm[nodes[3]] * lambda3;
	}

	//computes the value of the (lowest order) Whitney 1-form (interpolant of discreteForm) at evaluationPoint
	Vector3 evaluate1Form(const DelaunayMesh& mesh, const Buffer<double>& discreteForm, const Vector3& evaluationPoint) {

		//find the element containing evaluationPoint
		uint element = 0;
		if (!mesh.findElement(Vector4(evaluationPoint, 0.0), element))
			return Vector3(0.0, 0.0, 0.0);

		//find the opposite faces for each node of the element
		const Buffer<uint> nodes = mesh.getBodyNodes(element);
		Buffer<uint> edges(6);
		edges[0] = mesh.findEdge(nodes[0], nodes[1]);
		edges[1] = mesh.findEdge(nodes[0], nodes[2]);
		edges[2] = mesh.findEdge(nodes[0], nodes[3]);
		edges[3] = mesh.findEdge(nodes[1], nodes[2]);
		edges[4] = mesh.findEdge(nodes[1], nodes[3]);
		edges[5] = mesh.findEdge(nodes[2], nodes[3]);
		Buffer<uint> faces(4);
		Buffer<uint> e(3);
		e[0] = edges[0];
		e[1] = edges[1];
		e[2] = edges[3];
		faces[0] = mesh.findFace(e); //face 0 opposite to node 3
		e[0] = edges[0];
		e[1] = edges[2];
		e[2] = edges[4];
		faces[1] = mesh.findFace(e); //face 1 opposite to node 2
		e[0] = edges[1];
		e[1] = edges[2];
		e[2] = edges[5];
		faces[2] = mesh.findFace(e); //face 2 opposite to node 1
		e[0] = edges[3];
		e[1] = edges[4];
		e[2] = edges[5];
		faces[3] = mesh.findFace(e); //face 3 opposite to node 0

		//find barycentric coordinates and their gradients at evaluationPoint
		double lambda0, lambda1, lambda2, lambda3;
		Vector4 gradLambda_0 = mesh.getFaceDeviation(faces[3], mesh.getNodePosition(nodes[0]));
		double lensq0 = gradLambda_0.lensq();
		lambda0 = mesh.getFaceDeviation(faces[3], Vector4(evaluationPoint, 0)).len() / std::sqrt(lensq0);
		gradLambda_0 /= lensq0;
		Vector4 gradLambda_1 = mesh.getFaceDeviation(faces[2], mesh.getNodePosition(nodes[1]));
		double lensq1 = gradLambda_1.lensq();
		lambda1 = mesh.getFaceDeviation(faces[2], Vector4(evaluationPoint, 0)).len() / std::sqrt(lensq1);
		gradLambda_1 /= lensq1;
		Vector4 gradLambda_2 = mesh.getFaceDeviation(faces[1], mesh.getNodePosition(nodes[2]));
		double lensq2 = gradLambda_2.lensq();
		lambda2 = mesh.getFaceDeviation(faces[1], Vector4(evaluationPoint, 0)).len() / std::sqrt(lensq2);
		gradLambda_2 /= lensq2;
		Vector4 gradLambda_3 = mesh.getFaceDeviation(faces[0], mesh.getNodePosition(nodes[3]));
		double lensq3 = gradLambda_3.lensq();
		lambda3 = mesh.getFaceDeviation(faces[0], Vector4(evaluationPoint, 0)).len() / std::sqrt(lensq3);
		gradLambda_3 /= lensq3;

		//compute the result
		return (mesh.getEdgeIncidence(edges[0], nodes[1]) * (lambda0 * gradLambda_1 - lambda1 * gradLambda_0) * discreteForm[edges[0]]
			+ mesh.getEdgeIncidence(edges[1], nodes[2]) * (lambda0 * gradLambda_2 - lambda2 * gradLambda_0) * discreteForm[edges[1]]
			+ mesh.getEdgeIncidence(edges[2], nodes[3]) * (lambda0 * gradLambda_3 - lambda3 * gradLambda_0) * discreteForm[edges[2]]
			+ mesh.getEdgeIncidence(edges[3], nodes[2]) * (lambda1 * gradLambda_2 - lambda2 * gradLambda_1) * discreteForm[edges[3]]
			+ mesh.getEdgeIncidence(edges[4], nodes[3]) * (lambda1 * gradLambda_3 - lambda3 * gradLambda_1) * discreteForm[edges[4]]
			+ mesh.getEdgeIncidence(edges[5], nodes[3]) * (lambda2 * gradLambda_3 - lambda3 * gradLambda_2) * discreteForm[edges[5]]).toVector3();
	}

	//computes the value of the (lowest order) Whitney 2-form (interpolant of discreteForm) at evaluationPoint
	Vector3 evaluate2Form(const DelaunayMesh& mesh, const Buffer<double>& discreteForm, const Vector3& evaluationPoint) {

		//find the element containing evaluationPoint
		uint element = 0;
		if (!mesh.findElement(Vector4(evaluationPoint, 0.0), element))
			return Vector3(0.0, 0.0, 0.0);

		//get the faces of the element in expected order and determine incidence numbers (here +1 if face normal vector points inward, -1 if outward)
		const Buffer<uint> nodes = mesh.getBodyNodes(element);
		Buffer<uint> faces(4);
		faces[0] = mesh.findFace(findEdges(nodes[0], nodes[1], nodes[2], mesh));
		faces[1] = mesh.findFace(findEdges(nodes[0], nodes[1], nodes[3], mesh));
		faces[2] = mesh.findFace(findEdges(nodes[0], nodes[2], nodes[3], mesh));
		faces[3] = mesh.findFace(findEdges(nodes[1], nodes[2], nodes[3], mesh));
		double incidence0 = (ThreeVector3(mesh.getFaceVector3(faces[0]), mesh.getNodePosition3(nodes[3]) - mesh.getNodePosition3(nodes[0])).xyz > 0 ? 1.0 : -1.0);
		double incidence1 = (ThreeVector3(mesh.getFaceVector3(faces[1]), mesh.getNodePosition3(nodes[2]) - mesh.getNodePosition3(nodes[0])).xyz > 0 ? 1.0 : -1.0);
		double incidence2 = (ThreeVector3(mesh.getFaceVector3(faces[2]), mesh.getNodePosition3(nodes[1]) - mesh.getNodePosition3(nodes[0])).xyz > 0 ? 1.0 : -1.0);
		double incidence3 = (ThreeVector3(mesh.getFaceVector3(faces[3]), mesh.getNodePosition3(nodes[0]) - mesh.getNodePosition3(nodes[1])).xyz > 0 ? 1.0 : -1.0);

		//compute the result
		std::vector<Vector3> pos{ {mesh.getNodePosition3(nodes[0]), mesh.getNodePosition3(nodes[1]), mesh.getNodePosition3(nodes[2]), mesh.getNodePosition3(nodes[3])} };
		double c = 1 / (0.5 * ThreeVector3(pos[1] - pos[0], pos[2] - pos[1], pos[3] - pos[2]).determinant());
		return std::abs(c) * (incidence0 * discreteForm[faces[0]] * (pos[3] - evaluationPoint)
			+ incidence1 * discreteForm[faces[1]] * (pos[2] - evaluationPoint)
			+ incidence2 * discreteForm[faces[2]] * (pos[1] - evaluationPoint)
			+ incidence3 * discreteForm[faces[3]] * (pos[0] - evaluationPoint));
	}

	/*
	Computes the value of the (lowest order) Whitney 3-form (interpolant of discreteForm) at evaluationPoint.
	N.B. discreteForm[i] is interpreted with the positive sign as if getBodyVector(i).xyz > 0
	*/
	double evaluate3Form(const DelaunayMesh& mesh, const Buffer<double>& discreteForm, const Vector3& evaluationPoint) {

		//find the element containing evaluationPoint
		uint element = 0;
		if (!mesh.findElement(Vector4(evaluationPoint, 0.0), element))
			return 0.0;

		//the Whitney 3-form corresponding to element is the reciprocal of the element volume
		return discreteForm[element] / std::abs(mesh.getBodyVector3(element).determinant());
	}

	//Computes the value of the product of lambda_i ^ exponent[i] at evaluationPoint.
	double evaluateBarycentricProduct(const DelaunayMesh& mesh, const mi_t& exponent, const Vector3& evaluationPoint) {

		//find the element containing evaluationPoint
		uint element = 0;
		if (!mesh.findElement(Vector4(evaluationPoint, 0.0), element))
			return 0.0;

		//find the barycentric coordinates at evaluationPoint
		std::vector<double> barycentricCoords(4);
		Buffer<uint> nodes = mesh.getBodyNodes(element);
		findBarycentricCoordinates(mesh.getNodePosition3(nodes[0]), mesh.getNodePosition3(nodes[1]), mesh.getNodePosition3(nodes[2]), mesh.getNodePosition3(nodes[3]),
			evaluationPoint, barycentricCoords[0], barycentricCoords[1], barycentricCoords[2], barycentricCoords[3]);

		//compute the result
		return pow(barycentricCoords, exponent);
	}

	//Returns the value of the product of lambda_i ^ exponent[i] at evaluationPoint of 2D reference simplex.
	double evaluateBarycentricProductRefSimplex(const mi_t& exponent, Vector2& evaluationPoint) {
		std::vector<double> barycentricCoords{ 1 - evaluationPoint.x - evaluationPoint.y, evaluationPoint.x, evaluationPoint.y };
		return pow(barycentricCoords, exponent);
	}

	//Returns the value of the product of lambda_i ^ exponent[i] at evaluationPoint of 3D reference simplex.
	double evaluateBarycentricProductRefSimplex(const mi_t& exponent, Vector3& evaluationPoint) {
		std::vector<double> barycentricCoords{ 1 - evaluationPoint.x - evaluationPoint.y - evaluationPoint.z, evaluationPoint.x, evaluationPoint.y, evaluationPoint.z};
		return pow(barycentricCoords, exponent);
	}

	//Returns the mean of the product of lambda_i ^ exponent[i] over the simplex of corresponding dimension.
	double integrateBarycentricProduct(const mi_t& exponent) {
		uint sum = gfd::sum(exponent);
		if (sum == 0) //we are integrating the constant function 1
			return 1.0;
		double numerator = factorial(exponent) * factorial(exponent.size() - 1);
		double denominator = factorial(exponent.size() - 1 + sum);
		return numerator / denominator;
	}

	//Returns the mean of the product of lambda_i ^ exponent[i] over the given subsimplex.
	double integrateBarycentricProduct(const mi_t& exponent, const mi_t& subsimplex) {
		mi_t exponentRestricted{};
		for (uint i = 0; i < subsimplex.size(); ++i) {
			if (subsimplex[i] == 1)
				exponentRestricted.push_back(exponent[i]);
			else if (exponent[i] != 0)
				return 0.0; //subsimplex[i] == 0 && exponent[i] != 0, so the integrand vanishes on subsimplex
		}
		return integrateBarycentricProduct(exponentRestricted);
	}
}