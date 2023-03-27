#include "SmallSimplexPartition.hpp"
#include "NumericalIntegration.hpp"
#include "WhitneyForm.hpp"
#include <iostream>

namespace gfd {

	SmallSimplexPartition::SmallSimplexPartition(uint ord) : order{ord == 0 ? 1 : ord} {} //1st order is the same as lowest order

	//returns the product of lambda_i ^ multiIndex[i] in the 2D reference simplex as 2D function
	std::function<double(Vector2)> SmallSimplexPartition::getBarycentricProductRefSimplex2D(const mi_t& multiIndex) {
		return [&](Vector2 p) -> double {
			return std::pow(1 - p.x - p.y, multiIndex[0]) * std::pow(p.x, multiIndex[1]) * std::pow(p.y, multiIndex[2]);
		};
	}

	//returns the product of lambda_i ^ multiIndex[i] in the 3D reference simplex as 3D function
	std::function<double(Vector3)> SmallSimplexPartition::getBarycentricProductRefSimplex3D(const mi_t& multiIndex) {
		return [&](Vector3 p) -> double {
			return std::pow(1 - p.x - p.y - p.z, multiIndex[0]) * (std::pow(p.x, multiIndex[1])) * std::pow(p.y, multiIndex[2]) * std::pow(p.z, multiIndex[3]);
		};
	}

	/*
	Computes the integral of (lowest order) Whitney form over a small simplex.
	The simplex of the Whitney form are given in parameter whitneyForm.
	The small simplex corresponds to multiIndex and nodeIndices.
	All parameters should have the same size.
	*/
	double SmallSimplexPartition::integrateWhitneyForm(const mi_t& multiIndex, const mi_t& nodeIndices, const mi_t& whitneyForm) {
		double order = sum(multiIndex);
		Buffer<uint> subsimplexNodes;
		nonzeroIndices(nodeIndices, subsimplexNodes);
		if (subsimplexNodes.size() == 1) { //evaluation of 0-form at 0-simplex
			uint whitneyFormIndex = gfd::indexOf(whitneyForm, 1);
			return (multiIndex[whitneyFormIndex] + nodeIndices[whitneyFormIndex]) / (order + 1.0);
		}
		else if (multiIndex.size() == 2) { //1-simplex in 1D
			return 1.0 / (order + 1.0);
		}
		else if (multiIndex.size() == 3) {
			if (subsimplexNodes.size() == 2) { //1-simplex in 2D
				Vector2 subsimplexVec = (refSimplexNodes2D[subsimplexNodes[1]] - refSimplexNodes2D[subsimplexNodes[0]]) / (order + 1.0);
				//choose the image of subsimplexNodes[0] as evaluationPoint
				Vector2 evaluationPoint = getNodeImageRefSimplex2D(multiIndex, subsimplexNodes[0], order);
				//compute the value of whitneyForm at evaluationPoint
				Vector2 whitneyFormValue = evaluateWhitney1FormRefSimplex(evaluationPoint, whitneyForm);
				return whitneyFormValue.dot(subsimplexVec);
			}
			else { //2-simplex in 2D
				return 1.0 / std::pow(order + 1, 2);
			}
		}
		else { //assume multiIndex.size() == 4
			if (subsimplexNodes.size() == 2) { //1-simplex in 3D
				Vector3 subsimplexVec = (refSimplexNodes3D[subsimplexNodes[1]] - refSimplexNodes3D[subsimplexNodes[0]]) / (order + 1.0);
				//choose the image of subsimplexNodes[0] as evaluationPoint
				Vector3 evaluationPoint = getNodeImageRefSimplex3D(multiIndex, subsimplexNodes[0], order);
				//compute the value of whitneyForm at evaluationPoint
				Vector3 whitneyFormValue = evaluateWhitney1FormRefSimplex(evaluationPoint, whitneyForm);
				return whitneyFormValue.dot(subsimplexVec);
			}
			else if (subsimplexNodes.size() == 3) { //2-simplex in 3D
				TwoVector3 subsimplexVec = TwoVector3(refSimplexNodes3D[subsimplexNodes[1]] - refSimplexNodes3D[subsimplexNodes[0]],
					refSimplexNodes3D[subsimplexNodes[2]] - refSimplexNodes3D[subsimplexNodes[0]]) / (2.0 * std::pow(order + 1, 2));
				//choose the image of subsimplexNodes[0] as evaluationPoint
				Vector3 evaluationPoint = getNodeImageRefSimplex3D(multiIndex, subsimplexNodes[0], order);
				//compute the value of whitneyForm at evaluationPoint
				Vector3 whitneyFormValue = evaluateWhitney2FormRefSimplex(evaluationPoint, whitneyForm);
				return whitneyFormValue.dot(subsimplexVec.dual());
			}
			else { //3-simplex in 3D
				return 1.0 / std::pow(order + 1, 3);
			}
		}
	}

	//Computes the integral of the (higher order) Whitney form (parameters miIntegrand and niIntegrand) over the small simplex (parameters miSimplex and niSimplex).
	double SmallSimplexPartition::integrateWhitneyForm(const mi_t& miSimplex, const mi_t& niSimplex, const mi_t& miIntegrand, const mi_t& niIntegrand) {
		return integralAverage(miSimplex, niSimplex, miIntegrand) * integrateWhitneyForm(miSimplex, niSimplex, niIntegrand);
	}

	/*
	Computes the integral of the (higher order) Whitney form (parameters miIntegrand and niIntegrand) over the small simplex (parameters miSimplex and niSimplex)
	using numerical integration rules for simplices to compute the mean of the barycentric product. Computations are done in the reference simplex.
	*/
	double SmallSimplexPartition::integrateWhitneyFormNumerically(const mi_t& miSimplex, const mi_t& niSimplex, const mi_t& miIntegrand, const mi_t& niIntegrand) {
		Buffer<uint> simplexNodes;
		nonzeroIndices(niSimplex, simplexNodes);
		double order = sum(miSimplex);
		uint i;
		if (simplexNodes.size() == 1) { //evaluation of 0-form at 0-simplex
			double res = 1.0;
			for (i = 0; i < miSimplex.size(); ++i) {
				res *= std::pow((miSimplex[i] + niSimplex[i]) / (order + 1.0), miIntegrand[i] + niIntegrand[i]);
			}
			return res;
		}
		else if (miSimplex.size() == 2) { //1-simplex in 1D
			double x0 = miSimplex[1] / (order + 1.0);
			double x1 = (miSimplex[1] + 1) / (order + 1.0);
			std::function<double(double)> fn{ //the barycentric product
			[&](double x) -> double {
			 return std::pow(1 - x, miIntegrand[0]) * std::pow(x, miIntegrand[1]);
			}
			};
			return gfd::integralAverage(fn, x0, x1) / (order + 1.0);
		}
		else if (miSimplex.size() == 3) { //2D
			if (simplexNodes.size() == 2) { //1-simplex in 2D
				Vector2 node0 = getNodeImageRefSimplex2D(miSimplex, simplexNodes[0], order);
				Vector2 node1 = getNodeImageRefSimplex2D(miSimplex, simplexNodes[1], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex2D(miIntegrand), node0, node1) * integrateWhitneyForm(miSimplex, niSimplex, niIntegrand);
			}
			else { //2-simplex in 2D
				Vector2 node0 = getNodeImageRefSimplex2D(miSimplex, simplexNodes[0], order);
				Vector2 node1 = getNodeImageRefSimplex2D(miSimplex, simplexNodes[1], order);
				Vector2 node2 = getNodeImageRefSimplex2D(miSimplex, simplexNodes[2], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex2D(miIntegrand), node0, node1, node2) / std::pow(order + 1, 2);
			}
		}
		else { //3D
			if (simplexNodes.size() == 2) { //1-simplex in 3D
				Vector3 node0 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[0], order);
				Vector3 node1 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[1], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex3D(miIntegrand), node0, node1) * integrateWhitneyForm(miSimplex, niSimplex, niIntegrand);
			}
			else if (simplexNodes.size() == 3) { //2-simplex in 3D
				Vector3 node0 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[0], order);
				Vector3 node1 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[1], order);
				Vector3 node2 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[2], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex3D(miIntegrand), node0, node1, node2) * integrateWhitneyForm(miSimplex, niSimplex, niIntegrand);
			}
			else { //3-simplex in 3D
				Vector3 node0 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[0], order);
				Vector3 node1 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[1], order);
				Vector3 node2 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[2], order);
				Vector3 node3 = getNodeImageRefSimplex3D(miSimplex, simplexNodes[3], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex3D(miIntegrand), node0, node1, node2, node3) / std::pow(order + 1, 3);
			}
		}
	}

	//Computes the mean of barycentricProduct over smallSimplex. The small simplex is the image of the corresponding big simplex.
	double SmallSimplexPartition::integralAverage(const mi_t& smallSimplex, const mi_t& barycentricProduct) {
		if (smallSimplex.size() == 1)
			return 1.0;
		uint orderSimplex = sum(smallSimplex);
		if (orderSimplex == 0)
			return integrateBarycentricProduct(barycentricProduct); //simpler but not necessary, as the general case would give the correct result anyway
		std::vector<mi_t> multiIndices;
		getSmallerMultiIndices(barycentricProduct, multiIndices);
		uint orderProduct = sum(barycentricProduct);
		double multiplier = 1.0 / std::pow(1 + orderSimplex, orderProduct);
		double sum = 0;
		for (uint i = 0; i < multiIndices.size(); ++i) {
			const mi_t& r = multiIndices[i];
			sum += binCoef(barycentricProduct, r) * pow(smallSimplex, subtract(barycentricProduct, r)) * integrateBarycentricProduct(r);
		}
		return multiplier * sum;
	}

	//Computes the mean of barycentricProduct over smallSimplex. The small simplex is the image of the subsimplex given in nodeIndices.
	double SmallSimplexPartition::integralAverage(const mi_t& multiIndex, const mi_t& nodeIndices, const mi_t& barycentricProduct) {
		double order = sum(multiIndex);
		mi_t multiIndexRestricted;
		mi_t barycentricProductRestricted;
		double constantPart = 1.0; //this part of the product stays constant on the small simplex in question
		double scalingFactor = 1.0; //this is the sum of the barycentric functions on the small simplex in question
		for (uint i = 0; i < multiIndex.size(); ++i) {
			if (nodeIndices[i] == 1) {
				multiIndexRestricted.push_back(multiIndex[i]);
				barycentricProductRestricted.push_back(barycentricProduct[i]);
			}
			else {
				double lambda_i = multiIndex[i] / (order + 1.0);
				constantPart *= std::pow(lambda_i, barycentricProduct[i]);
				scalingFactor -= lambda_i;
			}
		}
		uint orderIntegrand = sum(barycentricProductRestricted);
		if (orderIntegrand == 0)
			return constantPart; //the whole product is constant on the given small simplex
		return constantPart * std::pow(scalingFactor, orderIntegrand) * SmallSimplexPartition::integralAverage(multiIndexRestricted, barycentricProductRestricted);

		//Alternatively the following code could be used:
		/*double multiplier = std::pow(1.0 / (sum(multiIndex) + 1.0), sum(barycentricProduct)); // 1/(k'+1)^k
		mi_t multiIndexRestricted;
		mi_t barycentricProductRestricted;
		double constantPart = 1.0; //this part of the product stays constant on the small simplex in question
		for (uint i = 0; i < multiIndex.size(); ++i) {
			if (nodeIndices[i] == 1) {
				multiIndexRestricted.push_back(multiIndex[i]);
				barycentricProductRestricted.push_back(barycentricProduct[i]);
			}
			else {
				constantPart *= std::pow(multiIndex[i], barycentricProduct[i]);
			}
		}

		if (multiIndexRestricted.size() == 1) //mean of 0-form over 0-simplex is simply the value
			return multiplier * constantPart * std::pow(multiIndexRestricted[0] + 1, barycentricProductRestricted[0]);

		std::vector<mi_t> multiIndices;
		getSmallerMultiIndices(barycentricProductRestricted, multiIndices);
		double sum = 0;
		for (uint i = 0; i < multiIndices.size(); ++i) {
			const mi_t& r = multiIndices[i];
			sum += binCoef(barycentricProductRestricted, r) * pow(multiIndexRestricted, subtract(barycentricProductRestricted, r)) * integrateBarycentricProduct(r);
		}
		return multiplier * constantPart * sum;*/
	}

	/*
	Computes the mean of barycentricProduct over smallSimplex. The small simplex is the image of the subsimplex given in nodeIndices.
	Computations are done in the reference simplex using numerical integration rules for simplices.
	*/
	double SmallSimplexPartition::integralAverageNumerically(const mi_t& multiIndex, const mi_t& nodeIndices, const mi_t& barycentricProduct) {
		Buffer<uint> simplexNodes;
		nonzeroIndices(nodeIndices, simplexNodes);
		double order = sum(multiIndex);
		uint i;
		if (simplexNodes.size() == 1) { //evaluation of 0-form at 0-simplex
			double res = 1.0;
			for (i = 0; i < multiIndex.size(); ++i) {
				res *= std::pow((multiIndex[i] + nodeIndices[i]) / (order + 1.0), barycentricProduct[i]);
			}
			return res;
		}
		else if (multiIndex.size() == 2) { //1-simplex in 1D
			double x0 = multiIndex[1] / (order + 1.0);
			double x1 = (multiIndex[1] + 1) / (order + 1.0);
			std::function<double(double)> fn{ //the barycentric product
			[&](double x) -> double {
			 return std::pow(1 - x, barycentricProduct[0]) * std::pow(x, barycentricProduct[1]);
			}
			};
			return gfd::integralAverage(fn, x0, x1);
		}
		else if (multiIndex.size() == 3) { //2D
			if (simplexNodes.size() == 2) { //1-simplex in 2D
				Vector2 node0 = getNodeImageRefSimplex2D(multiIndex, simplexNodes[0], order);
				Vector2 node1 = getNodeImageRefSimplex2D(multiIndex, simplexNodes[1], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex2D(barycentricProduct), node0, node1);
			}
			else { //2-simplex in 2D
				Vector2 node0 = getNodeImageRefSimplex2D(multiIndex, simplexNodes[0], order);
				Vector2 node1 = getNodeImageRefSimplex2D(multiIndex, simplexNodes[1], order);
				Vector2 node2 = getNodeImageRefSimplex2D(multiIndex, simplexNodes[2], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex2D(barycentricProduct), node0, node1, node2);
			}
		}
		else if (multiIndex.size() == 4) { //3D
			if (simplexNodes.size() == 2) { //1-simplex in 3D
				Vector3 node0 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[0], order);
				Vector3 node1 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[1], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex3D(barycentricProduct), node0, node1);
			}
			else if (simplexNodes.size() == 3) { //2-simplex in 3D
				Vector3 node0 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[0], order);
				Vector3 node1 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[1], order);
				Vector3 node2 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[2], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex3D(barycentricProduct), node0, node1, node2);
			}
			else { //3-simplex in 3D
				Vector3 node0 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[0], order);
				Vector3 node1 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[1], order);
				Vector3 node2 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[2], order);
				Vector3 node3 = getNodeImageRefSimplex3D(multiIndex, simplexNodes[3], order);
				return gfd::integralAverage(getBarycentricProductRefSimplex3D(barycentricProduct), node0, node1, node2, node3);
			}
		}
	}

	//returns true if multiIndex maps subsimplex to the interior of the big simplex
	bool SmallSimplexPartition::mapsToInterior(const mi_t& multiIndex, const mi_t& subsimplex) {
		for (uint i = 0; i < multiIndex.size(); ++i) {
			if (multiIndex[i] == 0 && subsimplex[i] == 0)
				return false;
		}
		return true;
	}

	//returns the image of multiIndex when we map the given node in 2D reference simplex
	Vector2 SmallSimplexPartition::getNodeImageRefSimplex2D(const mi_t& multiIndex, uint node, uint sum) {
		if (sum == 0)
			sum = gfd::sum(multiIndex);
		Vector2 image(multiIndex[1], multiIndex[2]);
		if (node == 1)
			image.x += 1.0;
		else if (node == 2)
			image.y += 1.0;
		return image / (sum + 1.0);
	}

	//returns the image of multiIndex when we map the given node in 3D reference simplex
	Vector3 SmallSimplexPartition::getNodeImageRefSimplex3D(const mi_t& multiIndex, uint node, uint sum) {
		if (sum == 0)
			sum = gfd::sum(multiIndex);
		Vector3 image(multiIndex[1], multiIndex[2], multiIndex[3]);
		if (node == 1)
			image.x += 1.0;
		else if (node == 2)
			image.y += 1.0;
		else if (node == 3)
			image.z += 1.0;
		return image / (sum + 1.0);
	}

	//returns the image of multiIndex when we map the given node
	Vector4 SmallSimplexPartition::getNodeImage(const mi_t& multiIndex, const Mesh& mesh, const Buffer<uint>& nodes, const mi_t& node) const {
		Vector4 image(0, 0, 0, 0);
		for (uint i = 0; i < multiIndex.size(); ++i) {
			image += (node[i] + multiIndex[i]) * mesh.getNodePosition(nodes[i]);
		}
		return image / order;
	}

	//returns the index of the multiIndex in multiIndexList
	uint SmallSimplexPartition::indexOf(const Buffer<mi_t>& multiIndexList, const mi_t& multiIndex) {
		for (uint i = 0; i < multiIndexList.size(); ++i) {
			if (equals(multiIndexList[i], multiIndex))
				return i;
		}
		return NONE;
	}

	//returns the index of the first small simplex with the given multiIndex in smallSimplexList
	uint SmallSimplexPartition::indexOf(const Buffer<SmallSimplex>& smallSimplexList, const mi_t& multiIndex) {
		for (uint i = 0; i < smallSimplexList.size(); ++i) {
			if (equals(*smallSimplexList[i].multiIndex, multiIndex))
				return i;
		}
		return NONE;
	}

	//returns the index of the small simplex in smallSimplexList by comparing pointers (use with caution)
	uint SmallSimplexPartition::indexOf(const Buffer<SmallSimplex>& smallSimplexList, const SmallSimplex& smallSimplex) {
		for (uint i = 0; i < smallSimplexList.size(); ++i) {
			if (smallSimplexList[i].multiIndex == smallSimplex.multiIndex && smallSimplexList[i].nodeIndices == smallSimplex.nodeIndices)
				return i;
		}
		return NONE;
	}

	//prints the small simplices in ssList to the console
	void SmallSimplexPartition::printSmallSimplexList(const Buffer<SmallSimplex>& ssList) {
		for (uint i = 0; i < ssList.size(); ++i) {
			std::cout << "Multi-index: ";
			printMultiIndex(*ssList[i].multiIndex);
			std::cout << ". Node indices: ";
			printMultiIndex(*ssList[i].nodeIndices);
			std::cout << ".\n";
		}
	}
}