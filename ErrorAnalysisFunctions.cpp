#include "ErrorAnalysisFunctions.hpp"
#include "LinearAlgebraFunctions.hpp"
#include "MeshHelperFunctions.hpp"
#include "NumericalIntegration.hpp"

namespace gfd {

	//test Whitney 0-forms of different orders with polynomial test functions
	void polynomialTestWhitney0Forms2D() {

		std::function<double(Vector2)> testFunction1{ //testFunction1(x,y) = 0.25
					[](Vector2 p) -> double {
						return 0.25;
					}
		};
		std::function<double(Vector2)> testFunction2{ //testFunction2(x,y) = 1/56 * x^4y - 1/45 * y^5
					[](Vector2 p) -> double {
						return p.x * p.x * p.x * p.x * p.y / 56 - p.y * p.y * p.y * p.y * p.y / 45;
					}
		};
		std::function<double(Vector2)> testFunction3{ //testFunction3(x,y,z) = 1/4800 * x^9y - 1/1300 * y^10
					[](Vector2 p) -> double {
						return p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y / 4800 - p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y / 1300;
					}
		};
		std::cout.precision(1);
		std::cout << std::scientific;
		Text output;
		output.precision(1);
		output << std::scientific;

		//create a simplicial mesh
		const double r = 2.0;
		BuilderMesh mesh(2);
		mesh.createTriangleGrid(Vector2(-r, -r), Vector2(r, r), 2.0); //mesh grain is 2.0

		for (uint k = 1; k <= 12; ++k) {
			//refine the mesh into small simplices of order k
			uint order = k;
			BuilderMesh mesh_new(2);
			SmallSimplexPartition2D ssp(order);
			ssp.refineMesh(mesh, mesh_new);

			//discretise testFunctions
			Buffer<double> discreteForm1(mesh_new.getNodeSize());
			discretise0Form(mesh_new, discreteForm1, testFunction1);
			Buffer<double> discreteForm2(mesh_new.getNodeSize());
			discretise0Form(mesh_new, discreteForm2, testFunction2);
			Buffer<double> discreteForm3(mesh_new.getNodeSize());
			discretise0Form(mesh_new, discreteForm3, testFunction3);

			//compute L2 norm of the differences of testFunctions and the interpolants of discreteForms
			std::cout << "Order: " << order << '\n';
			output << "Order: " << order << '\n';
			double diff1 = computeDifferenceL2Norm0Form(mesh, ssp, discreteForm1, testFunction1);
			std::cout << "testFunction1, error L2 norm: " << diff1 << '\n';
			output << "testFunction1, error L2 norm: " << diff1 << '\n';
			double diff2 = computeDifferenceL2Norm0Form(mesh, ssp, discreteForm2, testFunction2);
			std::cout << "testFunction2, error L2 norm: " << diff2 << '\n';
			output << "testFunction2, error L2 norm: " << diff2 << '\n';
			double diff3 = computeDifferenceL2Norm0Form(mesh, ssp, discreteForm3, testFunction3);
			std::cout << "testFunction3, error L2 norm: " << diff3 << '\n';
			output << "testFunction3, error L2 norm: " << diff3 << '\n';
			output.save("Files/polynomialTestWhitney0Forms2D.txt");
		}
	}

	//test Whitney 0-forms of different orders with polynomial test functions
	void polynomialTestWhitney0Forms3D() {
		std::function<double(Vector3)> testFunction1{ //testFunction1(x,y,z) = 0.25
					[](Vector3 p) -> double {
						return 0.25;
					}
		};
		std::function<double(Vector3)> testFunction2{ //testFunction2(x,y,z) = 64/75 * x^2y^2z - 8/75 * z^5
					[](Vector3 p) -> double {
						return (64.0 / 75.0) * p.x * p.x * p.y * p.y * p.z - (8.0 / 75.0) * p.z * p.z * p.z * p.z * p.z;
					}
		};
		std::function<double(Vector3)> testFunction3{ //testFunction3(x,y,z) = 32/11 * x^4y^4z^2 - 1/176 * z^10
					[](Vector3 p) -> double {
						return  (32.0 / 11.0) * p.x * p.x * p.x * p.x * p.y * p.y * p.y * p.y * p.z * p.z - p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z / 176.0;
					}
		};
		std::cout.precision(1);
		std::cout << std::scientific;
		Text output;
		output.precision(1);
		output << std::scientific;

		//create a simplicial mesh in Omega
		const double r = 2.0125;
		BuilderMesh mesh(3);
		mesh.createBccGrid(Vector3(-r, -r, -r), Vector3(r, r, r), 2.0); //mesh grain is 2.0
		for (uint i = mesh.getNodeSize(); i-- > 0; )
		{
			const Vector4 p = mesh.getNodePosition(i);
			if (p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.x > r - 1e-5) mesh.removeNode(i);
		}

		for (uint k = 1; k <= 12; ++k) {
			//refine the mesh into small simplices of order k
			uint order = k;
			BuilderMesh mesh_new(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh, mesh_new);

			//discretise testFunctions
			Buffer<double> discreteForm1(mesh_new.getNodeSize());
			discretise0Form(mesh_new, discreteForm1, testFunction1);
			Buffer<double> discreteForm2(mesh_new.getNodeSize());
			discretise0Form(mesh_new, discreteForm2, testFunction2);
			Buffer<double> discreteForm3(mesh_new.getNodeSize());
			discretise0Form(mesh_new, discreteForm3, testFunction3);

			//compute L2 norm of the differences of testFunctions and the interpolants of discreteForms
			std::cout << "Order: " << order << '\n';
			output << "Order: " << order << '\n';
			double diff1 = computeDifferenceL2Norm0Form(mesh, ssp, discreteForm1, testFunction1);
			std::cout << "testFunction1, error L2 norm: " << diff1 << '\n';
			output << "testFunction1, error L2 norm: " << diff1 << '\n';
			double diff2 = computeDifferenceL2Norm0Form(mesh, ssp, discreteForm2, testFunction2);
			std::cout << "testFunction2, error L2 norm: " << diff2 << '\n';
			output << "testFunction2, error L2 norm: " << diff2 << '\n';
			double diff3 = computeDifferenceL2Norm0Form(mesh, ssp, discreteForm3, testFunction3);
			std::cout << "testFunction3, error L2 norm: " << diff3 << '\n';
			output << "testFunction3, error L2 norm: " << diff3 << '\n';
			output.save("Files/polynomialTestWhitney0Forms3D.txt");
		}
	}

	//test Whitney 1-forms of different orders with polynomial test functions
	void polynomialTestWhitney1Forms2D() {
		std::function<Vector2(Vector2)> testFunction1{ //testFunction1(x,y) = 3/14 dx - 1/14 dy
					[](Vector2 p) -> Vector2 {
						return Vector2(3.0 / 14.0, -1.0 / 14.0);
					}
		};
		std::function<Vector2(Vector2)> testFunction2{ //testFunction2(x,y) = 1/55 * (x^4y dx + xy^4 dy)
					[](Vector2 p) -> Vector2 {
						return (1.0 / 55.0) * Vector2(p.x * p.x * p.x * p.x * p.y, p.x * p.y * p.y * p.y * p.y);
					}
		};
		std::function<Vector2(Vector2)> testFunction3{ //testFunction3(x,y) = 1/3400 * (x^9y dx + xy^9 dy)
					[](Vector2 p) -> Vector2 {
						return (1.0 / 3400.0) * Vector2(p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y, p.x * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y);
					}
		};
		std::cout.precision(1);
		std::cout << std::scientific;
		Text output;
		output.precision(1);
		output << std::scientific;

		//create a simplicial mesh
		const double r = 2.0;
		BuilderMesh mesh(2);
		mesh.createTriangleGrid(Vector2(-r, -r), Vector2(r, r), 2.0); //mesh grain is 2.0

		for (uint k = 1; k <= 12; ++k) {
			//refine the mesh into small simplices of order k
			uint order = k;
			BuilderMesh mesh_new(2);
			SmallSimplexPartition2D ssp(order);
			ssp.refineMesh(mesh, mesh_new);

			//discretise testFunctions
			Buffer<double> discreteForm1(mesh_new.getEdgeSize());
			discretise1Form(mesh_new, discreteForm1, testFunction1);
			Buffer<double> discreteForm2(mesh_new.getEdgeSize());
			discretise1Form(mesh_new, discreteForm2, testFunction2);
			Buffer<double> discreteForm3(mesh_new.getEdgeSize());
			discretise1Form(mesh_new, discreteForm3, testFunction3);

			//compute L2 norm of the differences of testFunctions and the interpolants of discreteForms
			std::cout << "Order: " << order << '\n';
			output << "Order: " << order << '\n';
			double diff1 = computeDifferenceL2Norm1Form(mesh, ssp, discreteForm1, testFunction1);
			std::cout << "testFunction1, error L2 norm: " << diff1 << '\n';
			output << "testFunction1, error L2 norm: " << diff1 << '\n';
			double diff2 = computeDifferenceL2Norm1Form(mesh, ssp, discreteForm2, testFunction2);
			std::cout << "testFunction2, error L2 norm: " << diff2 << '\n';
			output << "testFunction2, error L2 norm: " << diff2 << '\n';
			double diff3 = computeDifferenceL2Norm1Form(mesh, ssp, discreteForm3, testFunction3);
			std::cout << "testFunction3, error L2 norm: " << diff3 << '\n';
			output << "testFunction3, error L2 norm: " << diff3 << '\n';
			output.save("Files/polynomialTestWhitney1Forms2D.txt");
		}
	}

	//test Whitney 1-forms of different orders with polynomial test functions
	void polynomialTestWhitney1Forms3D() {
		std::function<Vector3(Vector3)> testFunction1{ //testFunction1(x,y,z) = 30/128 dx - 10/128 dy + 10/252 dz
					[](Vector3 p) -> Vector3 {
						return Vector3(30.0 / 128.0, -10.0 / 128.0, 10.0 / 252.0);
					}
		};
		std::function<Vector3(Vector3)> testFunction2{ //testFunction2(x,y,z) = x^2y^2z dx + x^2yz^2 dy + xy^2z^2 dz
					[](Vector3 p) -> Vector3 {
						return Vector3(p.x * p.x * p.y * p.y * p.z, p.x * p.x * p.y * p.z * p.z, p.x * p.y * p.y * p.z * p.z);
					}
		};
		std::function<Vector3(Vector3)> testFunction3{ //testFunction3(x,y,z) = 20/9 * (x^2y^4z^4 dx + x^4y^2z^4 dy + x^4y^4z^2 dz)
					[](Vector3 p) -> Vector3 {
						return (20.0 / 9.0) * Vector3(p.x * p.x * p.y * p.y * p.y * p.y * p.z * p.z * p.z * p.z, p.x * p.x * p.x * p.x * p.y * p.y * p.z * p.z * p.z * p.z,
							p.x * p.x * p.x * p.x * p.y * p.y * p.y * p.y * p.z * p.z);
					}
		};
		std::cout.precision(1);
		std::cout << std::scientific;
		Text output;
		output.precision(1);
		output << std::scientific;

		//create a simplicial mesh in Omega
		const double r = 2.0125;
		BuilderMesh mesh(3);
		mesh.createBccGrid(Vector3(-r, -r, -r), Vector3(r, r, r), 2.0); //mesh grain is 2.0
		for (uint i = mesh.getNodeSize(); i-- > 0; )
		{
			const Vector4 p = mesh.getNodePosition(i);
			if (p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.x > r - 1e-5) mesh.removeNode(i);
		}

		for (uint k = 1; k <= 12; ++k) {
			//refine the mesh into small simplices of order k
			uint order = k;
			BuilderMesh mesh_new(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh, mesh_new);

			//discretise testFunctions
			Buffer<double> discreteForm1(mesh_new.getEdgeSize());
			discretise1Form(mesh_new, discreteForm1, testFunction1);
			Buffer<double> discreteForm2(mesh_new.getEdgeSize());
			discretise1Form(mesh_new, discreteForm2, testFunction2);
			Buffer<double> discreteForm3(mesh_new.getEdgeSize());
			discretise1Form(mesh_new, discreteForm3, testFunction3);

			//compute L2 norm of the differences of testFunctions and the interpolants of discreteForms
			std::cout << "Order: " << order << '\n';
			output << "Order: " << order << '\n';
			double diff1 = computeDifferenceL2Norm1Form(mesh, ssp, discreteForm1, testFunction1);
			std::cout << "testFunction1, error L2 norm: " << diff1 << '\n';
			output << "testFunction1, error L2 norm: " << diff1 << '\n';
			double diff2 = computeDifferenceL2Norm1Form(mesh, ssp, discreteForm2, testFunction2);
			std::cout << "testFunction2, error L2 norm: " << diff2 << '\n';
			output << "testFunction2, error L2 norm: " << diff2 << '\n';
			double diff3 = computeDifferenceL2Norm1Form(mesh, ssp, discreteForm3, testFunction3);
			std::cout << "testFunction3, error L2 norm: " << diff3 << '\n';
			output << "testFunction3, error L2 norm: " << diff3 << '\n';
			output.save("Files/polynomialTestWhitney1Forms3D.txt");
		}
	}

	//test Whitney 2-forms of different orders with polynomial test functions
	void polynomialTestWhitney2Forms2D() {
		std::function<double(Vector2)> testFunction1{ //testFunction1(x,y) = 0.25
					[](Vector2 p) -> double {
						return 0.25;
					}
		};
		std::function<double(Vector2)> testFunction2{ //testFunction2(x,y) = 1/56 * x^4y - 1/45 * y^5
					[](Vector2 p) -> double {
						return p.x * p.x * p.x * p.x * p.y / 56 - p.y * p.y * p.y * p.y * p.y / 45;
					}
		};
		std::function<double(Vector2)> testFunction3{ //testFunction3(x,y,z) = 1/4800 * x^9y - 1/1300 * y^10
					[](Vector2 p) -> double {
						return p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y / 4800 - p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y / 1300;
					}
		};
		std::cout.precision(1);
		std::cout << std::scientific;
		Text output;
		output.precision(1);
		output << std::scientific;

		//create a simplicial mesh
		const double r = 2.0;
		BuilderMesh mesh(2);
		mesh.createTriangleGrid(Vector2(-r, -r), Vector2(r, r), 2.0); //mesh grain is 2.0

		for (uint k = 1; k <= 12; ++k) {
			//refine the mesh into small simplices of order k
			uint order = k;
			BuilderMesh mesh_new(2);
			SmallSimplexPartition2D ssp(order);
			ssp.refineMesh(mesh, mesh_new);

			//discretise testFunctions
			Buffer<double> discreteForm1(mesh_new.getFaceSize());
			discretise2Form(mesh_new, discreteForm1, testFunction1);
			Buffer<double> discreteForm2(mesh_new.getFaceSize());
			discretise2Form(mesh_new, discreteForm2, testFunction2);
			Buffer<double> discreteForm3(mesh_new.getFaceSize());
			discretise2Form(mesh_new, discreteForm3, testFunction3);

			//compute L2 norm of the differences of testFunctions and the interpolants of discreteForms
			std::cout << "Order: " << order << '\n';
			output << "Order: " << order << '\n';
			double diff1 = computeDifferenceL2Norm2Form(mesh, ssp, discreteForm1, testFunction1);
			std::cout << "testFunction1, error L2 norm: " << diff1 << '\n';
			output << "testFunction1, error L2 norm: " << diff1 << '\n';
			double diff2 = computeDifferenceL2Norm2Form(mesh, ssp, discreteForm2, testFunction2);
			std::cout << "testFunction2, error L2 norm: " << diff2 << '\n';
			output << "testFunction2, error L2 norm: " << diff2 << '\n';
			double diff3 = computeDifferenceL2Norm2Form(mesh, ssp, discreteForm3, testFunction3);
			std::cout << "testFunction3, error L2 norm: " << diff3 << '\n';
			output << "testFunction3, error L2 norm: " << diff3 << '\n';
			output.save("Files/polynomialTestWhitney2Forms2D.txt");
		}
	}

	//test Whitney 2-forms of different orders with polynomial test functions
	void polynomialTestWhitney2Forms3D() {
		std::function<Vector3(Vector3)> testFunction1{ //testFunction1(x,y,z) = 30/128 dydz - 10/128 dzdx + 10/252 dxdy
					[](Vector3 p) -> Vector3 {
						return Vector3(30.0 / 128.0, -10.0 / 128.0, 10.0 / 252.0);
					}
		};
		std::function<Vector3(Vector3)> testFunction2{ //testFunction2(x,y,z) = x^2y^2z dydz + x^2yz^2 dzdx + xy^2z^2 dxdy
					[](Vector3 p) -> Vector3 {
						return Vector3(p.x * p.x * p.y * p.y * p.z, p.x * p.x * p.y * p.z * p.z, p.x * p.y * p.y * p.z * p.z);
					}
		};
		std::function<Vector3(Vector3)> testFunction3{ //testFunction3(x,y,z) = 20/9 * (x^2y^4z^4 dydz + x^4y^2z^4 dzdx + x^4y^4z^2 dxdy)
					[](Vector3 p) -> Vector3 {
						return (20.0 / 9.0) * Vector3(p.x * p.x * p.y * p.y * p.y * p.y * p.z * p.z * p.z * p.z, p.x * p.x * p.x * p.x * p.y * p.y * p.z * p.z * p.z * p.z,
							p.x * p.x * p.x * p.x * p.y * p.y * p.y * p.y * p.z * p.z);
					}
		};
		std::cout.precision(1);
		std::cout << std::scientific;
		Text output;
		output.precision(1);
		output << std::scientific;

		//create a simplicial mesh in Omega
		const double r = 2.0125;
		BuilderMesh mesh(3);
		mesh.createBccGrid(Vector3(-r, -r, -r), Vector3(r, r, r), 2.0); //mesh grain is 2.0
		for (uint i = mesh.getNodeSize(); i-- > 0; )
		{
			const Vector4 p = mesh.getNodePosition(i);
			if (p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.x > r - 1e-5) mesh.removeNode(i);
		}

		for (uint k = 1; k <= 12; ++k) {
			//refine the mesh into small simplices of order k
			uint order = k;
			BuilderMesh mesh_new(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh, mesh_new);

			//discretise testFunctions
			Buffer<double> discreteForm1(mesh_new.getFaceSize());
			discretise2Form(mesh_new, discreteForm1, testFunction1);
			Buffer<double> discreteForm2(mesh_new.getFaceSize());
			discretise2Form(mesh_new, discreteForm2, testFunction2);
			Buffer<double> discreteForm3(mesh_new.getFaceSize());
			discretise2Form(mesh_new, discreteForm3, testFunction3);

			//compute L2 norm of the differences of testFunctions and the interpolants of discreteForms
			std::cout << "Order: " << order << '\n';
			output << "Order: " << order << '\n';
			double diff1 = computeDifferenceL2Norm2Form(mesh, ssp, discreteForm1, testFunction1);
			std::cout << "testFunction1, error L2 norm: " << diff1 << '\n';
			output << "testFunction1, error L2 norm: " << diff1 << '\n';
			double diff2 = computeDifferenceL2Norm2Form(mesh, ssp, discreteForm2, testFunction2);
			std::cout << "testFunction2, error L2 norm: " << diff2 << '\n';
			output << "testFunction2, error L2 norm: " << diff2 << '\n';
			double diff3 = computeDifferenceL2Norm2Form(mesh, ssp, discreteForm3, testFunction3);
			std::cout << "testFunction3, error L2 norm: " << diff3 << '\n';
			output << "testFunction3, error L2 norm: " << diff3 << '\n';
			output.save("Files/polynomialTestWhitney2Forms3D.txt");
		}
	}

	//test Whitney 3-forms of different orders with polynomial test functions
	void polynomialTestWhitney3Forms3D() {
		std::function<double(Vector3)> testFunction1{ //testFunction1(x,y,z) = 0.25 dxdydz
					[](Vector3 p) -> double {
						return 0.25;
					}
		};
		std::function<double(Vector3)> testFunction2{ //testFunction2(x,y,z) = (64/75 * x^2y^2z - 8/75 * z^5) dxdydz
					[](Vector3 p) -> double {
						return (64.0 / 75.0) * p.x * p.x * p.y * p.y * p.z - (8.0 / 75.0) * p.z * p.z * p.z * p.z * p.z;
					}
		};
		std::function<double(Vector3)> testFunction3{ //testFunction3(x,y,z) = (32/11 * x^4y^4z^2 - 1/176 * z^10) dxdydz
					[](Vector3 p) -> double {
						return  (32.0 / 11.0) * p.x * p.x * p.x * p.x * p.y * p.y * p.y * p.y * p.z * p.z - p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z / 176.0;
					}
		};
		std::cout.precision(1);
		std::cout << std::scientific;
		Text output;
		output.precision(1);
		output << std::scientific;

		//create a simplicial mesh in Omega
		const double r = 2.0125;
		BuilderMesh mesh(3);
		mesh.createBccGrid(Vector3(-r, -r, -r), Vector3(r, r, r), 2.0); //mesh grain is 2.0
		for (uint i = mesh.getNodeSize(); i-- > 0; )
		{
			const Vector4 p = mesh.getNodePosition(i);
			if (p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.x - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.y > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.y > r - 1e-5) mesh.removeNode(i);
			else if (p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (p.z - p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z + p.x > r - 1e-5) mesh.removeNode(i);
			else if (-p.z - p.x > r - 1e-5) mesh.removeNode(i);
		}

		for (uint k = 1; k <= 12; ++k) {
			//refine the mesh into small simplices of order k
			uint order = k;
			BuilderMesh mesh_new(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh, mesh_new);

			//discretise testFunctions
			Buffer<double> discreteForm1(mesh_new.getBodySize());
			discretise3Form(mesh_new, discreteForm1, testFunction1);
			Buffer<double> discreteForm2(mesh_new.getBodySize());
			discretise3Form(mesh_new, discreteForm2, testFunction2);
			Buffer<double> discreteForm3(mesh_new.getBodySize());
			discretise3Form(mesh_new, discreteForm3, testFunction3);

			//compute L2 norm of the differences of testFunctions and the interpolants of discreteForms
			std::cout << "Order: " << order << '\n';
			output << "Order: " << order << '\n';
			double diff1 = computeDifferenceL2Norm3Form(mesh, ssp, discreteForm1, testFunction1);
			std::cout << "testFunction1, error L2 norm: " << diff1 << '\n';
			output << "testFunction1, error L2 norm: " << diff1 << '\n';
			double diff2 = computeDifferenceL2Norm3Form(mesh, ssp, discreteForm2, testFunction2);
			std::cout << "testFunction2, error L2 norm: " << diff2 << '\n';
			output << "testFunction2, error L2 norm: " << diff2 << '\n';
			double diff3 = computeDifferenceL2Norm3Form(mesh, ssp, discreteForm3, testFunction3);
			std::cout << "testFunction3, error L2 norm: " << diff3 << '\n';
			output << "testFunction3, error L2 norm: " << diff3 << '\n';
			output.save("Files/polynomialTestWhitney3Forms3D.txt");
		}
	}

	//test Whitney 1-forms of different orders on four meshes
	void testWhitney1Forms3D() {	
		std::function<Vector3(Vector3)> testFunction{ //testFunction(x,y,z) = 1/4 * (sin(2y)cos(2z)exp(x^2/4)dx + sin(2z)cos(2x)exp(y^2/4)dy + sin(2x)cos(2y)exp(z^2/4)dz)
			[](Vector3 p) -> Vector3 {
				return Vector3(std::sin(2 * p.y) * std::cos(2 * p.z) * std::exp(p.x * p.x / 4), std::sin(2 * p.z) * std::cos(2 * p.x) * std::exp(p.y * p.y / 4), std::sin(2 * p.x) * std::cos(2 * p.y) * std::exp(p.z * p.z / 4)) / 4.0;
			}
		};
		std::cout.precision(4);
		std::cout << std::scientific;
		Text output;
		output.precision(4);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 8; meshNumber *= 2) {
			std::cout << "mesh grain: " << 2.0 / meshNumber << '\n';
			output << "mesh grain: " << 2.0 / meshNumber << '\n';

			//create a simplicial mesh in Omega
			const double r = 2.0125;
			BuilderMesh mesh(3);
			mesh.createBccGrid(Vector3(-r, -r, -r), Vector3(r, r, r), 2.0 / meshNumber); //mesh grain is 2.0, 1.0, 0.5, and 0.25 with meshNumber 1, 2, 4, and 8 respectively
			for (uint i = mesh.getNodeSize(); i-- > 0; )
			{
				const Vector4 p = mesh.getNodePosition(i);
				if (p.x + p.y > r - 1e-5) mesh.removeNode(i);
				else if (p.x - p.y > r - 1e-5) mesh.removeNode(i);
				else if (-p.x + p.y > r - 1e-5) mesh.removeNode(i);
				else if (-p.x - p.y > r - 1e-5) mesh.removeNode(i);
				else if (p.z + p.y > r - 1e-5) mesh.removeNode(i);
				else if (p.z - p.y > r - 1e-5) mesh.removeNode(i);
				else if (-p.z + p.y > r - 1e-5) mesh.removeNode(i);
				else if (-p.z - p.y > r - 1e-5) mesh.removeNode(i);
				else if (p.z + p.x > r - 1e-5) mesh.removeNode(i);
				else if (p.z - p.x > r - 1e-5) mesh.removeNode(i);
				else if (-p.z + p.x > r - 1e-5) mesh.removeNode(i);
				else if (-p.z - p.x > r - 1e-5) mesh.removeNode(i);
			}

			for (uint k = 1; k <= 12; ++k) {
				//refine the mesh into small simplices of order k
				uint order = k;
				BuilderMesh mesh_new(3);
				SmallSimplexPartition3D ssp(order);
				ssp.refineMesh(mesh, mesh_new);

				//discretise testFunction
				Buffer<double> discreteForm(mesh_new.getEdgeSize());
				discretise1Form(mesh_new, discreteForm, testFunction);

				//compute L2 norm of the difference of testFunction and the interpolant of discreteForm
				double diff = computeDifferenceL2Norm1Form(mesh, ssp, discreteForm, testFunction);
				std::cout << "Order: " << order << '\n';
				std::cout << "Error L2 norm: " << diff << '\n';
				output << "Order: " << order << '\n';
				output << "Error L2 norm: " << diff << '\n';
				output.save("Files/testWhitney1Forms3D.txt");
			}
		}
	}

	//test cubical 1-forms of different orders on four meshes
	void testCubical1Forms2D() {
		std::function<Vector2(Vector2)> testFunction{ //testFunction(x,y,z) = 1/4 * (sin(2y)exp(x^2/4)dx + cos(2x)exp(y^2/4)dy)
			[](Vector2 p) -> Vector2 {
				return Vector2(std::sin(2 * p.y) * std::exp(p.x * p.x / 4), std::cos(2 * p.x) * std::exp(p.y * p.y / 4)) / 4.0;
			}
		};
		std::cout.precision(4);
		std::cout << std::scientific;
		Text output;
		output.precision(4);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 8; meshNumber *= 2) {
			std::cout << "mesh edge length: " << 2.0 / meshNumber << '\n';
			output << "mesh edge length: " << 2.0 / meshNumber << '\n';

			BuilderMesh mesh(2);
			createCartesianMesh(mesh, 2, 2, meshNumber, meshNumber); //mesh edge length is 2.0, 1.0, 0.5, and 0.25 with meshNumber 1, 2, 4, and 8 respectively
			for (uint i = 0; i < mesh.getNodeSize(); ++i)
				mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(1, 1, 1, 0));

			for (uint k = 1; k <= 8; ++k) {
				//refine the mesh into small cubes of order k
				uint order = k;
				BuilderMesh mesh_new(2);
				SmallCubePartition2D scp(order);
				scp.refineMesh(mesh, mesh_new);

				//discretise testFunction
				Buffer<double> discreteForm(mesh_new.getEdgeSize());
				discretise1Form(mesh_new, discreteForm, testFunction);

				//compute L2 norm of the difference of testFunction and the interpolant of discreteForm
				double diff = computeDifferenceL2Norm1Form(mesh, scp, discreteForm, testFunction, 0.25);
				std::cout << "Order: " << order << '\n';
				std::cout << "Error L2 norm: " << diff << '\n';
				output << "Order: " << order << '\n';
				output << "Error L2 norm: " << diff << '\n';
				output.save("Files/testCubical1Forms2D.txt");
			}
		}
	}

	//test cubical 1-forms of different orders on four meshes
	void testCubical1Forms3D() {
		std::function<Vector3(Vector3)> testFunction{ //testFunction(x,y,z) = 1/4 * (sin(2y)cos(2z)exp(x^2/4)dx + sin(2z)cos(2x)exp(y^2/4)dy + sin(2x)cos(2y)exp(z^2/4)dz)
			[](Vector3 p) -> Vector3 {
				return Vector3(std::sin(2 * p.y) * std::cos(2 * p.z) * std::exp(p.x * p.x / 4), std::sin(2 * p.z) * std::cos(2 * p.x) * std::exp(p.y * p.y / 4), std::sin(2 * p.x) * std::cos(2 * p.y) * std::exp(p.z * p.z / 4)) / 4.0;
			}
		};
		std::cout.precision(4);
		std::cout << std::scientific;
		Text output;
		output.precision(4);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 8; meshNumber *= 2) {
			std::cout << "mesh edge length: " << 2.0 / meshNumber << '\n';
			output << "mesh edge length: " << 2.0 / meshNumber << '\n';

			BuilderMesh mesh(3);
			createCartesianMesh(mesh, 2, 2, 2, meshNumber, meshNumber, meshNumber); //mesh edge length is 2.0, 1.0, 0.5, and 0.25 with meshNumber 1, 2, 4, and 8 respectively
			for (uint i = 0; i < mesh.getNodeSize(); ++i)
				mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(1, 1, 1, 0));

			for (uint k = 1; k <= 8; ++k) {
				//refine the mesh into small cubes of order k
				uint order = k;
				BuilderMesh mesh_new(3);
				SmallCubePartition3D scp(order);
				scp.refineMesh(mesh, mesh_new);

				//discretise testFunction
				Buffer<double> discreteForm(mesh_new.getEdgeSize());
				discretise1Form(mesh_new, discreteForm, testFunction);

				//compute L2 norm of the difference of testFunction and the interpolant of discreteForm
				double diff = computeDifferenceL2Norm1Form(mesh, scp, discreteForm, testFunction, 0.25);
				std::cout << "Order: " << order << '\n';
				std::cout << "Error L2 norm: " << diff << '\n';
				output << "Order: " << order << '\n';
				output << "Error L2 norm: " << diff << '\n';
				output.save("Files/testCubical1Forms3D.txt");
			}
		}
	}

	/*
	Study the consistency error of the discrete Hodge defined using higher order Whitney forms in the quadrilateral with vertices (-2,-2), (2,-2), (4,2), and (-4,2) with test function
	f(x,y) = (exp(x / 4) * cos(y) * (sin(x + y + 2) / 4 + cos(x + y + 2)), exp(x / 4) * (cos(y) * cos(x + y + 2) - sin(y) * sin(x + y + 2))).
	*/
	void hodgeTestWhitney1Forms2D() {
		std::function<Vector2(Vector2)> input{
						[&](Vector2 p) -> Vector2 {
							return Vector2(std::exp(p.x / 4) * std::cos(p.y) * (std::sin(p.x + p.y + 2) / 4 + std::cos(p.x + p.y + 2)), 
								std::exp(p.x / 4) * (std::cos(p.y) * std::cos(p.x + p.y + 2) - std::sin(p.y) * std::sin(p.x + p.y + 2)));
						}
		};
		std::function<Vector2(Vector2)> solution{
							[&](Vector2 p) -> Vector2 {
								return Vector2(-std::exp(p.x / 4) * (std::cos(p.y) * std::cos(p.x + p.y + 2) - std::sin(p.y) * std::sin(p.x + p.y + 2)), 
									std::exp(p.x / 4) * std::cos(p.y) * (std::sin(p.x + p.y + 2) / 4 + std::cos(p.x + p.y + 2)));
							}
		};

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 7; ++meshNumber) { //compute the consistency error on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//The following code divides the quadrilateral with vertices (-2,-2), (2,-2), (4,2), and (-4,2) into three triangles, creating a simplicial mesh.
			BuilderMesh mesh_old(2);
			const double r = 2.0;
			mesh_old.createTriangleGrid(Vector2(-r, -r), Vector2(r, r), 4.0);

			//Refine the mesh meshNumber - 1 times.
			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(2);
				SmallSimplexPartition2D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			//Compute the consistency error for orders 1-8.
			for (uint order = 1; order <= 8 && (meshNumber <= 4 || order <= 4) && (meshNumber <= 5 || order <= 2) && (meshNumber <= 6 || order <= 1); ++order) {
				BuilderMesh mesh(2);
				SmallSimplexPartition2D ssp(order);
				ssp.refineMesh(mesh_old, mesh);
				Buffer<std::unordered_map<uint, double>> star;
				ssp.formHodgeMatrix1FormsCustom(star, "DefaultTriangle", false, false);

				Buffer<double> exactSolution(mesh.getEdgeSize());
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					exactSolution[i] = integrateDual1Cell(mesh, i, solution, false);
				}

				Buffer<double> approximateSolution(mesh.getEdgeSize(), 0.0);
				Buffer<double> C_input(mesh.getEdgeSize());
				discretise1Form(mesh, C_input, input);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint j = 0; j < mesh.getEdgeSize(); ++j) {
						auto it = row.find(j);
						if (it != row.end())
							approximateSolution[i] += it->second * C_input[j];
					}
				}

				std::cout << "mesh " << meshNumber << " order " << order << '\n';
				output << "mesh " << meshNumber << " order " << order << '\n';
				double maxDiff = 0.0;
				double averageDiff = 0.0;
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					double diff = std::abs(exactSolution[i] - approximateSolution[i]) / dualEdgeLength(i, mesh, false);
					averageDiff += diff;
					if (diff > maxDiff)
						maxDiff = diff;
				}
				averageDiff /= mesh.getEdgeSize();
				std::cout << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output.save("Files/hodgeTestWhitney1Forms2D.txt");
			}
		}
	}

	/*
	Study the consistency error of the discrete Hodge defined using higher order Whitney forms in the octahedron with vertices (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,-1), and (0,0,1)
	with test function f(x,y,z) = exp(x + z) * ((cos(y * z) * (sin(x + y + 2) + cos(x + y + 2)),
								cos(y * z) * cos(x + y + 2) - z * sin(y * z) * sin(x + y + 2), 
								sin(x + y + 2) * (cos(y * z) - y * sin(y * z)))).
	*/
	void hodgeTestWhitney1Forms3D() {
		std::function<Vector3(Vector3)> input{
						[&](Vector3 p) -> Vector3 {
							return std::exp(p.x + p.z) * Vector3(std::cos(p.y * p.z) * (std::sin(p.x + p.y + 2) + std::cos(p.x + p.y + 2)),
							std::cos(p.y * p.z) * std::cos(p.x + p.y + 2) - p.z * std::sin(p.y * p.z) * std::sin(p.x + p.y + 2), std::sin(p.x + p.y + 2) * (std::cos(p.y * p.z) - p.y * std::sin(p.y * p.z)));
						}
		};

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 6; ++meshNumber) { //compute the consistency error on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//The following code divides the octahedron with vertices (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,-1), and (0,0,1) into four tetrahedra, creating a simplicial mesh.
			BuilderMesh mesh_old(3);
			const double r = 2.0125;
			createRhombicDodecahedronMesh(mesh_old, r, 2.0);
			for (uint i = mesh_old.getNodeSize(); i-- > 0; )
			{
				const Vector4 p = mesh_old.getNodePosition(i);
				if (p.z < 0 - 1e-5) mesh_old.removeNode(i);
				else if (p.x > p.z + 1e-5) mesh_old.removeNode(i);
				else if (-p.x > p.z + 1e-5) mesh_old.removeNode(i);
				else if (p.y > p.z + 1e-5) mesh_old.removeNode(i);
				else if (-p.y > p.z + 1e-5) mesh_old.removeNode(i);
				else mesh_old.setNodePosition(i, Vector4(p.x, p.y, p.z - 1.0, 0.0));
			}

			//Refine the mesh meshNumber - 1 times.
			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(3);
				SmallSimplexPartition3D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			//Compute the consistency error for orders 1-8.
			for (uint order = 1; order <= 8 && (meshNumber <= 3 || order <= 4) && (meshNumber <= 4 || order <= 2) && (meshNumber <= 5 || order <= 1); ++order) {

				BuilderMesh mesh(3);
				SmallSimplexPartition3D ssp(order);
				ssp.refineMesh(mesh_old, mesh);
				Buffer<std::unordered_map<uint, double>> star;
				ssp.formHodgeMatrix1FormsCustom(star, "BccTetrahedron", false, false);

				Buffer<double> exactSolution(mesh.getEdgeSize());
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					exactSolution[i] = integrateDual2Cell(mesh, i, input, false);
				}

				Buffer<double> approximateSolution(mesh.getEdgeSize(), 0.0);
				Buffer<double> C_input(mesh.getEdgeSize());
				discretise1Form(mesh, C_input, input);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint j = 0; j < mesh.getEdgeSize(); ++j) {
						auto it = row.find(j);
						if (it != row.end())
							approximateSolution[i] += it->second * C_input[j];
					}
				}

				std::cout << "mesh " << meshNumber << " order " << order << '\n';
				output << "mesh " << meshNumber << " order " << order << '\n';
				double maxDiff = 0.0;
				double averageDiff = 0.0;
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					double diff = std::abs(exactSolution[i] - approximateSolution[i]) / dualFaceArea(i, mesh, false);
					averageDiff += diff;
					if (diff > maxDiff)
						maxDiff = diff;
				}
				averageDiff /= mesh.getEdgeSize();
				std::cout << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output.save("Files/hodgeTestWhitney1Forms3D.txt");
			}
		}
	}

	/*
	Study the consistency error of the laplacian defined using higher order Whitney forms in the quadrilateral with vertices (-2,-2), (2,-2), (4,2), and (-4,2).
	Test function is exp(x/4)sin(x+y+2)cos(y).
	*/
	void laplacianTestWhitney2D() {
		std::function<double(Vector2)> input{
						[&](Vector2 p) -> double {
							return std::exp(p.x / 4) * std::sin(p.x + p.y + 2) * std::cos(p.y);
						}
		};
		std::function<double(Vector2)> laplacian{
							[&](Vector2 p) -> double {
								return -std::exp(p.x / 4) * (32 * std::sin(p.y) * std::cos(p.x + p.y + 2) + std::cos(p.y) * (47 * std::sin(p.x + p.y + 2) - 8 * std::cos(p.x + p.y + 2))) / 16;
							}
		};

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 7; ++meshNumber) { //compute the consistency error on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//The following code divides the quadrilateral with vertices (-2,-2), (2,-2), (4,2), and (-4,2) into three triangles, creating a simplicial mesh.
			BuilderMesh mesh_old(2);
			const double r = 2.0;
			mesh_old.createTriangleGrid(Vector2(-r, -r), Vector2(r, r), 4.0);

			//Refine the mesh meshNumber - 1 times.
			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(2);
				SmallSimplexPartition2D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			//Compute the consistency error for orders 1-8.
			for (uint order = 1; order <= 8 && (meshNumber <= 4 || order <= 4) && (meshNumber <= 5 || order <= 2) && (meshNumber <= 6 || order <= 1); ++order) {
				BuilderMesh mesh(2);
				SmallSimplexPartition2D ssp(order);
				ssp.refineMesh(mesh_old, mesh);
				Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
				{
					Buffer<std::unordered_map<uint, double>> star;
					ssp.formHodgeMatrix1FormsCustom(star, "DefaultTriangle", false, false);
					for (uint j = 0; j < mesh.getNodeSize(); ++j) {
						Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
						const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
						for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
							std::unordered_map<uint, double>& row = star[i];
							for (uint e = 0; e < edges_j.size(); ++e) {
								auto it = row.find(edges_j[e]);
								if (it != row.end())
									col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
							}
						}
						for (uint i = 0; i < mesh.getNodeSize(); ++i) {
							const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
							for (uint e = 0; e < edges_i.size(); ++e) {
								if (col_j_of_star_d[edges_i[e]] != 0.0)
									system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}

				Buffer<double> exactSolution(mesh.getNodeSize());
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					//exactSolution[i] = integrateDual2Cell(mesh, i, laplacian);
					exactSolution[i] = integrateDual2CellInParts(mesh, i, laplacian, 0.25);
				}

				Buffer<double> approximateSolution(mesh.getNodeSize(), 0.0);
				Buffer<double> C_input(mesh.getNodeSize());
				discretise0Form(mesh, C_input, input);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						approximateSolution[i] += it->second * C_input[it->first];
					}
				}

				std::cout << "mesh " << meshNumber << " order " << order << '\n';
				output << "mesh " << meshNumber << " order " << order << '\n';
				double maxDiff = 0.0;
				double averageDiff = 0.0;
				uint interiorNodes = 0;
				mesh.fillBoundaryFlags(1);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						continue;
					double diff = std::abs(exactSolution[i] - approximateSolution[i]) / dualFaceArea(i, mesh, false); //choose dualFaceBoundaryLength or dualFaceArea
					averageDiff += diff;
					if (diff > maxDiff)
						maxDiff = diff;
					++interiorNodes;
				}
				averageDiff /= interiorNodes;
				std::cout << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output.save("Files/laplacianTestWhitney2D.txt");
			}
		}
	}

	/*
	Study the consistency error of the laplacian defined using higher order Whitney forms in the octahedron with vertices (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,-1), and (0,0,1).
	Test function is exp(x+z)sin(x+y+2)cos(y*z).
	*/
	void laplacianTestWhitney3D() {
		std::function<double(Vector3)> input{
						[&](Vector3 p) -> double {
							return std::exp(p.x + p.z) * std::sin(p.x + p.y + 2) * std::cos(p.y * p.z);
						}
		};
		std::function<double(Vector3)> laplacian{
							[&](Vector3 p) -> double {
								return std::exp(p.x + p.z) * (std::cos(p.y * p.z) * (-(p.z * p.z + p.y * p.y) * std::sin(p.x + p.y + 2) + 2 * std::cos(p.x + p.y + 2))
								- 2 * std::sin(p.y * p.z) * (p.y * std::sin(p.x + p.y + 2) + p.z * std::cos(p.x + p.y + 2)));
							}
		};

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 6; ++meshNumber) { //compute the consistency error on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//The following code divides the octahedron with vertices (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,-1), and (0,0,1) into four tetrahedra, creating a simplicial mesh.
			BuilderMesh mesh_old(3);
			const double r = 2.0125;
			createRhombicDodecahedronMesh(mesh_old, r, 2.0);
			for (uint i = mesh_old.getNodeSize(); i-- > 0; )
			{
				const Vector4 p = mesh_old.getNodePosition(i);
				if (p.z < 0 - 1e-5) mesh_old.removeNode(i);
				else if (p.x > p.z + 1e-5) mesh_old.removeNode(i);
				else if (-p.x > p.z + 1e-5) mesh_old.removeNode(i);
				else if (p.y > p.z + 1e-5) mesh_old.removeNode(i);
				else if (-p.y > p.z + 1e-5) mesh_old.removeNode(i);
				else mesh_old.setNodePosition(i, Vector4(p.x, p.y, p.z - 1.0, 0.0));
			}

			//Refine the mesh meshNumber - 1 times.
			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(3);
				SmallSimplexPartition3D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			//Compute the consistency error for orders 1-8.
			for (uint order = 1; order <= 8 && (meshNumber <= 3 || order <= 4) && (meshNumber <= 4 || order <= 2) && (meshNumber <= 5 || order <= 1); ++order) {

				BuilderMesh mesh(3);
				SmallSimplexPartition3D ssp(order);
				ssp.refineMesh(mesh_old, mesh);
				Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
				{
					Buffer<std::unordered_map<uint, double>> star;
					ssp.formHodgeMatrix1FormsCustom(star, "BccTetrahedron", false, false);
					for (uint j = 0; j < mesh.getNodeSize(); ++j) {
						Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
						const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
						for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
							std::unordered_map<uint, double>& row = star[i];
							for (uint e = 0; e < edges_j.size(); ++e) {
								auto it = row.find(edges_j[e]);
								if (it != row.end())
									col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
							}
						}
						for (uint i = 0; i < mesh.getNodeSize(); ++i) {
							const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
							for (uint e = 0; e < edges_i.size(); ++e) {
								if (col_j_of_star_d[edges_i[e]] != 0.0)
									system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}

				Buffer<double> exactSolution(mesh.getNodeSize());
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					//exactSolution[i] = integrateDual3Cell(mesh, i, laplacian);
					exactSolution[i] = integrateDual3CellInParts(mesh, i, laplacian, 0.25);
				}

				Buffer<double> approximateSolution(mesh.getNodeSize(), 0.0);
				Buffer<double> C_input(mesh.getNodeSize());
				discretise0Form(mesh, C_input, input);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						approximateSolution[i] += it->second * C_input[it->first];
					}
				}

				std::cout << "mesh " << meshNumber << " order " << order << '\n';
				output << "mesh " << meshNumber << " order " << order << '\n';
				double maxDiff = 0.0;
				double averageDiff = 0.0;
				uint interiorNodes = 0;
				mesh.fillBoundaryFlags(1);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						continue;
					double diff = std::abs(exactSolution[i] - approximateSolution[i]) / dualCellVolume(i, mesh, false); //choose dualCellBoundaryArea or dualCellVolume
					averageDiff += diff;
					if (diff > maxDiff)
						maxDiff = diff;
					++interiorNodes;
				}
				averageDiff /= interiorNodes;
				std::cout << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output.save("Files/laplacianTestWhitney3D.txt");
			}
		}
	}

	/*
	Study the consistency error of the operator d*d defined using higher order Whitney forms in Minkowski space [0, L] x [0, T] with L=2 and T=2.
	Test function is exp(x/4)sin(x+t+2)cos(t).
	*/
	void minkowskiLaplacianTestWhitney2D() {
		const double L = 2.0;
		const double T = 2.0;
		std::function<double(Vector2)> input{
						[&](Vector2 p) -> double {
							return std::exp(p.x / 4) * std::sin(p.x + p.y + 2) * std::cos(p.y);
						}
		};
		std::function<double(Vector2)> laplacian{
							[&](Vector2 p) -> double {
								return std::exp(p.x / 4) * (0.5 * std::cos(p.y) * std::cos(2 + p.x + p.y) + 2 * std::cos(2 + p.x + p.y) * std::sin(p.y) + 17.0 / 16 * std::cos(p.y) * std::sin(2 + p.x + p.y));
							}
		};

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 7; ++meshNumber) { //compute the consistency error on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//Create triangle mesh in [0, 2] x [0, 2].
			BuilderMesh mesh_old(2);
			createSpacetimeMesh(mesh_old, L, T, 2, 2);

			//Refine the mesh meshNumber - 1 times.
			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(2);
				SmallSimplexPartition2D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			//Compute the consistency error for orders 1-6.
			for (uint order = 1; order <= 6 && (meshNumber <= 4 || order <= 5) && (meshNumber <= 5 || order <= 2) && (meshNumber <= 6 || order <= 1); ++order) {
				BuilderMesh mesh(2);
				SmallSimplexPartition2D ssp(order);
				ssp.refineMesh(mesh_old, mesh);
				Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
				{
					Buffer<std::unordered_map<uint, double>> star;
					Buffer<uint> refElementIndices(1, 0);
					Buffer <Buffer<uint>> refElements(1);
					refElements[0].resize(mesh_old.getFaceSize());
					for (uint i = 0; i < mesh_old.getFaceSize(); ++i)
						refElements[0][i] = i;
					ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "RightTriangle", false, true, false);
					for (uint j = 0; j < mesh.getNodeSize(); ++j) {
						Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
						const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
						for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
							std::unordered_map<uint, double>& row = star[i];
							for (uint e = 0; e < edges_j.size(); ++e) {
								auto it = row.find(edges_j[e]);
								if (it != row.end())
									col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
							}
						}
						for (uint i = 0; i < mesh.getNodeSize(); ++i) {
							const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
							for (uint e = 0; e < edges_i.size(); ++e) {
								if (col_j_of_star_d[edges_i[e]] != 0.0)
									system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}

				Buffer<double> exactSolution(mesh.getNodeSize());
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					exactSolution[i] = integrateDual2Cell(mesh, i, laplacian);
					//exactSolution[i] = integrateDual2CellInParts(mesh, i, laplacian, 0.25);
				}

				Buffer<double> approximateSolution(mesh.getNodeSize(), 0.0);
				Buffer<double> C_input(mesh.getNodeSize());
				discretise0Form(mesh, C_input, input);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						approximateSolution[i] += it->second * C_input[it->first];
					}
				}

				std::cout << "mesh " << meshNumber << " order " << order << '\n';
				output << "mesh " << meshNumber << " order " << order << '\n';
				double maxDiff = 0.0;
				double averageDiff = 0.0;
				uint interiorNodes = 0;
				mesh.fillBoundaryFlags(1);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						continue;
					double diff = std::abs(exactSolution[i] - approximateSolution[i]) / dualFaceArea(i, mesh, false); //choose dualFaceBoundaryLength or dualFaceArea
					averageDiff += diff;
					if (diff > maxDiff)
						maxDiff = diff;
					++interiorNodes;
				}
				averageDiff /= interiorNodes;
				std::cout << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output.save("Files/minkowskiLaplacianTestWhitney2D.txt");
			}
		}
	}

	/*
	Study the consistency error of the operator d*d defined using higher order Whitney forms in Minkowski space [0, Lx] x [0, Ly] x [0, T] with Lx=2, Ly=2, and T=2.
	Test function is exp(x+t)sin(x+y+2)cos(y*t).
	*/
	void minkowskiLaplacianTestWhitney3D() {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double T = 2.0;
		std::function<double(Vector3)> input{
						[&](Vector3 p) -> double {
							return std::exp(p.x + p.z) * std::sin(p.x + p.y + 2) * std::cos(p.y * p.z);
						}
		};
		std::function<double(Vector3)> laplacian{
							[&](Vector3 p) -> double {
								return std::exp(p.x + p.z) * (std::sin(p.x + p.y + 2) * ((p.y * p.y - p.z * p.z - 2) * std::cos(p.y * p.z) + 2 * p.y * std::sin(p.y * p.z))
									+ 2 * std::cos(p.x + p.y + 2) * (std::cos(p.y * p.z) - p.z * std::sin(p.y * p.z)));
							}
		};

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 6; ++meshNumber) { //compute the consistency error on meshes of different grain
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//Create tetrahedral mesh in [0, 2] x [0, 2] x [0, 2].
			BuilderMesh mesh_old(3);
			uint refs = std::pow(2, meshNumber - 1);
			createSpacetimeMesh(mesh_old, Lx, Ly, T, refs, refs, refs);

			//Compute the consistency error for orders 1-5.
			for (uint order = 1; order <= 5 && (meshNumber <= 3 || order <= 4) && (meshNumber <= 4 || order <= 2) && (meshNumber <= 5 || order <= 1); ++order) {

				BuilderMesh mesh(3);
				SmallSimplexPartition3D ssp(order);
				ssp.refineMesh(mesh_old, mesh);
				Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
				{
					Buffer<std::unordered_map<uint, double>> star;
					Buffer<uint> refElementIndices(6);
					Buffer<Buffer<uint>> refElements(6);
					for (uint i = 0; i < 6; ++i) {
						refElementIndices[i] = i;
						refElements[i].resize(mesh_old.getBodySize() / 6);
					}

					for (uint i = 0; i < mesh_old.getBodySize(); ++i)
						refElements[i % 6][i / 6] = i;
					ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "CubeTetrahedron", false, true, false);
					for (uint j = 0; j < mesh.getNodeSize(); ++j) {
						Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
						const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
						for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
							std::unordered_map<uint, double>& row = star[i];
							for (uint e = 0; e < edges_j.size(); ++e) {
								auto it = row.find(edges_j[e]);
								if (it != row.end())
									col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
							}
						}
						for (uint i = 0; i < mesh.getNodeSize(); ++i) {
							const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
							for (uint e = 0; e < edges_i.size(); ++e) {
								if (col_j_of_star_d[edges_i[e]] != 0.0)
									system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}

				Buffer<double> exactSolution(mesh.getNodeSize());
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					exactSolution[i] = integrateDual3Cell(mesh, i, laplacian);
					//exactSolution[i] = integrateDual3CellInParts(mesh, i, laplacian, 0.25);
				}

				Buffer<double> approximateSolution(mesh.getNodeSize(), 0.0);
				Buffer<double> C_input(mesh.getNodeSize());
				discretise0Form(mesh, C_input, input);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						approximateSolution[i] += it->second * C_input[it->first];
					}
				}

				std::cout << "mesh " << meshNumber << " order " << order << '\n';
				output << "mesh " << meshNumber << " order " << order << '\n';
				double maxDiff = 0.0;
				double averageDiff = 0.0;
				uint interiorNodes = 0;
				mesh.fillBoundaryFlags(1);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						continue;
					double diff = std::abs(exactSolution[i] - approximateSolution[i]) / dualCellVolume(i, mesh, false); //choose dualCellBoundaryArea or dualCellVolume
					averageDiff += diff;
					if (diff > maxDiff)
						maxDiff = diff;
					++interiorNodes;
				}
				averageDiff /= interiorNodes;
				std::cout << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output << "Max consistency error " << maxDiff << " (average " << averageDiff << ')' << '\n' << '\n';
				output.save("Files/minkowskiLaplacianTestWhitney3D.txt");
			}
		}
	}

	/*
	Solves Poisson's equation in the quadrilateral with vertices (-2,-2), (2,-2), (4,2), and (-4,2).
	Dirichlet boundary condition is imposed on the line segment from (-2,-2) to (2,-2).
	Neumann boundary condition is imposed elsewhere on the boundary.
	The exact solution is the test function a), b), or c) chosen with parameter testFunc.
	This implies boundary conditions and the source term.
	The second parameter tells whether to use circumcentric (if true) or barycentric duals.
	*/
	void polynomialPoissonExampleWhitney2D(const char testFunc, bool circumcentric) {
		std::function<double(Vector2)> testFunction;
		std::function<double(Vector2)> sourceTerm;
		std::function<Vector2(Vector2)> testFunctionGradient;
		if (testFunc == 'a') { //testFunction a) - 1st order Whitney form
			testFunction = [&](Vector2 p) -> double {
				return 2.5 * p.x - p.y;
			};
			sourceTerm = [&](Vector2 p) -> double {
				return 0.0;
			};
			testFunctionGradient = [&](Vector2 p) -> Vector2 {
				return Vector2(2.5, -1.0);
			};
		}
		else if (testFunc == 'b') { //testFunction b) - 5th order Whitney form
			testFunction = [&](Vector2 p) -> double {
				return p.x * p.x * p.x * p.y * p.y / 8.0 - p.y * p.y * p.y * p.y / 2.0;
			};
			sourceTerm = [&](Vector2 p) -> double {
				return (3.0 * p.x * p.y * p.y + p.x * p.x * p.x) / 4.0 - 6.0 * p.y * p.y;
			};
			testFunctionGradient = [&](Vector2 p) -> Vector2 {
				return Vector2((3.0 / 8.0) * p.x * p.x * p.y * p.y, p.x * p.x * p.x * p.y / 4.0 - 2.0 * p.y * p.y * p.y);
			};
		}
		else if (testFunc == 'c') { //testFunction c) - 10th order Whitney form
			testFunction = [&](Vector2 p) -> double {
				return p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y * p.y * p.y / 2048.0 + p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y / 64.0;
			};
			sourceTerm = [&](Vector2 p) -> double {
				return (42.0 / 2048.0) * p.x * p.x * p.x * p.x * p.x * p.y * p.y * p.y + (6.0 / 2048.0) * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y
					+ (90.0 / 64.0) * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y;
			};
			testFunctionGradient = [&](Vector2 p) -> Vector2 {
				return Vector2((7.0 / 2048.0) * p.x * p.x * p.x * p.x * p.x * p.x * p.y * p.y * p.y,
					(3.0 / 2048.0) * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y * p.y + (10.0 / 64.0) * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y);
			};
		}
		else {
			std::cout << "polynomialPoissonExampleWhitney2D a), b), or c)" << '\n';
			return;
		}

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text filepath;
		filepath << "Files/polynomialPoissonExampleWhitney2D_" << testFunc << (circumcentric ? "_circumcentric" : "_barycentric") << ".txt";

		//The following code divides the quadrilateral with vertices (-2,-2), (2,-2), (4,2), and (-4,2) into three triangles, creating a simplicial mesh.
		BuilderMesh mesh_old(2);
		const double r = 2.0;
		mesh_old.createTriangleGrid(Vector2(-r, -r), Vector2(r, r), 4.0);

		//mark boundary nodes and edges to indicate whether Dirichlet (flag = 1) or Neumann (flag = 2) boundary condition is imposed
		mesh_old.fillBoundaryFlags(1);
		for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
			if (mesh_old.getNodeFlag(i) == 1 && mesh_old.getNodePosition2(i).y > -2.0 + 1.0e-14)
				mesh_old.setNodeFlag(i, 2);
		}
		for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
			if (mesh_old.getEdgeFlag(i) == 1 && mesh_old.getEdgeAverage2(i).y > -2.0 + 1.0e-14)
				mesh_old.setEdgeFlag(i, 2);
		}

		//solve the problem using Whitney forms of orders 1-10
		for (uint order = 1; order <= 10; ++order) {

			//refine the mesh into small simplices
			BuilderMesh mesh(2);
			SmallSimplexPartition2D ssp(order);
			ssp.refineMesh(mesh_old, mesh);

			//similarly mark boundary nodes and edges of refined mesh
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 && mesh.getNodePosition2(i).y > -2.0 + 1.0e-14)
					mesh.setNodeFlag(i, 2);
			}
			for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
				if (mesh.getEdgeFlag(i) == 1 && mesh.getEdgeAverage2(i).y > -2.0 + 1.0e-14)
					mesh.setEdgeFlag(i, 2);
			}

			//compute discrete solution using Hodge operator of given order
			Buffer<double> dec_sol;
			getDECSolutionPoissonEq(mesh, ssp, sourceTerm, testFunction, testFunctionGradient, dec_sol, circumcentric, true);

			//error analysis
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//maximum and average error
			double maxError = 0.0;
			double averageError = 0.0;
			uint interiorNodes = 0;
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1)
					continue;
				++interiorNodes;
				double diff = std::abs(dec_sol[i] - testFunction(mesh.getNodePosition2(i)));
				if (maxError < diff)
					maxError = diff;
				averageError += diff;
			}
			averageError /= interiorNodes;
			std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

			//L2 error
			double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, dec_sol, testFunction, 1000);
			std::cout << "L2 error: " << l2error1 << '\n';
			output << "L2 error: " << l2error1 << '\n';

			//H1 error
			Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
			for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
				const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
				for (uint j = 0; j < nodes.size(); ++j)
					dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * dec_sol[nodes[j]];
			}
			double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
			std::cout << "Gradient L2 error: " << graderror1 << '\n';
			output << "Gradient L2 error: " << graderror1 << '\n';
			double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
			std::cout << "H1 error: " << h1error1 << '\n' << '\n';
			output << "H1 error: " << h1error1 << '\n' << '\n';
			output.save(filepath.str());
		}
	}

	/*
	Solves Poisson's equation with inhomogenous Dirichlet boundary condition in the quadrilateral with vertices (-2,-2), (2,-2), (4,2), and (-4,2).
	The exact solution is chosen to be exp(x/4)sin(x+y+2)cos(y).
	This implies the boundary condition and the source term.
	The parameter tells whether to use circumcentric (if true) or barycentric duals.
	*/
	void poissonExampleWhitney2D(bool circumcentric) {
		std::function<double(Vector2)> testFunction{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * std::sin(p.x + p.y + 2) * std::cos(p.y); } };
		std::function<double(Vector2)> sourceTerm{ [&](Vector2 p) -> double { return -std::exp(p.x / 4) * (32 * std::sin(p.y) * std::cos(p.x + p.y + 2)
			+ std::cos(p.y) * (47 * std::sin(p.x + p.y + 2) - 8 * std::cos(p.x + p.y + 2))) / 16; } };
		std::function<Vector2(Vector2)> testFunctionGradient{ [&](Vector2 p) -> Vector2 { return Vector2(std::exp(p.x / 4) * std::cos(p.y) * (std::sin(p.x + p.y + 2) / 4 + std::cos(p.x + p.y + 2)),
			std::exp(p.x / 4) * (std::cos(p.y) * std::cos(p.x + p.y + 2) - std::sin(p.y) * std::sin(p.x + p.y + 2))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text filepath;
		filepath << "Files/poissonExampleWhitney2D" << (circumcentric ? "_circumcentric" : "_barycentric") << ".txt";

		for (uint meshNumber = 1; meshNumber <= 7; ++meshNumber) { //solve the problem on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//The following code divides the quadrilateral with vertices (-2,-2), (2,-2), (4,2), and (-4,2) into three triangles, creating a simplicial mesh.
			BuilderMesh mesh_old(2);
			const double r = 2.0;
			mesh_old.createTriangleGrid(Vector2(-r, -r), Vector2(r, r), 4.0);

			//Refine the mesh meshNumber - 1 times.
			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(2);
				SmallSimplexPartition2D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			//mark boundary nodes and edges to indicate Dirichlet boundary condition (flag = 1)
			mesh_old.fillBoundaryFlags(1);

			//solve the problem using Whitney forms of orders 1-8
			for (uint order = 1; order <= 8 && (meshNumber <= 4 || order <= 4) && (meshNumber <= 5 || order <= 2) && (meshNumber <= 6 || order <= 1); ++order) {

				//refine the mesh into small simplices
				BuilderMesh mesh(2);
				SmallSimplexPartition2D ssp(order);
				ssp.refineMesh(mesh_old, mesh);

				//similarly mark boundary nodes and edges of refined mesh
				mesh.fillBoundaryFlags(1);

				//compute discrete solution using Hodge operator of given order
				Buffer<double> dec_sol;
				getDECSolutionPoissonEq(mesh, ssp, sourceTerm, testFunction, testFunctionGradient, dec_sol, circumcentric, true);

				//error analysis
				std::cout << "DEC solution of order " << order << '\n';
				output << "DEC solution of order " << order << '\n';

				//maximum and average error
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						continue;
					++interiorNodes;
					double diff = std::abs(dec_sol[i] - testFunction(mesh.getNodePosition2(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageError /= interiorNodes;
				std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
				output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

				//L2 error
				double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, dec_sol, testFunction, 1000);
				std::cout << "L2 error: " << l2error1 << '\n';
				output << "L2 error: " << l2error1 << '\n';

				//H1 error
				Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
					for (uint j = 0; j < nodes.size(); ++j)
						dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * dec_sol[nodes[j]];
				}
				double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
				std::cout << "Gradient L2 error: " << graderror1 << '\n';
				output << "Gradient L2 error: " << graderror1 << '\n';
				double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
				std::cout << "H1 error: " << h1error1 << '\n' << '\n';
				output << "H1 error: " << h1error1 << '\n' << '\n';

			}
			std::cout << '\n';
			output << '\n';
			output.save(filepath.str());
		}
	}

	/*
	Solves Poisson's equation in the octahedron with vertices (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,-1), and (0,0,1).
	Dirichlet boundary condition is imposed on the part of the boundary where z >= 0.
	Neumann boundary condition is imposed elsewhere on the boundary.
	The exact solution is the test function a), b), or c) chosen with parameter testFunc.
	This implies boundary conditions and the source term.
	The second parameter tells whether to use circumcentric (if true) or barycentric duals.
	*/
	void polynomialPoissonExampleWhitney3D(const char testFunc, bool circumcentric) {
		std::function<double(Vector3)> testFunction;
		std::function<double(Vector3)> sourceTerm;
		std::function<Vector3(Vector3)> testFunctionGradient;
		if (testFunc == 'a') { //testFunction a) - 1st order Whitney form
			testFunction = [&](Vector3 p) -> double {
				return 2.5 * p.x - p.y + 1.25 * p.z;
			};
			sourceTerm = [&](Vector3 p) -> double {
				return 0.0;
			};
			testFunctionGradient = [&](Vector3 p) -> Vector3 {
				return Vector3(2.5, -1.0, 1.25);
			};
		}
		else if (testFunc == 'b') { //testFunction b) - 5th order Whitney form
			testFunction = [&](Vector3 p) -> double {
				return 10 * p.x * p.x * p.y * p.y * p.y - 40 * p.z * p.z * p.z * p.z * p.z;
			};
			sourceTerm = [&](Vector3 p) -> double {
				return 20 * p.y * p.y * p.y + 60 * p.x * p.x * p.y - 800 * p.z * p.z * p.z;
			};
			testFunctionGradient = [&](Vector3 p) -> Vector3 {
				return Vector3(20 * p.x * p.y * p.y * p.y, 30 * p.x * p.x * p.y * p.y, -200 * p.z * p.z * p.z * p.z);
			};
		}
		else if (testFunc == 'c') { //testFunction c) - 10th order Whitney form
			testFunction = [&](Vector3 p) -> double {
				return 20 * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x - 1000 * p.x * p.y * p.y * p.y * p.y * p.y * p.y * p.z * p.z * p.z;
			};
			sourceTerm = [&](Vector3 p) -> double {
				return 1800 * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x - 30000 * p.x * p.y * p.y * p.y * p.y * p.z * p.z * p.z - 6000 * p.x * p.y * p.y * p.y * p.y * p.y * p.y * p.z;
			};
			testFunctionGradient = [&](Vector3 p) -> Vector3 {
				return Vector3(200 * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x - 1000 * p.y * p.y * p.y * p.y * p.y * p.y * p.z * p.z * p.z,
					-6000 * p.x * p.y * p.y * p.y * p.y * p.y * p.z * p.z * p.z, -3000 * p.x * p.y * p.y * p.y * p.y * p.y * p.y * p.z * p.z);
			};
		}
		else {
			std::cout << "Choose polynomialPoissonExampleWhitney3D a), b), or c)" << '\n';
			return;
		}

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text filepath;
		filepath << "Files/polynomialPoissonExampleWhitney3D_" << testFunc << (circumcentric ? "_circumcentric" : "_barycentric") << ".txt";

		//The following code divides the octahedron with vertices (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,-1), and (0,0,1) into four tetrahedra, creating a simplicial mesh.
		BuilderMesh mesh_old(3);
		const double r = 2.0125;
		createRhombicDodecahedronMesh(mesh_old, r, 2.0);
		for (uint i = mesh_old.getNodeSize(); i-- > 0; )
		{
			const Vector4 p = mesh_old.getNodePosition(i);
			if (p.z < 0 - 1e-5) mesh_old.removeNode(i);
			else if (p.x > p.z + 1e-5) mesh_old.removeNode(i);
			else if (-p.x > p.z + 1e-5) mesh_old.removeNode(i);
			else if (p.y > p.z + 1e-5) mesh_old.removeNode(i);
			else if (-p.y > p.z + 1e-5) mesh_old.removeNode(i);
			else mesh_old.setNodePosition(i, Vector4(p.x, p.y, p.z - 1.0, 0.0));
		}

		//mark boundary nodes, edges, and faces to indicate whether Dirichlet (flag = 1) or Neumann (flag = 2) boundary condition is imposed
		mesh_old.fillBoundaryFlags(1);
		for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
			if (mesh_old.getNodeFlag(i) == 1 && mesh_old.getNodePosition3(i).z < -1.0e-14)
				mesh_old.setNodeFlag(i, 2);
		}
		for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
			if (mesh_old.getEdgeFlag(i) == 1 && mesh_old.getEdgeAverage3(i).z < -1.0e-14)
				mesh_old.setEdgeFlag(i, 2);
		}

		//solve the problem using Whitney forms of orders 1-10
		for (uint order = 1; order <= 10; ++order) {

			//refine the mesh into small simplices
			BuilderMesh mesh(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh_old, mesh);

			//similarly mark boundary nodes and edges of refined mesh
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 && mesh.getNodePosition3(i).z < -1.0e-14)
					mesh.setNodeFlag(i, 2);
			}
			for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
				if (mesh.getEdgeFlag(i) == 1 && mesh.getEdgeAverage3(i).z < -1.0e-14)
					mesh.setEdgeFlag(i, 2);
			}

			//compute discrete solution using Hodge operator of given order
			Buffer<double> dec_sol;
			getDECSolutionPoissonEq(mesh, ssp, sourceTerm, testFunction, testFunctionGradient, dec_sol, circumcentric, true);

			//error analysis
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//maximum and average error
			double maxError = 0.0;
			double averageError = 0.0;
			uint interiorNodes = 0;
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1)
					continue;
				++interiorNodes;
				double diff = std::abs(dec_sol[i] - testFunction(mesh.getNodePosition3(i)));
				if (maxError < diff)
					maxError = diff;
				averageError += diff;
			}
			averageError /= interiorNodes;
			std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

			//L2 error
			double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, dec_sol, testFunction, 1000);
			std::cout << "L2 error: " << l2error1 << '\n';
			output << "L2 error: " << l2error1 << '\n';

			//H1 error
			Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
			for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
				const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
				for (uint j = 0; j < nodes.size(); ++j)
					dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * dec_sol[nodes[j]];
			}
			double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
			std::cout << "Gradient L2 error: " << graderror1 << '\n';
			output << "Gradient L2 error: " << graderror1 << '\n';
			double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
			std::cout << "H1 error: " << h1error1 << '\n' << '\n';
			output << "H1 error: " << h1error1 << '\n' << '\n';
			output.save(filepath.str());
		}
	}

	/*
	Solves Poisson's equation with inhomogenous Dirichlet boundary condition in the octahedron with vertices (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,-1), and (0,0,1).
	The exact solution is chosen to be exp(x+z)sin(x+y+2)cos(y*z).
	This implies the boundary condition and the source term.
	The parameter tells whether to use circumcentric (if true) or barycentric duals.
	*/
	void poissonExampleWhitney3D(bool circumcentric) {
		std::function<double(Vector3)> testFunction{ [&](Vector3 p) -> double { return std::exp(p.x + p.z) * std::sin(p.x + p.y + 2) * std::cos(p.y * p.z); } };
		std::function<double(Vector3)> sourceTerm{ [&](Vector3 p) -> double { return std::exp(p.x + p.z) * (std::cos(p.y * p.z) * (-(p.z * p.z + p.y * p.y) * std::sin(p.x + p.y + 2) + 2 * std::cos(p.x + p.y + 2))
			- 2 * std::sin(p.y * p.z) * (p.y * std::sin(p.x + p.y + 2) + p.z * std::cos(p.x + p.y + 2))); } };
		std::function<Vector3(Vector3)> testFunctionGradient{ [&](Vector3 p) -> Vector3 { return std::exp(p.x + p.z) * Vector3(std::cos(p.y * p.z) * (std::sin(p.x + p.y + 2) + std::cos(p.x + p.y + 2)),
			std::cos(p.y * p.z) * std::cos(p.x + p.y + 2) - p.z * std::sin(p.y * p.z) * std::sin(p.x + p.y + 2), std::sin(p.x + p.y + 2) * (std::cos(p.y * p.z) - p.y * std::sin(p.y * p.z))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text filepath;
		filepath << "Files/poissonExampleWhitney3D" << (circumcentric ? "_circumcentric" : "_barycentric") << ".txt";

		for (uint meshNumber = 1; meshNumber <= 6; ++meshNumber) { //solve the problem on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//The following code divides the octahedron with vertices (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (0,0,-1), and (0,0,1) into four tetrahedra, creating a simplicial mesh.
			BuilderMesh mesh_old(3);
			const double r = 2.0125;
			createRhombicDodecahedronMesh(mesh_old, r, 2.0);
			for (uint i = mesh_old.getNodeSize(); i-- > 0; )
			{
				const Vector4 p = mesh_old.getNodePosition(i);
				if (p.z < 0 - 1e-5) mesh_old.removeNode(i);
				else if (p.x > p.z + 1e-5) mesh_old.removeNode(i);
				else if (-p.x > p.z + 1e-5) mesh_old.removeNode(i);
				else if (p.y > p.z + 1e-5) mesh_old.removeNode(i);
				else if (-p.y > p.z + 1e-5) mesh_old.removeNode(i);
				else mesh_old.setNodePosition(i, Vector4(p.x, p.y, p.z - 1.0, 0.0));
			}

			//Refine the mesh meshNumber - 1 times.
			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(3);
				SmallSimplexPartition3D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			//mark boundary nodes and edges to indicate Dirichlet boundary condition (flag = 1)
			mesh_old.fillBoundaryFlags(1);

			//solve the problem using Whitney forms of orders 1-8
			for (uint order = 1; order <= 8 && (meshNumber <= 3 || order <= 4) && (meshNumber <= 4 || order <= 2) && (meshNumber <= 5 || order <= 1); ++order) {

				//refine the mesh into small simplices
				BuilderMesh mesh(3);
				SmallSimplexPartition3D ssp(order);
				ssp.refineMesh(mesh_old, mesh);

				//similarly mark boundary nodes and edges of refined mesh
				mesh.fillBoundaryFlags(1);

				//compute discrete solution using Hodge operator of given order
				Buffer<double> dec_sol;
				getDECSolutionPoissonEq(mesh, ssp, sourceTerm, testFunction, testFunctionGradient, dec_sol, circumcentric, true);

				//error analysis
				std::cout << "DEC solution of order " << order << '\n';
				output << "DEC solution of order " << order << '\n';

				//maximum and average error
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						continue;
					++interiorNodes;
					double diff = std::abs(dec_sol[i] - testFunction(mesh.getNodePosition3(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageError /= interiorNodes;
				std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
				output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

				//L2 error
				double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, dec_sol, testFunction, 1000);
				std::cout << "L2 error: " << l2error1 << '\n';
				output << "L2 error: " << l2error1 << '\n';

				//H1 error
				Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
					for (uint j = 0; j < nodes.size(); ++j)
						dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * dec_sol[nodes[j]];
				}
				double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
				std::cout << "Gradient L2 error: " << graderror1 << '\n';
				output << "Gradient L2 error: " << graderror1 << '\n';
				double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
				std::cout << "H1 error: " << h1error1 << '\n' << '\n';
				output << "H1 error: " << h1error1 << '\n' << '\n';
			}
			std::cout << '\n';
			output << '\n';
			output.save(filepath.str());
		}
	}

	/*
	Solves inhomogenous wave equation in domain [0, L] x [0, T] with L=2, T=2, using higher order Whitney forms and discrete Hodge operators in spacetime.
	The exact solution is the test function a), b), or c) chosen with parameter testFunc.
	This implies boundary conditions and the source term.
	*/
	void polynomialWaveExampleWhitney2D(const char testFunc) {
		const double L = 2.0;
		const double T = 2.0;
		std::function<double(Vector2)> testFunction;
		std::function<double(Vector2)> sourceTerm;
		std::function<Vector2(Vector2)> testFunctionGradient;
		if (testFunc == 'a') { //testFunction a) - 1st order Whitney form
			testFunction = [&](Vector2 p) -> double {
				return 2.5 * p.x - p.y;
			};
			sourceTerm = [&](Vector2 p) -> double {
				return 0.0;
			};
			testFunctionGradient = [&](Vector2 p) -> Vector2 {
				return Vector2(2.5, -1.0);
			};
		}
		else if (testFunc == 'b') { //testFunction b) - 5th order Whitney form
			testFunction = [&](Vector2 p) -> double {
				return p.x * p.x * p.x * p.y * p.y / 8.0 - p.y * p.y * p.y * p.y / 2.0;
			};
			sourceTerm = [&](Vector2 p) -> double {
				return (3.0 * p.x * p.y * p.y - p.x * p.x * p.x) / 4.0 + 6.0 * p.y * p.y;
			};
			testFunctionGradient = [&](Vector2 p) -> Vector2 {
				return Vector2((3.0 / 8.0) * p.x * p.x * p.y * p.y, p.x * p.x * p.x * p.y / 4.0 - 2.0 * p.y * p.y * p.y);
			};
		}
		else if (testFunc == 'c') { //testFunction c) - 10th order Whitney form
			testFunction = [&](Vector2 p) -> double {
				return p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y * p.y * p.y / 2048.0 + p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y / 64.0;
			};
			sourceTerm = [&](Vector2 p) -> double {
				return (42.0 / 2048.0) * p.x * p.x * p.x * p.x * p.x * p.y * p.y * p.y - (6.0 / 2048.0) * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y
					- (90.0 / 64.0) * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y;
			};
			testFunctionGradient = [&](Vector2 p) -> Vector2 {
				return Vector2((7.0 / 2048.0) * p.x * p.x * p.x * p.x * p.x * p.x * p.y * p.y * p.y,
					(3.0 / 2048.0) * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.y * p.y + (10.0 / 64.0) * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y);
			};
		}
		else {
			std::cout << "Choose polynomialWaveExampleWhitney2D a), b), or c)" << '\n';
			return;
		}

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text filepath;
		filepath << "Files/polynomialWaveExampleWhitney2D_" << testFunc << ".txt";

		//Create triangle mesh in [0, 2] x [0, 2].
		BuilderMesh mesh_old(2);
		createSpacetimeMesh(mesh_old, L, T, 2, 2);

		//mark boundary nodes to indicate boundary conditions: flag is 2 for nodes in ]0,2[ x {0}, 3 for nodes in ]0,2[ x {2}, and 1 for other boundary nodes
		mesh_old.fillBoundaryFlags(1);
		for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
			if (mesh_old.getNodeFlag(i) == 1 && mesh_old.getNodePosition2(i).y < 1e-14 && mesh_old.getNodePosition2(i).x > 1e-14 && mesh_old.getNodePosition2(i).x < L - 1e-14)
				mesh_old.setNodeFlag(i, 2);
			else if (mesh_old.getNodeFlag(i) == 1 && mesh_old.getNodePosition2(i).y > T - 1e-14 && mesh_old.getNodePosition2(i).x > 1e-14 && mesh_old.getNodePosition2(i).x < L - 1e-14)
				mesh_old.setNodeFlag(i, 3);
		}

		//solve the problem using Whitney forms of orders 1-10
		for (uint order = 1; order <= 10; ++order) {

			//refine the mesh into small simplices
			BuilderMesh mesh(2);
			SmallSimplexPartition2D ssp(order);
			ssp.refineMesh(mesh_old, mesh);

			//similarly mark boundary nodes of refined mesh
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 && mesh.getNodePosition2(i).y < 1e-14 && mesh.getNodePosition2(i).x > 1e-14 && mesh.getNodePosition2(i).x < L - 1e-14)
					mesh.setNodeFlag(i, 2);
				else if (mesh.getNodeFlag(i) == 1 && mesh.getNodePosition2(i).y > T - 1e-14 && mesh.getNodePosition2(i).x > 1e-14 && mesh.getNodePosition2(i).x < L - 1e-14)
					mesh.setNodeFlag(i, 3);
			}

			//compute discrete solution using Hodge operator of given order
			Buffer<double> dec_sol;
			getDECSolutionWaveEq(mesh_old, mesh, ssp, sourceTerm, testFunction, testFunctionGradient, dec_sol);

			//error analysis
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//maximum and average error
			double maxError = 0.0;
			double averageError = 0.0;
			uint interiorNodes = 0;
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
					continue;
				++interiorNodes;
				double diff = std::abs(dec_sol[i] - testFunction(mesh.getNodePosition2(i)));
				if (maxError < diff)
					maxError = diff;
				averageError += diff;
			}
			averageError /= interiorNodes;
			std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

			//L2 error
			double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, dec_sol, testFunction, 1000);
			std::cout << "L2 error: " << l2error1 << '\n';
			output << "L2 error: " << l2error1 << '\n';

			//H1 error
			Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
			for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
				const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
				for (uint j = 0; j < nodes.size(); ++j)
					dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * dec_sol[nodes[j]];
			}
			double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
			std::cout << "Gradient L2 error: " << graderror1 << '\n';
			output << "Gradient L2 error: " << graderror1 << '\n';
			double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
			std::cout << "H1 error: " << h1error1 << '\n' << '\n';
			output << "H1 error: " << h1error1 << '\n' << '\n';
			output.save(filepath.str());
		}
	}

	/*
	Solves inhomogenous wave equation in domain [0, L] x [0, T] with L=2, T=2, using higher order Whitney forms and discrete Hodge operators in spacetime.
	The exact solution is chosen to be exp(x/4)sin(x+t+2)cos(t).
	This implies the boundary condition and the source term.
	Here we solve the whole system at once; see waveExampleWhitney2D_TimeStepping for time stepping.
	*/
	void waveExampleWhitney2D() {
		const double L = 2.0;
		const double T = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector2)> testFunction{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * std::sin(p.x + p.y + 2) * std::cos(p.y); } };
		std::function<double(Vector2)> sourceTerm{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * (0.5 * std::cos(p.y) * std::cos(2 + p.x + p.y) + 2 * std::cos(2 + p.x + p.y) * std::sin(p.y)
			+ 17.0 / 16 * std::cos(p.y) * std::sin(2 + p.x + p.y)); } };
		std::function<Vector2(Vector2)> testFunctionGradient{ [&](Vector2 p) -> Vector2 { return Vector2(std::exp(p.x / 4) * std::cos(p.y) * (std::sin(p.x + p.y + 2) / 4 + std::cos(p.x + p.y + 2)),
			std::exp(p.x / 4) * (std::cos(p.y) * std::cos(p.x + p.y + 2) - std::sin(p.y) * std::sin(p.x + p.y + 2))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 8; ++meshNumber) { //solve the problem on meshes of different grain
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//Create triangle mesh in [0, 2] x [0, 2].
			BuilderMesh mesh_old(2);
			uint refs = std::pow(2, meshNumber - 1);
			createSpacetimeMesh(mesh_old, L, T, refs, refs);

			//mark boundary nodes to indicate boundary conditions: flag is 2 for nodes in ]0,2[ x {0}, 3 for nodes in ]0,2[ x {2}, and 1 for other boundary nodes
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				if (mesh_old.getNodeFlag(i) == 1 && mesh_old.getNodePosition2(i).y < eps && mesh_old.getNodePosition2(i).x > eps && mesh_old.getNodePosition2(i).x < L - eps)
					mesh_old.setNodeFlag(i, 2);
				else if (mesh_old.getNodeFlag(i) == 1 && mesh_old.getNodePosition2(i).y > T - eps && mesh_old.getNodePosition2(i).x > eps && mesh_old.getNodePosition2(i).x < L - eps)
					mesh_old.setNodeFlag(i, 3);
			}

			//solve the problem using Whitney forms of orders 1-6
			for (uint order = 1; order <= 6 && (meshNumber <= 5 || order <= 5) && (meshNumber <= 6 || order <= 2) && (meshNumber <= 7 || order <= 1); ++order) {

				//refine the mesh into small simplices
				BuilderMesh mesh(2);
				SmallSimplexPartition2D ssp(order);
				ssp.refineMesh(mesh_old, mesh);

				//similarly mark boundary nodes of refined mesh
				mesh.fillBoundaryFlags(1);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 && mesh.getNodePosition2(i).y < eps && mesh.getNodePosition2(i).x > eps && mesh.getNodePosition2(i).x < L - eps)
						mesh.setNodeFlag(i, 2);
					else if (mesh.getNodeFlag(i) == 1 && mesh.getNodePosition2(i).y > T - eps && mesh.getNodePosition2(i).x > eps && mesh.getNodePosition2(i).x < L - eps)
						mesh.setNodeFlag(i, 3);
				}

				//compute discrete solution using Hodge operator of given order
				Buffer<double> dec_sol;
				getDECSolutionWaveEq(mesh_old, mesh, ssp, sourceTerm, testFunction, testFunctionGradient, dec_sol);

				//error analysis
				std::cout << "DEC solution of order " << order << '\n';
				output << "DEC solution of order " << order << '\n';

				//maximum and average error
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(dec_sol[i] - testFunction(mesh.getNodePosition2(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageError /= interiorNodes;
				std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
				output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

				//L2 error and H1 error
				double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, dec_sol, testFunction, 1000);
				std::cout << "L2 error: " << l2error1 << '\n';
				output << "L2 error: " << l2error1 << '\n';
				Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
					for (uint j = 0; j < nodes.size(); ++j)
						dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * dec_sol[nodes[j]];
				}
				double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
				std::cout << "Gradient L2 error: " << graderror1 << '\n';
				output << "Gradient L2 error: " << graderror1 << '\n';
				double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
				std::cout << "H1 error: " << h1error1 << '\n';
				output << "H1 error: " << h1error1 << '\n';
				std::cout << '\n';
				output << '\n';
				output.save("Files/waveExampleWhitney2D.txt");
			}
		}
	}

	/*
	Solves inhomogenous wave equation in domain [0, Lx] x [0, Ly] x [0, T] with Lx=2, Ly=2, and T=2, using higher order Whitney forms and discrete Hodge operators in spacetime.
	The exact solution is the test function a), b), or c) chosen with parameter testFunc.
	This implies boundary conditions and the source term.
	*/
	void polynomialWaveExampleWhitney3D(const char testFunc) {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double T = 2.0;
		std::function<double(Vector3)> testFunction;
		std::function<double(Vector3)> sourceTerm;
		std::function<Vector3(Vector3)> testFunctionGradient;
		if (testFunc == 'a') { //testFunction a) - 1st order Whitney form
			testFunction = [&](Vector3 p) -> double {
				return 2.5 * p.x - p.y + 1.25 * p.z;
			};
			sourceTerm = [&](Vector3 p) -> double {
				return 0.0;
			};
			testFunctionGradient = [&](Vector3 p) -> Vector3 {
				return Vector3(2.5, -1.0, 1.25);
			};
		}
		else if (testFunc == 'b') { //testFunction b) - 5th order Whitney form
			testFunction = [&](Vector3 p) -> double {
				return 10 * p.x * p.x * p.y * p.y * p.y - 40 * p.z * p.z * p.z * p.z * p.z;
			};
			sourceTerm = [&](Vector3 p) -> double {
				return 20 * p.y * p.y * p.y + 60 * p.x * p.x * p.y + 800 * p.z * p.z * p.z;
			};
			testFunctionGradient = [&](Vector3 p) -> Vector3 {
				return Vector3(20 * p.x * p.y * p.y * p.y, 30 * p.x * p.x * p.y * p.y, -200 * p.z * p.z * p.z * p.z);
			};
		}
		else if (testFunc == 'c') { //testFunction c) - 10th order Whitney form
			testFunction = [&](Vector3 p) -> double {
				return 20 * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x - 1000 * p.x * p.y * p.y * p.y * p.y * p.y * p.y * p.z * p.z * p.z;
			};
			sourceTerm = [&](Vector3 p) -> double {
				return 1800 * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x - 30000 * p.x * p.y * p.y * p.y * p.y * p.z * p.z * p.z + 6000 * p.x * p.y * p.y * p.y * p.y * p.y * p.y * p.z;
			};
			testFunctionGradient = [&](Vector3 p) -> Vector3 {
				return Vector3(200 * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x - 1000 * p.y * p.y * p.y * p.y * p.y * p.y * p.z * p.z * p.z,
					-6000 * p.x * p.y * p.y * p.y * p.y * p.y * p.z * p.z * p.z, -3000 * p.x * p.y * p.y * p.y * p.y * p.y * p.y * p.z * p.z);
			};
		}
		else {
			std::cout << "Choose polynomialWaveExampleWhitney3D a), b), or c)" << '\n';
			return;
		}

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text filepath;
		filepath << "Files/polynomialWaveExampleWhitney3D_" << testFunc << ".txt";

		//Create tetrahedral mesh in [0, 2] x [0, 2] x [0, 2].
		BuilderMesh mesh_old(3);
		createSpacetimeMesh(mesh_old, Lx, Ly, T, 2, 2, 2);

		//mark boundary nodes to indicate boundary conditions: flag is 2 for nodes in ]0,2[ x ]0,2[ x {0}, 3 for nodes in ]0,2[ x ]0,2[ x {2}, and 1 for other boundary nodes
		mesh_old.fillBoundaryFlags(1);
		for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
			if (mesh_old.getNodeFlag(i) == 1 && mesh_old.getNodePosition3(i).z < 1e-14 && mesh_old.getNodePosition3(i).x > 1e-14 && mesh_old.getNodePosition3(i).x < Lx - 1e-14
				&& mesh_old.getNodePosition3(i).y > 1e-14 && mesh_old.getNodePosition3(i).y < Ly - 1e-14)
				mesh_old.setNodeFlag(i, 2);
			else if (mesh_old.getNodeFlag(i) == 1 && mesh_old.getNodePosition3(i).z > T - 1e-14 && mesh_old.getNodePosition3(i).x > 1e-14 && mesh_old.getNodePosition3(i).x < Lx - 1e-14
				&& mesh_old.getNodePosition3(i).y > 1e-14 && mesh_old.getNodePosition3(i).y < Ly - 1e-14)
				mesh_old.setNodeFlag(i, 3);
		}

		//solve the problem using Whitney forms of orders 1-10
		for (uint order = 1; order <= 10; ++order) {

			//refine the mesh into small simplices
			BuilderMesh mesh(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh_old, mesh);

			//similarly mark boundary nodes of refined mesh
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 && mesh.getNodePosition3(i).z < 1e-14 && mesh.getNodePosition3(i).x > 1e-14 && mesh.getNodePosition3(i).x < Lx - 1e-14
					&& mesh.getNodePosition3(i).y > 1e-14 && mesh.getNodePosition3(i).y < Ly - 1e-14)
					mesh.setNodeFlag(i, 2);
				else if (mesh.getNodeFlag(i) == 1 && mesh.getNodePosition3(i).z > T - 1e-14 && mesh.getNodePosition3(i).x > 1e-14 && mesh.getNodePosition3(i).x < Lx - 1e-14
					&& mesh.getNodePosition3(i).y > 1e-14 && mesh.getNodePosition3(i).y < Ly - 1e-14)
					mesh.setNodeFlag(i, 3);
			}

			//compute discrete solution using Hodge operator of given order
			Buffer<double> dec_sol;
			getDECSolutionWaveEq(mesh_old, mesh, ssp, sourceTerm, testFunction, testFunctionGradient, dec_sol);

			//error analysis
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//maximum and average error
			double maxError = 0.0;
			double averageError = 0.0;
			uint interiorNodes = 0;
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
					continue;
				++interiorNodes;
				double diff = std::abs(dec_sol[i] - testFunction(mesh.getNodePosition3(i)));
				if (maxError < diff)
					maxError = diff;
				averageError += diff;
			}
			averageError /= interiorNodes;
			std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

			//L2 error
			double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, dec_sol, testFunction, 1000);
			std::cout << "L2 error: " << l2error1 << '\n';
			output << "L2 error: " << l2error1 << '\n';

			//H1 error
			Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
			for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
				const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
				for (uint j = 0; j < nodes.size(); ++j)
					dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * dec_sol[nodes[j]];
			}
			double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
			std::cout << "Gradient L2 error: " << graderror1 << '\n';
			output << "Gradient L2 error: " << graderror1 << '\n';
			double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
			std::cout << "H1 error: " << h1error1 << '\n' << '\n';
			output << "H1 error: " << h1error1 << '\n' << '\n';

			output.save(filepath.str());
		}
	}

	/*
	Solves inhomogenous wave equation in domain [0, Lx] x [0, Ly] x [0, T] with Lx=2, Ly=2, and T=2, using higher order Whitney forms and discrete Hodge operators in spacetime.
	The exact solution is chosen to be exp(x/4)sin(x+y+t+1)cos(y)(2-0.003t-0.00002t^2).
	This implies boundary conditions and the source term.
	Here we solve the whole system at once; see waveExampleWhitney3D_TimeStepping for time stepping.
	*/
	void waveExampleWhitney3D() {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double T = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector3)> testFunction{ [&](Vector3 p) -> double { return std::exp(0.25 * p.x) * std::sin(p.x + p.y + p.z + 1) * std::cos(p.y) * (2 - 0.003 * p.z - 0.00002 * p.z * p.z); } };
		std::function<double(Vector3)> sourceTerm{ [&](Vector3 p) -> double { return -0.00001 * std::exp(0.25 * p.x) * (((-4 * p.z - 600) * p.z + 400000) * std::sin(p.y) * std::cos(p.z + p.x + p.y + 1)
			+ std::cos(p.y) * (((-3.875 * p.z - 581.25) * p.z + 387496) * std::sin(p.z + p.x + p.y + 1) + (p.z * (p.z + 142) - 100600) * std::cos(p.z + p.x + p.y + 1))); } };
		std::function<Vector3(Vector3)> testFunctionGradient{ [&](Vector3 p) -> Vector3 { return std::exp(0.25 * p.x) * Vector3(0.25 * (2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.y) * (4 * std::cos(p.x + p.y + p.z + 1) + std::sin(p.x + p.y + p.z + 1)),
			(2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + 2 * p.y + p.z + 1), std::cos(p.y) * ((2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + p.y + p.z + 1) - (0.003 + 0.00004 * p.z) * std::sin(p.x + p.y + p.z + 1))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		for (uint meshNumber = 1; meshNumber <= 6; ++meshNumber) { //solve the problem on meshes of different grain
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//Create tetrahedral mesh in [0, 2] x [0, 2] x [0, 2].
			BuilderMesh mesh_old(3);
			uint refs = std::pow(2, meshNumber - 1);
			createSpacetimeMesh(mesh_old, Lx, Ly, T, refs, refs, refs);

			//mark boundary nodes to indicate boundary conditions: flag is 2 for nodes in ]0,2[ x ]0,2[ x {0}, 3 for nodes in ]0,2[ x ]0,2[ x {2}, and 1 for other boundary nodes
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				Vector3 nodePos = mesh_old.getNodePosition3(i);
				if (mesh_old.getNodeFlag(i) == 1 && nodePos.z < eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
					mesh_old.setNodeFlag(i, 2);
				else if (mesh_old.getNodeFlag(i) == 1 && nodePos.z > T - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
					mesh_old.setNodeFlag(i, 3);
			}

			//solve the problem using Whitney forms of orders 1-5
			for (uint order = 1; order <= 5 && (meshNumber <= 3 || order <= 4) && (meshNumber <= 4 || order <= 2) && (meshNumber <= 5 || order <= 1); ++order) {

				//refine the mesh into small simplices
				BuilderMesh mesh(3);
				SmallSimplexPartition3D ssp(order);
				ssp.refineMesh(mesh_old, mesh);

				//similarly mark boundary nodes of refined mesh
				mesh.fillBoundaryFlags(1);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					Vector3 nodePos = mesh.getNodePosition3(i);
					if (mesh.getNodeFlag(i) == 1 && nodePos.z < eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
						mesh.setNodeFlag(i, 2);
					else if (mesh.getNodeFlag(i) == 1 && nodePos.z > T - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
						mesh.setNodeFlag(i, 3);
				}

				//compute discrete solution using Hodge operator of given order
				Buffer<double> dec_sol;
				getDECSolutionWaveEq(mesh_old, mesh, ssp, sourceTerm, testFunction, testFunctionGradient, dec_sol);

				//error analysis
				std::cout << "DEC solution of order " << order << '\n';
				output << "DEC solution of order " << order << '\n';

				//maximum and average error
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(dec_sol[i] - testFunction(mesh.getNodePosition3(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageError /= interiorNodes;
				std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
				output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

				//L2 error and H1 error
				double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, dec_sol, testFunction, 1000);
				std::cout << "L2 error: " << l2error1 << '\n';
				output << "L2 error: " << l2error1 << '\n';
				Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
					for (uint j = 0; j < nodes.size(); ++j)
						dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * dec_sol[nodes[j]];
				}
				double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
				std::cout << "Gradient L2 error: " << graderror1 << '\n';
				output << "Gradient L2 error: " << graderror1 << '\n';
				double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
				std::cout << "H1 error: " << h1error1 << '\n';
				output << "H1 error: " << h1error1 << '\n';
				std::cout << '\n';
				output << '\n';

				output.save("Files/waveExampleWhitney3D.txt");
			}
		}
	}

	/*
	Solves inhomogenous wave equation in domain [0, L] x [0, T * T_REPEATS], using higher order Whitney forms and discrete Hodge operators in spacetime.
	The exact solution is chosen to be exp(x/4)sin(x+t+2)cos(t).
	This implies boundary conditions and the source term.
	Instead of solving the full system, we form a system for a single time step (of length T) and update the solution one time step at a time (for T_REPEATS steps).
	*/
	void waveExampleWhitney2D_TimeStepping(uint order, const double T, const uint T_REPEATS) {
		const double L = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector2)> testFunction{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * std::sin(p.x + p.y + 2) * std::cos(p.y); } };
		std::function<double(Vector2)> sourceTerm{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * (0.5 * std::cos(p.y) * std::cos(2 + p.x + p.y) + 2 * std::cos(2 + p.x + p.y) * std::sin(p.y)
			+ 17.0 / 16 * std::cos(p.y) * std::sin(2 + p.x + p.y)); } };
		std::function<Vector2(Vector2)> testFunctionGradient{ [&](Vector2 p) -> Vector2 { return Vector2(std::exp(p.x / 4) * std::cos(p.y) * (std::sin(p.x + p.y + 2) / 4 + std::cos(p.x + p.y + 2)),
			std::exp(p.x / 4) * (std::cos(p.y) * std::cos(p.x + p.y + 2) - std::sin(p.y) * std::sin(p.x + p.y + 2))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleWhitney2D_TimeStepping_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		uint maxMeshNumber = 8;
		if (order == 2)
			maxMeshNumber = 7;
		else if (order == 3 || order == 4 || order == 5)
			maxMeshNumber = 6;
		else if (order >= 6)
			maxMeshNumber = 5;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) { //solve the problem on meshes of different grain
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;

			//solve the problem using Whitney forms of chosen order
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//save maximum error and maximum stepwise average
			double maxErrorFull = 0.0;
			double maxAverage = 0.0;
			double averageSum = 0.0;
			uint averageNum = 0;

			//Create the mesh covering one time step only and its refinement and mark cells with appropriate flags.
			BuilderMesh mesh_old(2);
			createSpacetimeMesh(mesh_old, L, dtime, refs, 1);
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				Vector2 nodePos = mesh_old.getNodePosition2(i);
				if (nodePos.y > dtime - eps && nodePos.x > eps && nodePos.x < L - eps)
					mesh_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
				Vector2 edgePos = mesh_old.getEdgeAverage2(i);
				if (edgePos.y > dtime - eps && edgePos.x > eps && edgePos.x < L - eps)
					mesh_old.setEdgeFlag(i, 0);
			}
			BuilderMesh mesh(2);
			SmallSimplexPartition2D ssp(order);
			ssp.refineMesh(mesh_old, mesh);
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1) {
					Vector2 nodePos = mesh.getNodePosition2(i);
					if (nodePos.y < eps && nodePos.x > eps && nodePos.x < L - eps) {
						mesh.setNodeFlag(i, 2);
					}
					else if (nodePos.y > dtime - eps && nodePos.x > eps && nodePos.x < L - eps)
						mesh.setNodeFlag(i, 3);
				}
			}

			//form the system matrix
			Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
			{
				Buffer<std::unordered_map<uint, double>> star;
				Buffer<uint> refElementIndices(2);
				Buffer <Buffer<uint>> refElements(2);
				for (uint i = 0; i < 2; ++i) {
					refElementIndices[i] = i;
					refElements[i].resize(mesh_old.getFaceSize() / 2);
				}
				for (uint i = 0; i < mesh_old.getFaceSize(); ++i)
					refElements[i % 2][i / 2] = i;
				ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "RightTriangle", false, true, false);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//restrict to active nodes
			uint boundary2Nodes = order * refs - 1;
			uint activeNodes = order * boundary2Nodes;
			Buffer<Buffer<double>> system_matrix_interior(activeNodes);
			for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
					continue;
				std::unordered_map<uint, double>& row = system_matrix[i];
				system_matrix_interior[ii] = Buffer<double>(activeNodes, 0.0);
				for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
					if (mesh.getNodeFlag(j) == 1 || mesh.getNodeFlag(j) == 2)
						continue;
					auto it = row.find(j);
					if (it != row.end() && std::abs(it->second) > 1.0e-14) {
						system_matrix_interior[ii][jj] = it->second;
					}
					++jj;
				}
				++ii;
			}
			Buffer<uint> p;
			decomposeLUP(system_matrix_interior, p, 1e-14);

			//initial condition for update step
			Buffer<double> rhs_prepared(boundary2Nodes);
			Buffer<double> initial_values(boundary2Nodes);
			for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 2) {
					initial_values[j] = testFunction(mesh.getNodePosition2(i));
					rhs_prepared[j] = -integrateNodeDualBoundary(mesh, i, testFunctionGradient, false, true);
					++j;
				}
			}

			uint stepNumber = 1;
			while (true) {
				//form the right-hand side
				Buffer<double> rhs(activeNodes);
				Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						boundaryValues[i] = testFunction(mesh.getNodePosition2(i));
					else if (mesh.getNodeFlag(i) == 2) {
						boundaryValues[i] = initial_values[j++];
					}
				}
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
						continue;
					rhs[j] = integrateDual2Cell(mesh, i, sourceTerm);
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						rhs[j] -= it->second * boundaryValues[it->first];
					}
					if (mesh.getNodeFlag(i) == 2)
						rhs[j] += rhs_prepared[k++];
					++j;
				}

				//get the solution on active nodes by solving the system and on other nodes from the boundary condition
				Buffer<double> sol_interior;
				solveLUP(system_matrix_interior, p, rhs, sol_interior);
				Buffer<double> sol;
				sol.resize(mesh.getNodeSize());
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						sol[i] = testFunction(mesh.getNodePosition2(i));
					else if (mesh.getNodeFlag(i) == 2)
						sol[i] = initial_values[k++];
					else
						sol[i] = sol_interior[j++];
				}

				//stepwise error analysis
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(sol[i] - testFunction(mesh.getNodePosition2(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageSum += averageError;
				averageNum += interiorNodes;
				averageError /= interiorNodes;
				if (stepNumber % 10 == 0) { //choose how often to print the result
					std::cout << "Time elapsed: " << stepNumber * dtime << '\n';
					output << "Time elapsed: " << stepNumber * dtime << '\n';
					std::cout << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					output << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
				}
				if (maxError > maxErrorFull)
					maxErrorFull = maxError;
				if (averageError > maxAverage)
					maxAverage = averageError;

				//stop if done
				if (++stepNumber > T_REPEATS * refs)
					break;

				//prepare initial condition for the next step
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 3) {
						rhs_prepared[j] = integrateDual2Cell(mesh, i, sourceTerm);
						std::unordered_map<uint, double>& row = system_matrix[i];
						for (auto it = row.begin(); it != row.end(); ++it) {
							rhs_prepared[j] -= it->second * sol[it->first];
						}
						initial_values[j] = sol[i];
						++j;
					}
				}

				//move the mesh
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					mesh_old.setNodePosition(i, mesh_old.getNodePosition(i) + Vector4(0.0, dtime, 0.0, 0.0));
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(0.0, dtime, 0.0, 0.0));
				}
			}

			//error analysis for the full solution
			std::cout << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";
			output << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";

			//maximum and average error
			averageSum /= averageNum;
			std::cout << "Max error " << maxErrorFull << " (average " << averageSum << ", maximum stepwise average " << maxAverage << ')' << '\n';
			output << "Max error " << maxErrorFull << " (average " << averageSum << ", maximum stepwise average " << maxAverage << ')' << '\n';
			std::cout << '\n';
			output << '\n';

			output.save(outputFile.str());
		}
	}

	//same as waveExampleWhitney2D_TimeStepping but saves the full solution for the full mesh
	void waveExampleWhitney2D_TimeSteppingFullSolution(uint order, const double T, const uint T_REPEATS) {
		const double L = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector2)> testFunction{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * std::sin(p.x + p.y + 2) * std::cos(p.y); } };
		std::function<double(Vector2)> sourceTerm{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * (0.5 * std::cos(p.y) * std::cos(2 + p.x + p.y) + 2 * std::cos(2 + p.x + p.y) * std::sin(p.y)
			+ 17.0 / 16 * std::cos(p.y) * std::sin(2 + p.x + p.y)); } };
		std::function<Vector2(Vector2)> testFunctionGradient{ [&](Vector2 p) -> Vector2 { return Vector2(std::exp(p.x / 4) * std::cos(p.y) * (std::sin(p.x + p.y + 2) / 4 + std::cos(p.x + p.y + 2)),
			std::exp(p.x / 4) * (std::cos(p.y) * std::cos(p.x + p.y + 2) - std::sin(p.y) * std::sin(p.x + p.y + 2))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleWhitney2D_TimeSteppingFullSolution_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		uint maxMeshNumber = 8;
		if (order == 2)
			maxMeshNumber = 7;
		else if (order == 3 || order == 4 || order == 5)
			maxMeshNumber = 6;
		else if (order >= 6)
			maxMeshNumber = 5;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) { //solve the problem on meshes of different grain
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//Create the full mesh in [0, L] x [0, T * T_REPEATS] and mark cells with appropriate flags.
			uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;
			BuilderMesh meshFull_old(2);
			createSpacetimeMesh(meshFull_old, L, T * T_REPEATS, refs, refs * T_REPEATS);
			meshFull_old.fillBoundaryFlags(1);
			for (uint i = 0; i < meshFull_old.getNodeSize(); ++i) {
				Vector2 nodePos = meshFull_old.getNodePosition2(i);
				if (nodePos.y > T_REPEATS * T - eps && nodePos.x > eps && nodePos.x < L - eps)
					meshFull_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < meshFull_old.getEdgeSize(); ++i) {
				Vector2 edgePos = meshFull_old.getEdgeAverage2(i);
				if (edgePos.y > T_REPEATS * T - eps && edgePos.x > eps && edgePos.x < L - eps)
					meshFull_old.setEdgeFlag(i, 0);
			}

			//solve the problem using Whitney forms of chosen order
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//Refine the full mesh into small simplices and mark cells with appropriate flags.
			BuilderMesh meshFull(2);
			SmallSimplexPartition2D sspFull(order);
			sspFull.refineMesh(meshFull_old, meshFull);
			meshFull.fillBoundaryFlags(1);
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1) {
					Vector2 nodePos = meshFull.getNodePosition2(i);
					if (nodePos.y < eps && nodePos.x > eps && nodePos.x < L - eps) {
						meshFull.setNodeFlag(i, 2);
					}
					else if (nodePos.y > T_REPEATS * T - eps && nodePos.x > eps && nodePos.x < L - eps)
						meshFull.setNodeFlag(i, 3);
				}
			}

			//prepare memory and variables to store values in the full solution on each time step and set boundary values using boundary condition
			uint nodesUpdated = 0;
			uint edgesUpdated = 0;
			uint facesUpdated = 0;
			const uint smallNodesInFaces = binCoef(order - 1, 2);
			Buffer<double> solFull(meshFull.getNodeSize());
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1 || meshFull.getNodeFlag(i) == 2) {
					solFull[i] = testFunction(meshFull.getNodePosition2(i));
				}
			}

			//Create the mesh covering one time step only and its refinement and mark cells with appropriate flags.
			BuilderMesh mesh_old(2);
			createSpacetimeMesh(mesh_old, L, dtime, refs, 1);
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				Vector2 nodePos = mesh_old.getNodePosition2(i);
				if (nodePos.y > dtime - eps && nodePos.x > eps && nodePos.x < L - eps)
					mesh_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
				Vector2 edgePos = mesh_old.getEdgeAverage2(i);
				if (edgePos.y > dtime - eps && edgePos.x > eps && edgePos.x < L - eps)
					mesh_old.setEdgeFlag(i, 0);
			}
			BuilderMesh mesh(2);
			SmallSimplexPartition2D ssp(order);
			ssp.refineMesh(mesh_old, mesh);
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1) {
					Vector2 nodePos = mesh.getNodePosition2(i);
					if (nodePos.y < eps && nodePos.x > eps && nodePos.x < L - eps) {
						mesh.setNodeFlag(i, 2);
					}
					else if (nodePos.y > dtime - eps && nodePos.x > eps && nodePos.x < L - eps)
						mesh.setNodeFlag(i, 3);
				}
			}

			//form the system matrix
			Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
			{
				Buffer<std::unordered_map<uint, double>> star;
				Buffer<uint> refElementIndices(1, 0);
				Buffer <Buffer<uint>> refElements(1);
				refElements[0].resize(mesh_old.getFaceSize());
				for (uint i = 0; i < mesh_old.getFaceSize(); ++i)
					refElements[0][i] = i;
				ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "RightTriangle", false, true, false);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//restrict to active nodes
			uint boundary2Nodes = order * refs - 1;
			uint activeNodes = order * boundary2Nodes;
			Buffer<Buffer<double>> system_matrix_interior(activeNodes);
			for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
					continue;
				std::unordered_map<uint, double>& row = system_matrix[i];
				system_matrix_interior[ii] = Buffer<double>(activeNodes, 0.0);
				for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
					if (mesh.getNodeFlag(j) == 1 || mesh.getNodeFlag(j) == 2)
						continue;
					auto it = row.find(j);
					if (it != row.end() && std::abs(it->second) > 1.0e-14) {
						system_matrix_interior[ii][jj] = it->second;
					}
					++jj;
				}
				++ii;
			}
			Buffer<uint> p;
			decomposeLUP(system_matrix_interior, p, 1e-14);

			//initial condition for update step
			Buffer<double> rhs_prepared(boundary2Nodes);
			Buffer<double> initial_values(boundary2Nodes);
			for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 2) {
					initial_values[j] = testFunction(mesh.getNodePosition2(i));
					rhs_prepared[j] = -integrateNodeDualBoundary(mesh, i, testFunctionGradient, false, true);
					++j;
				}
			}

			uint stepNumber = 1;
			while (true) {
				//form the right-hand side
				Buffer<double> rhs(activeNodes);
				Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						boundaryValues[i] = testFunction(mesh.getNodePosition2(i));
					else if (mesh.getNodeFlag(i) == 2) {
						boundaryValues[i] = initial_values[j++];
					}
				}
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
						continue;
					rhs[j] = integrateDual2Cell(mesh, i, sourceTerm);
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						rhs[j] -= it->second * boundaryValues[it->first];
					}
					if (mesh.getNodeFlag(i) == 2)
						rhs[j] += rhs_prepared[k++];
					++j;
				}

				//get the solution on active nodes by solving the system and on other nodes from the boundary condition
				Buffer<double> sol_interior;
				solveLUP(system_matrix_interior, p, rhs, sol_interior);
				Buffer<double> sol;
				sol.resize(mesh.getNodeSize());
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						sol[i] = testFunction(mesh.getNodePosition2(i));
					else if (mesh.getNodeFlag(i) == 2)
						sol[i] = initial_values[k++];
					else
						sol[i] = sol_interior[j++];
				}

				//store the values in full solution
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					if (mesh_old.getNodeFlag(i) == 1)
						continue;
					while (meshFull_old.getNodeFlag(nodesUpdated) == 1)
						++nodesUpdated;
					solFull[nodesUpdated++] = sol[i];
				}
				for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
					if (mesh_old.getEdgeFlag(i) == 1)
						continue;
					while (meshFull_old.getEdgeFlag(edgesUpdated) == 1)
						++edgesUpdated;
					for (uint j = 0; j < order - 1; ++j)
						solFull[sspFull.smallNodesEdgeList(edgesUpdated, j)] = sol[ssp.smallNodesEdgeList(i, j)];
					++edgesUpdated;
				}
				for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
					for (uint j = 0; j < smallNodesInFaces; ++j)
						solFull[sspFull.smallNodesFaceList(facesUpdated, j)] = sol[ssp.smallNodesFaceList(i, j)];
					++facesUpdated;
				}

				//stepwise error analysis				
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(sol[i] - testFunction(mesh.getNodePosition2(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageError /= interiorNodes;
				if (stepNumber % 10 == 0) { //choose how often to print the result
					std::cout << "Time elapsed: " << stepNumber * dtime << '\n';
					output << "Time elapsed: " << stepNumber * dtime << '\n';
					std::cout << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					output << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					/*double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, sol, testFunction, 1000);
					std::cout << "L2 error: " << l2error1 << '\n';
					output << "L2 error: " << l2error1 << '\n';
					Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
						for (uint j = 0; j < nodes.size(); ++j)
							dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * sol[nodes[j]];
					}
					double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
					std::cout << "Gradient L2 error: " << graderror1 << '\n';
					output << "Gradient L2 error: " << graderror1 << '\n';
					double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
					std::cout << "H1 error: " << h1error1 << '\n';
					output << "H1 error: " << h1error1 << '\n';*/
				}

				//stop if done
				if (++stepNumber > T_REPEATS * refs)
					break;

				//prepare initial condition for the next step
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 3) {
						rhs_prepared[j] = integrateDual2Cell(mesh, i, sourceTerm);
						std::unordered_map<uint, double>& row = system_matrix[i];
						for (auto it = row.begin(); it != row.end(); ++it) {
							rhs_prepared[j] -= it->second * sol[it->first];
						}
						initial_values[j] = sol[i];
						++j;
					}
				}

				//move the mesh
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					mesh_old.setNodePosition(i, mesh_old.getNodePosition(i) + Vector4(0.0, dtime, 0.0, 0.0));
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(0.0, dtime, 0.0, 0.0));
				}
			}

			//error analysis for the full solution
			std::cout << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";
			output << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";

			//maximum and average error
			double maxError = 0.0;
			double averageError = 0.0;
			uint interiorNodes = 0;
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1 || meshFull.getNodeFlag(i) == 2)
					continue;
				++interiorNodes;
				double diff = std::abs(solFull[i] - testFunction(meshFull.getNodePosition2(i)));
				if (maxError < diff)
					maxError = diff;
				averageError += diff;
			}
			averageError /= interiorNodes;
			std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

			//L2 error and H1 error
			/*double l2error1 = computeDifferenceL2Norm0Form(meshFull_old, sspFull, solFull, testFunction, 1000);
			std::cout << "L2 error: " << l2error1 << '\n';
			output << "L2 error: " << l2error1 << '\n';
			Buffer<double> dec_solFull_grad(meshFull.getEdgeSize(), 0.0);
			for (uint i = 0; i < meshFull.getEdgeSize(); ++i) {
				const Buffer<uint>& nodes = meshFull.getEdgeNodes(i);
				for (uint j = 0; j < nodes.size(); ++j)
					dec_solFull_grad[i] += meshFull.getEdgeIncidence(i, nodes[j]) * solFull[nodes[j]];
			}
			double graderror1 = computeDifferenceL2Norm1Form(meshFull_old, sspFull, dec_solFull_grad, testFunctionGradient, 1000);
			std::cout << "Gradient L2 error: " << graderror1 << '\n';
			output << "Gradient L2 error: " << graderror1 << '\n';
			double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
			std::cout << "H1 error: " << h1error1 << '\n';
			output << "H1 error: " << h1error1 << '\n';*/

			std::cout << '\n';
			output << '\n';
			output.save(outputFile.str());
		}
	}

	/*
	Solves inhomogenous wave equation in domain [0, Lx] x [0, Ly] x [0, T * T_REPEATS], using higher order Whitney forms and discrete Hodge operators in spacetime.
	The exact solution is chosen to be exp(x/4)sin(x+y+t+1)cos(y)(2-0.003t-0.00002t^2).
	This implies boundary conditions and the source term.
	Instead of solving the full system, we form a system for a single time step (of length T) and update the solution one time step at a time (for T_REPEATS steps).
	*/
	void waveExampleWhitney3D_TimeStepping(uint order, const double T, const uint T_REPEATS) {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector3)> testFunction{ [&](Vector3 p) -> double { return std::exp(0.25 * p.x) * std::sin(p.x + p.y + p.z + 1) * std::cos(p.y) * (2 - 0.003 * p.z - 0.00002 * p.z * p.z); } };
		std::function<double(Vector3)> sourceTerm{ [&](Vector3 p) -> double { return -0.00001 * std::exp(0.25 * p.x) * (((-4 * p.z - 600) * p.z + 400000) * std::sin(p.y) * std::cos(p.z + p.x + p.y + 1)
			+ std::cos(p.y) * (((-3.875 * p.z - 581.25) * p.z + 387496) * std::sin(p.z + p.x + p.y + 1) + (p.z * (p.z + 142) - 100600) * std::cos(p.z + p.x + p.y + 1))); } };
		std::function<Vector3(Vector3)> testFunctionGradient{ [&](Vector3 p) -> Vector3 { return std::exp(0.25 * p.x) * Vector3(0.25 * (2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.y) * (4 * std::cos(p.x + p.y + p.z + 1) + std::sin(p.x + p.y + p.z + 1)),
			(2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + 2 * p.y + p.z + 1), std::cos(p.y) * ((2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + p.y + p.z + 1) - (0.003 + 0.00004 * p.z) * std::sin(p.x + p.y + p.z + 1))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleWhitney3D_TimeStepping_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		uint maxMeshNumber = 7;
		if (order == 2)
			maxMeshNumber = 6;
		else if (order == 3 || order == 4)
			maxMeshNumber = 5;
		else if (order >= 5)
			maxMeshNumber = 4;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) { //solve the problem on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;

			//solve the problem using Whitney forms of chosen order
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//save maximum error and maximum stepwise average
			double maxErrorFull = 0.0;
			double maxAverage = 0.0;
			double averageSum = 0.0;
			uint averageNum = 0;

			//Create the mesh covering one time step only and its refinement and mark cells with appropriate flags.
			BuilderMesh mesh_old(3);
			createSpacetimeMesh(mesh_old, Lx, Ly, dtime, refs, refs, 1);
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				Vector3 nodePos = mesh_old.getNodePosition3(i);
				if (nodePos.z > dtime - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
					mesh_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
				Vector3 edgePos = mesh_old.getEdgeAverage3(i);
				if (edgePos.z > dtime - eps && edgePos.x > eps && edgePos.x < Lx - eps && edgePos.y > eps && edgePos.y < Ly - eps)
					mesh_old.setEdgeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
				Vector3 facePos = mesh_old.getFaceAverage3(i);
				if (facePos.z > dtime - eps && facePos.x > eps && facePos.x < Lx - eps && facePos.y > eps && facePos.y < Ly - eps)
					mesh_old.setFaceFlag(i, 0);
			}
			BuilderMesh mesh(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh_old, mesh);
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1) {
					Vector3 nodePos = mesh.getNodePosition3(i);
					if (nodePos.z < eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps) {
						mesh.setNodeFlag(i, 2);
					}
					else if (nodePos.z > dtime - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
						mesh.setNodeFlag(i, 3);
				}
			}

			//form the system matrix
			Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
			{
				Buffer<std::unordered_map<uint, double>> star;
				Buffer<uint> refElementIndices(6);
				Buffer<Buffer<uint>> refElements(6);
				for (uint i = 0; i < 6; ++i) {
					refElementIndices[i] = i;
					refElements[i].resize(mesh_old.getBodySize() / 6);
				}
				for (uint i = 0; i < mesh_old.getBodySize(); ++i)
					refElements[i % 6][i / 6] = i;
				ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "CubeTetrahedron", false, true, false);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//restrict to active nodes
			uint boundary2Nodes = (order * refs - 1) * (order * refs - 1);
			uint activeNodes = order * boundary2Nodes;
			Buffer<Buffer<double>> system_matrix_interior(activeNodes);
			for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
					continue;
				std::unordered_map<uint, double>& row = system_matrix[i];
				system_matrix_interior[ii] = Buffer<double>(activeNodes, 0.0);
				for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
					if (mesh.getNodeFlag(j) == 1 || mesh.getNodeFlag(j) == 2)
						continue;
					auto it = row.find(j);
					if (it != row.end() && std::abs(it->second) > 1.0e-14) {
						system_matrix_interior[ii][jj] = it->second;
					}
					++jj;
				}
				++ii;
			}
			Buffer<uint> p;
			decomposeLUP(system_matrix_interior, p, 1e-14);

			//initial condition for update step
			Buffer<double> rhs_prepared(boundary2Nodes);
			Buffer<double> initial_values(boundary2Nodes);
			for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 2) {
					initial_values[j] = testFunction(mesh.getNodePosition3(i));
					rhs_prepared[j] = -integrateNodeDualBoundary(mesh, i, testFunctionGradient, false, true);
					++j;
				}
			}

			uint stepNumber = 1;
			while (true) {
				//form the right-hand side
				Buffer<double> rhs(activeNodes);
				Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						boundaryValues[i] = testFunction(mesh.getNodePosition3(i));
					else if (mesh.getNodeFlag(i) == 2) {
						boundaryValues[i] = initial_values[j++];
					}
				}
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
						continue;
					rhs[j] = integrateDual3Cell(mesh, i, sourceTerm);
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						rhs[j] -= it->second * boundaryValues[it->first];
					}
					if (mesh.getNodeFlag(i) == 2)
						rhs[j] += rhs_prepared[k++];
					++j;
				}

				//get the solution on active nodes by solving the system and on other nodes from the boundary condition
				Buffer<double> sol_interior;
				solveLUP(system_matrix_interior, p, rhs, sol_interior);
				Buffer<double> sol;
				sol.resize(mesh.getNodeSize());
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						sol[i] = testFunction(mesh.getNodePosition3(i));
					else if (mesh.getNodeFlag(i) == 2)
						sol[i] = initial_values[k++];
					else
						sol[i] = sol_interior[j++];
				}


				//stepwise error analysis
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(sol[i] - testFunction(mesh.getNodePosition3(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageSum += averageError;
				averageNum += interiorNodes;
				averageError /= interiorNodes;
				if (maxError > maxErrorFull)
					maxErrorFull = maxError;
				if (averageError > maxAverage)
					maxAverage = averageError;
				if (stepNumber % 10 == 0) { //choose how often to print the result
					std::cout << "Time elapsed: " << stepNumber * dtime << '\n';
					output << "Time elapsed: " << stepNumber * dtime << '\n';
					std::cout << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					output << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
				}

				//stop if done
				if (++stepNumber > T_REPEATS * refs)
					break;

				//prepare initial condition for the next step
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 3) {
						rhs_prepared[j] = integrateDual3Cell(mesh, i, sourceTerm);
						std::unordered_map<uint, double>& row = system_matrix[i];
						for (auto it = row.begin(); it != row.end(); ++it) {
							rhs_prepared[j] -= it->second * sol[it->first];
						}
						initial_values[j] = sol[i];
						++j;
					}
				}

				//move the mesh
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					mesh_old.setNodePosition(i, mesh_old.getNodePosition(i) + Vector4(0.0, 0.0, dtime, 0.0));
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(0.0, 0.0, dtime, 0.0));
				}
			}

			//error analysis for the full solution
			std::cout << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";
			output << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";

			//maximum and average error
			averageSum /= averageNum;
			std::cout << "Max error " << maxErrorFull << " (average " << averageSum << ", maximum stepwise average " << maxAverage << ')' << '\n';
			output << "Max error " << maxErrorFull << " (average " << averageSum << ", maximum stepwise average " << maxAverage << ')' << '\n';
			std::cout << '\n';
			output << '\n';
			output.save(outputFile.str());
		}
	}

	//same as waveExampleWhitney3D_TimeStepping but saves the full solution for the full mesh
	void waveExampleWhitney3D_TimeSteppingFullSolution(uint order, const double T, const uint T_REPEATS) {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector3)> testFunction{ [&](Vector3 p) -> double { return std::exp(0.25 * p.x) * std::sin(p.x + p.y + p.z + 1) * std::cos(p.y) * (2 - 0.003 * p.z - 0.00002 * p.z * p.z); } };
		std::function<double(Vector3)> sourceTerm{ [&](Vector3 p) -> double { return -0.00001 * std::exp(0.25 * p.x) * (((-4 * p.z - 600) * p.z + 400000) * std::sin(p.y) * std::cos(p.z + p.x + p.y + 1)
			+ std::cos(p.y) * (((-3.875 * p.z - 581.25) * p.z + 387496) * std::sin(p.z + p.x + p.y + 1) + (p.z * (p.z + 142) - 100600) * std::cos(p.z + p.x + p.y + 1))); } };
		std::function<Vector3(Vector3)> testFunctionGradient{ [&](Vector3 p) -> Vector3 { return std::exp(0.25 * p.x) * Vector3(0.25 * (2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.y) * (4 * std::cos(p.x + p.y + p.z + 1) + std::sin(p.x + p.y + p.z + 1)),
			(2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + 2 * p.y + p.z + 1), std::cos(p.y) * ((2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + p.y + p.z + 1) - (0.003 + 0.00004 * p.z) * std::sin(p.x + p.y + p.z + 1))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleWhitney3D_TimeSteppingFullSolution_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		uint maxMeshNumber = 7;
		if (order == 2)
			maxMeshNumber = 6;
		else if (order == 3 || order == 4)
			maxMeshNumber = 5;
		else if (order >= 5)
			maxMeshNumber = 4;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) { //solve the problem on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//Create the full mesh in [0, Lx] x [0, Ly] x [0, T * T_REPEATS] and mark cells with appropriate flags.
			uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;
			BuilderMesh meshFull_old(3);
			createSpacetimeMesh(meshFull_old, Lx, Ly, T * T_REPEATS, refs, refs, refs * T_REPEATS);
			meshFull_old.fillBoundaryFlags(1);
			for (uint i = 0; i < meshFull_old.getNodeSize(); ++i) {
				Vector3 nodePos = meshFull_old.getNodePosition3(i);
				if (nodePos.z > T_REPEATS * T - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
					meshFull_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < meshFull_old.getEdgeSize(); ++i) {
				Vector3 edgePos = meshFull_old.getEdgeAverage3(i);
				if (edgePos.z > T_REPEATS * T - eps && edgePos.x > eps && edgePos.x < Lx - eps && edgePos.y > eps && edgePos.y < Ly - eps)
					meshFull_old.setEdgeFlag(i, 0);
			}
			for (uint i = 0; i < meshFull_old.getFaceSize(); ++i) {
				Vector3 facePos = meshFull_old.getFaceAverage3(i);
				if (facePos.z > T_REPEATS * T - eps && facePos.x > eps && facePos.x < Lx - eps && facePos.y > eps && facePos.y < Ly - eps)
					meshFull_old.setFaceFlag(i, 0);
			}

			//solve the problem using Whitney forms of chosen order
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//Refine the full mesh into small simplices and mark cells with appropriate flags.
			BuilderMesh meshFull(3);
			SmallSimplexPartition3D sspFull(order);
			sspFull.refineMesh(meshFull_old, meshFull);
			meshFull.fillBoundaryFlags(1);
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1) {
					Vector3 nodePos = meshFull.getNodePosition3(i);
					if (nodePos.z < eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps) {
						meshFull.setNodeFlag(i, 2);
					}
					else if (nodePos.z > T_REPEATS * T - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
						meshFull.setNodeFlag(i, 3);
				}
			}

			//prepare memory and variables to store values in the full solution on each time step and set boundary values using boundary condition
			uint nodesUpdated = 0;
			uint edgesUpdated = 0;
			uint facesUpdated = 0;
			uint bodiesUpdated = 0;
			const uint smallNodesInFaces = binCoef(order - 1, 2);
			const uint smallNodesInBodies = binCoef(order - 1, 3);
			Buffer<double> solFull(meshFull.getNodeSize());
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1 || meshFull.getNodeFlag(i) == 2) {
					solFull[i] = testFunction(meshFull.getNodePosition3(i));
				}
			}

			//Create the mesh covering one time step only and its refinement and mark cells with appropriate flags.
			BuilderMesh mesh_old(3);
			createSpacetimeMesh(mesh_old, Lx, Ly, dtime, refs, refs, 1);
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				Vector3 nodePos = mesh_old.getNodePosition3(i);
				if (nodePos.z > dtime - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
					mesh_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
				Vector3 edgePos = mesh_old.getEdgeAverage3(i);
				if (edgePos.z > dtime - eps && edgePos.x > eps && edgePos.x < Lx - eps && edgePos.y > eps && edgePos.y < Ly - eps)
					mesh_old.setEdgeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
				Vector3 facePos = mesh_old.getFaceAverage3(i);
				if (facePos.z > dtime - eps && facePos.x > eps && facePos.x < Lx - eps && facePos.y > eps && facePos.y < Ly - eps)
					mesh_old.setFaceFlag(i, 0);
			}
			BuilderMesh mesh(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh_old, mesh);
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1) {
					Vector3 nodePos = mesh.getNodePosition3(i);
					if (nodePos.z < eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps) {
						mesh.setNodeFlag(i, 2);
					}
					else if (nodePos.z > dtime - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
						mesh.setNodeFlag(i, 3);
				}
			}

			//form the system matrix
			Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
			{
				Buffer<std::unordered_map<uint, double>> star;
				Buffer<uint> refElementIndices(6);
				Buffer<Buffer<uint>> refElements(6);
				for (uint i = 0; i < 6; ++i) {
					refElementIndices[i] = i;
					refElements[i].resize(mesh_old.getBodySize() / 6);
				}
				for (uint i = 0; i < mesh_old.getBodySize(); ++i)
					refElements[i % 6][i / 6] = i;
				ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "CubeTetrahedron", false, true, false);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//restrict to active nodes
			uint boundary2Nodes = (order * refs - 1) * (order * refs - 1);
			uint activeNodes = order * boundary2Nodes;
			Buffer<Buffer<double>> system_matrix_interior(activeNodes);
			for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
					continue;
				std::unordered_map<uint, double>& row = system_matrix[i];
				system_matrix_interior[ii] = Buffer<double>(activeNodes, 0.0);
				for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
					if (mesh.getNodeFlag(j) == 1 || mesh.getNodeFlag(j) == 2)
						continue;
					auto it = row.find(j);
					if (it != row.end() && std::abs(it->second) > 1.0e-14) {
						system_matrix_interior[ii][jj] = it->second;
					}
					++jj;
				}
				++ii;
			}
			Buffer<uint> p;
			decomposeLUP(system_matrix_interior, p, 1e-14);

			//initial condition for update step
			Buffer<double> rhs_prepared(boundary2Nodes);
			Buffer<double> initial_values(boundary2Nodes);
			for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 2) {
					initial_values[j] = testFunction(mesh.getNodePosition3(i));
					rhs_prepared[j] = -integrateNodeDualBoundary(mesh, i, testFunctionGradient, false, true);
					++j;
				}
			}

			uint stepNumber = 1;
			while (true) {
				//form the right-hand side
				Buffer<double> rhs(activeNodes);
				Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						boundaryValues[i] = testFunction(mesh.getNodePosition3(i));
					else if (mesh.getNodeFlag(i) == 2) {
						boundaryValues[i] = initial_values[j++];
					}
				}
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
						continue;
					rhs[j] = integrateDual3Cell(mesh, i, sourceTerm);
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						rhs[j] -= it->second * boundaryValues[it->first];
					}
					if (mesh.getNodeFlag(i) == 2)
						rhs[j] += rhs_prepared[k++];
					++j;
				}

				//get the solution on active nodes by solving the system and on other nodes from the boundary condition
				Buffer<double> sol_interior;
				solveLUP(system_matrix_interior, p, rhs, sol_interior);
				Buffer<double> sol;
				sol.resize(mesh.getNodeSize());
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						sol[i] = testFunction(mesh.getNodePosition3(i));
					else if (mesh.getNodeFlag(i) == 2)
						sol[i] = initial_values[k++];
					else
						sol[i] = sol_interior[j++];
				}

				//store the values in full solution
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					if (mesh_old.getNodeFlag(i) == 1)
						continue;
					while (meshFull_old.getNodeFlag(nodesUpdated) == 1)
						++nodesUpdated;
					solFull[nodesUpdated++] = sol[i];
				}
				for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
					if (mesh_old.getEdgeFlag(i) == 1)
						continue;
					while (meshFull_old.getEdgeFlag(edgesUpdated) == 1)
						++edgesUpdated;
					for (uint j = 0; j < order - 1; ++j)
						solFull[sspFull.smallNodesEdgeList(edgesUpdated, j)] = sol[ssp.smallNodesEdgeList(i, j)];
					++edgesUpdated;
				}
				for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
					if (mesh_old.getFaceFlag(i) == 1)
						continue;
					while (meshFull_old.getFaceFlag(facesUpdated) == 1)
						++facesUpdated;
					for (uint j = 0; j < smallNodesInFaces; ++j)
						solFull[sspFull.smallNodesFaceList(facesUpdated, j)] = sol[ssp.smallNodesFaceList(i, j)];
					++facesUpdated;
				}
				for (uint i = 0; i < mesh_old.getBodySize(); ++i) {
					for (uint j = 0; j < smallNodesInBodies; ++j)
						solFull[sspFull.smallNodesBodyList(bodiesUpdated, j)] = sol[ssp.smallNodesBodyList(i, j)];
					++bodiesUpdated;
				}

				//stepwise error analysis
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(sol[i] - testFunction(mesh.getNodePosition3(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageError /= interiorNodes;
				if (stepNumber % 10 == 0) { //choose how often to print the result
					std::cout << "Time elapsed: " << stepNumber * dtime << '\n';
					output << "Time elapsed: " << stepNumber * dtime << '\n';
					std::cout << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					output << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					/*double l2error1 = computeDifferenceL2Norm0Form(mesh_old, ssp, sol, testFunction, 1000);
					std::cout << "L2 error: " << l2error1 << '\n';
					output << "L2 error: " << l2error1 << '\n';
					Buffer<double> dec_sol_grad(mesh.getEdgeSize(), 0.0);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						const Buffer<uint>& nodes = mesh.getEdgeNodes(i);
						for (uint j = 0; j < nodes.size(); ++j)
							dec_sol_grad[i] += mesh.getEdgeIncidence(i, nodes[j]) * sol[nodes[j]];
					}
					double graderror1 = computeDifferenceL2Norm1Form(mesh_old, ssp, dec_sol_grad, testFunctionGradient, 1000);
					std::cout << "Gradient L2 error: " << graderror1 << '\n';
					output << "Gradient L2 error: " << graderror1 << '\n';
					double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
					std::cout << "H1 error: " << h1error1 << '\n';
					output << "H1 error: " << h1error1 << '\n';*/
				}

				//stop if done
				if (++stepNumber > T_REPEATS * refs)
					break;

				//prepare initial condition for the next step
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 3) {
						rhs_prepared[j] = integrateDual3Cell(mesh, i, sourceTerm);
						std::unordered_map<uint, double>& row = system_matrix[i];
						for (auto it = row.begin(); it != row.end(); ++it) {
							rhs_prepared[j] -= it->second * sol[it->first];
						}
						initial_values[j] = sol[i];
						++j;
					}
				}

				//move the mesh
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					mesh_old.setNodePosition(i, mesh_old.getNodePosition(i) + Vector4(0.0, 0.0, dtime, 0.0));
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(0.0, 0.0, dtime, 0.0));
				}
			}

			//error analysis for the full solution
			std::cout << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";
			output << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";

			//maximum and average error
			double maxError = 0.0;
			double averageError = 0.0;
			uint interiorNodes = 0;
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1 || meshFull.getNodeFlag(i) == 2)
					continue;
				++interiorNodes;
				double diff = std::abs(solFull[i] - testFunction(meshFull.getNodePosition3(i)));
				if (maxError < diff)
					maxError = diff;
				averageError += diff;
			}
			averageError /= interiorNodes;
			std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			output << "Max error " << maxError << " (average " << averageError << ')' << '\n';

			//L2 error and H1 error
			/*double l2error1 = computeDifferenceL2Norm0Form(meshFull_old, sspFull, solFull, testFunction, 1000);
			std::cout << "L2 error: " << l2error1 << '\n';
			output << "L2 error: " << l2error1 << '\n';
			Buffer<double> dec_solFull_grad(meshFull.getEdgeSize(), 0.0);
			for (uint i = 0; i < meshFull.getEdgeSize(); ++i) {
				const Buffer<uint>& nodes = meshFull.getEdgeNodes(i);
				for (uint j = 0; j < nodes.size(); ++j)
					dec_solFull_grad[i] += meshFull.getEdgeIncidence(i, nodes[j]) * solFull[nodes[j]];
			}
			double graderror1 = computeDifferenceL2Norm1Form(meshFull_old, sspFull, dec_solFull_grad, testFunctionGradient, 1000);
			std::cout << "Gradient L2 error: " << graderror1 << '\n';
			output << "Gradient L2 error: " << graderror1 << '\n';
			double h1error1 = std::sqrt(l2error1 * l2error1 + graderror1 * graderror1);
			std::cout << "H1 error: " << h1error1 << '\n';
			output << "H1 error: " << h1error1 << '\n';*/

			std::cout << '\n';
			output << '\n';
			output.save(outputFile.str());
		}
	}

	//same as waveExampleWhitney2D_TimeStepping but with cubical forms
	void waveExampleCubical2D_TimeStepping(uint order, const double T, const uint T_REPEATS) {
		const double L = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector2)> testFunction{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * std::sin(p.x + p.y + 2) * std::cos(p.y); } };
		std::function<double(Vector2)> sourceTerm{ [&](Vector2 p) -> double { return std::exp(p.x / 4) * (0.5 * std::cos(p.y) * std::cos(2 + p.x + p.y) + 2 * std::cos(2 + p.x + p.y) * std::sin(p.y)
			+ 17.0 / 16 * std::cos(p.y) * std::sin(2 + p.x + p.y)); } };
		std::function<Vector2(Vector2)> testFunctionGradient{ [&](Vector2 p) -> Vector2 { return Vector2(std::exp(p.x / 4) * std::cos(p.y) * (std::sin(p.x + p.y + 2) / 4 + std::cos(p.x + p.y + 2)),
			std::exp(p.x / 4) * (std::cos(p.y) * std::cos(p.x + p.y + 2) - std::sin(p.y) * std::sin(p.x + p.y + 2))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleCubical2D_TimeStepping_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		uint maxMeshNumber = 8;
		if (order == 2)
			maxMeshNumber = 7;
		else if (order == 3 || order == 4)
			maxMeshNumber = 6;
		else if (order >= 5)
			maxMeshNumber = 5;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) { //solve the problem on meshes of different grain
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;

			//solve the problem using Whitney forms of chosen order
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//save maximum error and maximum stepwise average
			double maxErrorFull = 0.0;
			double maxAverage = 0.0;
			double averageSum = 0.0;
			uint averageNum = 0;

			//Create the mesh covering one time step only and its refinement and mark cells with appropriate flags.
			BuilderMesh mesh_old(2);
			createCartesianMesh(mesh_old, L, dtime, refs, 1);
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				Vector2 nodePos = mesh_old.getNodePosition2(i);
				if (nodePos.y > dtime - eps && nodePos.x > eps && nodePos.x < L - eps)
					mesh_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
				Vector2 edgePos = mesh_old.getEdgeAverage2(i);
				if (edgePos.y > dtime - eps && edgePos.x > eps && edgePos.x < L - eps)
					mesh_old.setEdgeFlag(i, 0);
			}
			BuilderMesh mesh(2);
			SmallCubePartition2D scp(order);
			scp.refineMesh(mesh_old, mesh);
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1) {
					Vector2 nodePos = mesh.getNodePosition2(i);
					if (nodePos.y < eps && nodePos.x > eps && nodePos.x < L - eps) {
						mesh.setNodeFlag(i, 2);
					}
					else if (nodePos.y > dtime - eps && nodePos.x > eps && nodePos.x < L - eps)
						mesh.setNodeFlag(i, 3);
				}
			}

			//form the system matrix
			Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
			{
				Buffer<std::unordered_map<uint, double>> star;
				scp.formHodgeMatrix1Forms(star, true);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//restrict to active nodes
			uint boundary2Nodes = order * refs - 1;
			uint activeNodes = order * boundary2Nodes;
			Buffer<Buffer<double>> system_matrix_interior(activeNodes);
			for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
					continue;
				std::unordered_map<uint, double>& row = system_matrix[i];
				system_matrix_interior[ii] = Buffer<double>(activeNodes, 0.0);
				for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
					if (mesh.getNodeFlag(j) == 1 || mesh.getNodeFlag(j) == 2)
						continue;
					auto it = row.find(j);
					if (it != row.end() && std::abs(it->second) > 1.0e-14) {
						system_matrix_interior[ii][jj] = it->second;
					}
					++jj;
				}
				++ii;
			}
			Buffer<uint> p;
			decomposeLUP(system_matrix_interior, p, 1e-14);

			//initial condition for update step
			Buffer<double> rhs_prepared(boundary2Nodes);
			Buffer<double> initial_values(boundary2Nodes);
			for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 2) {
					initial_values[j] = testFunction(mesh.getNodePosition2(i));
					rhs_prepared[j] = -integrateNodeDualBoundary(mesh, i, testFunctionGradient, false, true);
					++j;
				}
			}

			uint stepNumber = 1;
			while (true) {
				//form the right-hand side
				Buffer<double> rhs(activeNodes);
				Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						boundaryValues[i] = testFunction(mesh.getNodePosition2(i));
					else if (mesh.getNodeFlag(i) == 2) {
						boundaryValues[i] = initial_values[j++];
					}
				}
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
						continue;
					rhs[j] = integrateDual2Cell(mesh, i, sourceTerm);
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						rhs[j] -= it->second * boundaryValues[it->first];
					}
					if (mesh.getNodeFlag(i) == 2)
						rhs[j] += rhs_prepared[k++];
					++j;
				}

				//get the solution on active nodes by solving the system and on other nodes from the boundary condition
				Buffer<double> sol_interior;
				solveLUP(system_matrix_interior, p, rhs, sol_interior);
				Buffer<double> sol;
				sol.resize(mesh.getNodeSize());
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						sol[i] = testFunction(mesh.getNodePosition2(i));
					else if (mesh.getNodeFlag(i) == 2)
						sol[i] = initial_values[k++];
					else
						sol[i] = sol_interior[j++];
				}

				//stepwise error analysis
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(sol[i] - testFunction(mesh.getNodePosition2(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageSum += averageError;
				averageNum += interiorNodes;
				averageError /= interiorNodes;
				if (stepNumber % 10 == 0) { //choose how often to print the result
					std::cout << "Time elapsed: " << stepNumber * dtime << '\n';
					output << "Time elapsed: " << stepNumber * dtime << '\n';
					std::cout << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					output << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
				}
				if (maxError > maxErrorFull)
					maxErrorFull = maxError;
				if (averageError > maxAverage)
					maxAverage = averageError;

				//stop if done
				if (++stepNumber > T_REPEATS * refs)
					break;

				//prepare initial condition for the next step
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 3) {
						rhs_prepared[j] = integrateDual2Cell(mesh, i, sourceTerm);
						std::unordered_map<uint, double>& row = system_matrix[i];
						for (auto it = row.begin(); it != row.end(); ++it) {
							rhs_prepared[j] -= it->second * sol[it->first];
						}
						initial_values[j] = sol[i];
						++j;
					}
				}

				//move the mesh
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					mesh_old.setNodePosition(i, mesh_old.getNodePosition(i) + Vector4(0.0, dtime, 0.0, 0.0));
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(0.0, dtime, 0.0, 0.0));
				}
			}

			//error analysis for the full solution
			std::cout << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";
			output << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";

			//maximum and average error
			averageSum /= averageNum;
			std::cout << "Max error " << maxErrorFull << " (average " << averageSum << ", maximum stepwise average " << maxAverage << ')' << '\n';
			output << "Max error " << maxErrorFull << " (average " << averageSum << ", maximum stepwise average " << maxAverage << ')' << '\n';
			std::cout << '\n';
			output << '\n';
			output.save(outputFile.str());
		}
	}

	//same as waveExampleWhitney3D_TimeStepping but with cubical forms
	void waveExampleCubical3D_TimeStepping(uint order, const double T, const uint T_REPEATS) {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector3)> testFunction{ [&](Vector3 p) -> double { return std::exp(0.25 * p.x) * std::sin(p.x + p.y + p.z + 1) * std::cos(p.y) * (2 - 0.003 * p.z - 0.00002 * p.z * p.z); } };
		std::function<double(Vector3)> sourceTerm{ [&](Vector3 p) -> double { return -0.00001 * std::exp(0.25 * p.x) * (((-4 * p.z - 600) * p.z + 400000) * std::sin(p.y) * std::cos(p.z + p.x + p.y + 1)
			+ std::cos(p.y) * (((-3.875 * p.z - 581.25) * p.z + 387496) * std::sin(p.z + p.x + p.y + 1) + (p.z * (p.z + 142) - 100600) * std::cos(p.z + p.x + p.y + 1))); } };
		std::function<Vector3(Vector3)> testFunctionGradient{ [&](Vector3 p) -> Vector3 { return std::exp(0.25 * p.x) * Vector3(0.25 * (2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.y) * (4 * std::cos(p.x + p.y + p.z + 1) + std::sin(p.x + p.y + p.z + 1)),
			(2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + 2 * p.y + p.z + 1), std::cos(p.y) * ((2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + p.y + p.z + 1) - (0.003 + 0.00004 * p.z) * std::sin(p.x + p.y + p.z + 1))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleCubical3D_TimeStepping_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		uint maxMeshNumber = 7;
		if (order == 2)
			maxMeshNumber = 6;
		else if (order == 3 || order == 4)
			maxMeshNumber = 5;
		else if (order >= 5)
			maxMeshNumber = 4;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) { //solve the problem on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;

			//solve the problem using cubical forms of chosen order
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//save maximum error and maximum stepwise average
			double maxErrorFull = 0.0;
			double maxAverage = 0.0;
			double averageSum = 0.0;
			uint averageNum = 0;

			//Create the mesh covering one time step only and its refinement and mark cells with appropriate flags.
			BuilderMesh mesh_old(3);
			createCartesianMesh(mesh_old, Lx, Ly, dtime, refs, refs, 1);
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				Vector3 nodePos = mesh_old.getNodePosition3(i);
				if (nodePos.z > dtime - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
					mesh_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
				Vector3 edgePos = mesh_old.getEdgeAverage3(i);
				if (edgePos.z > dtime - eps && edgePos.x > eps && edgePos.x < Lx - eps && edgePos.y > eps && edgePos.y < Ly - eps)
					mesh_old.setEdgeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
				Vector3 facePos = mesh_old.getFaceAverage3(i);
				if (facePos.z > dtime - eps && facePos.x > eps && facePos.x < Lx - eps && facePos.y > eps && facePos.y < Ly - eps)
					mesh_old.setFaceFlag(i, 0);
			}
			BuilderMesh mesh(3);
			SmallCubePartition3D scp(order);
			scp.refineMesh(mesh_old, mesh);
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1) {
					Vector3 nodePos = mesh.getNodePosition3(i);
					if (nodePos.z < eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps) {
						mesh.setNodeFlag(i, 2);
					}
					else if (nodePos.z > dtime - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
						mesh.setNodeFlag(i, 3);
				}
			}

			//form the system matrix
			Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
			{
				Buffer<std::unordered_map<uint, double>> star;
				scp.formHodgeMatrix1Forms(star, true);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//restrict to active nodes
			uint boundary2Nodes = (order * refs - 1) * (order * refs - 1);
			uint activeNodes = order * boundary2Nodes;
			Buffer<Buffer<double>> system_matrix_interior(activeNodes);
			for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
					continue;
				std::unordered_map<uint, double>& row = system_matrix[i];
				system_matrix_interior[ii] = Buffer<double>(activeNodes, 0.0);
				for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
					if (mesh.getNodeFlag(j) == 1 || mesh.getNodeFlag(j) == 2)
						continue;
					auto it = row.find(j);
					if (it != row.end() && it->second != 0.0) {
						system_matrix_interior[ii][jj] = it->second;
					}
					++jj;
				}
				++ii;
			}
			Buffer<uint> p;
			decomposeLUP(system_matrix_interior, p, 1e-14);

			//initial condition for update step
			Buffer<double> rhs_prepared(boundary2Nodes);
			Buffer<double> initial_values(boundary2Nodes);
			for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 2) {
					initial_values[j] = testFunction(mesh.getNodePosition3(i));
					rhs_prepared[j] = -integrateNodeDualBoundaryCartesianMesh(mesh, i, testFunctionGradient, true);
					++j;
				}
			}

			uint stepNumber = 1;
			while (true) {
				//form the right-hand side
				Buffer<double> rhs(activeNodes);
				Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						boundaryValues[i] = testFunction(mesh.getNodePosition3(i));
					else if (mesh.getNodeFlag(i) == 2) {
						boundaryValues[i] = initial_values[j++];
					}
				}
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
						continue;
					rhs[j] = integrateDual3CellCartesianMesh(mesh, i, sourceTerm);
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						rhs[j] -= it->second * boundaryValues[it->first];
					}
					if (mesh.getNodeFlag(i) == 2)
						rhs[j] += rhs_prepared[k++];
					++j;
				}

				//get the solution on active nodes by solving the system and on other nodes from the boundary condition
				Buffer<double> sol_interior;
				solveLUP(system_matrix_interior, p, rhs, sol_interior);
				Buffer<double> sol;
				sol.resize(mesh.getNodeSize());
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						sol[i] = testFunction(mesh.getNodePosition3(i));
					else if (mesh.getNodeFlag(i) == 2)
						sol[i] = initial_values[k++];
					else
						sol[i] = sol_interior[j++];
				}

				//stepwise error analysis
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(sol[i] - testFunction(mesh.getNodePosition3(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageSum += averageError;
				averageNum += interiorNodes;
				averageError /= interiorNodes;
				if (stepNumber % 10 == 0) { //choose how often to print the result
					std::cout << "Time elapsed: " << stepNumber * dtime << '\n';
					output << "Time elapsed: " << stepNumber * dtime << '\n';
					std::cout << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					output << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
				}
				if (maxError > maxErrorFull)
					maxErrorFull = maxError;
				if (averageError > maxAverage)
					maxAverage = averageError;

				//stop if done
				if (++stepNumber > T_REPEATS * refs)
					break;

				//prepare initial condition for the next step
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 3) {
						rhs_prepared[j] = integrateDual3CellCartesianMesh(mesh, i, sourceTerm);
						std::unordered_map<uint, double>& row = system_matrix[i];
						for (auto it = row.begin(); it != row.end(); ++it) {
							rhs_prepared[j] -= it->second * sol[it->first];
						}
						initial_values[j] = sol[i];
						++j;
					}
				}

				//move the mesh
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					mesh_old.setNodePosition(i, mesh_old.getNodePosition(i) + Vector4(0.0, 0.0, dtime, 0.0));
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(0.0, 0.0, dtime, 0.0));
				}
			}

			//error analysis for the full solution
			std::cout << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";
			output << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";

			//maximum and average error
			averageSum /= averageNum;
			std::cout << "Max error " << maxErrorFull << " (average " << averageSum << ", maximum stepwise average " << maxAverage << ')' << '\n';
			output << "Max error " << maxErrorFull << " (average " << averageSum << ", maximum stepwise average " << maxAverage << ')' << '\n';
			std::cout << '\n';
			output << '\n';
			output.save(outputFile.str());
		}
	}

	//same as waveExampleCubical3D_TimeStepping but saves the full solution for the full mesh
	void waveExampleCubical3D_TimeSteppingFullSolution(uint order, const double T, const uint T_REPEATS) {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector3)> testFunction{ [&](Vector3 p) -> double { return std::exp(0.25 * p.x) * std::sin(p.x + p.y + p.z + 1) * std::cos(p.y) * (2 - 0.003 * p.z - 0.00002 * p.z * p.z); } };
		std::function<double(Vector3)> sourceTerm{ [&](Vector3 p) -> double { return -0.00001 * std::exp(0.25 * p.x) * (((-4 * p.z - 600) * p.z + 400000) * std::sin(p.y) * std::cos(p.z + p.x + p.y + 1)
			+ std::cos(p.y) * (((-3.875 * p.z - 581.25) * p.z + 387496) * std::sin(p.z + p.x + p.y + 1) + (p.z * (p.z + 142) - 100600) * std::cos(p.z + p.x + p.y + 1))); } };
		std::function<Vector3(Vector3)> testFunctionGradient{ [&](Vector3 p) -> Vector3 { return std::exp(0.25 * p.x) * Vector3(0.25 * (2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.y) * (4 * std::cos(p.x + p.y + p.z + 1) + std::sin(p.x + p.y + p.z + 1)),
			(2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + 2 * p.y + p.z + 1), std::cos(p.y) * ((2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + p.y + p.z + 1) - (0.003 + 0.00004 * p.z) * std::sin(p.x + p.y + p.z + 1))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleCubial3D_TimeSteppingFullSolution_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		uint maxMeshNumber = 7;
		if (order == 2)
			maxMeshNumber = 6;
		else if (order == 3 || order == 4)
			maxMeshNumber = 5;
		else if (order >= 5)
			maxMeshNumber = 4;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) { //solve the problem on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//Create the full mesh in [0, Lx] x [0, Ly] x [0, T * T_REPEATS] and mark cells with appropriate flags.
			uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;
			BuilderMesh meshFull_old(3);
			createCartesianMesh(meshFull_old, Lx, Ly, T * T_REPEATS, refs, refs, refs * T_REPEATS);
			meshFull_old.fillBoundaryFlags(1);
			for (uint i = 0; i < meshFull_old.getNodeSize(); ++i) {
				Vector3 nodePos = meshFull_old.getNodePosition3(i);
				if (nodePos.z > T_REPEATS * T - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
					meshFull_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < meshFull_old.getEdgeSize(); ++i) {
				Vector3 edgePos = meshFull_old.getEdgeAverage3(i);
				if (edgePos.z > T_REPEATS * T - eps && edgePos.x > eps && edgePos.x < Lx - eps && edgePos.y > eps && edgePos.y < Ly - eps)
					meshFull_old.setEdgeFlag(i, 0);
			}
			for (uint i = 0; i < meshFull_old.getFaceSize(); ++i) {
				Vector3 facePos = meshFull_old.getFaceAverage3(i);
				if (facePos.z > T_REPEATS * T - eps && facePos.x > eps && facePos.x < Lx - eps && facePos.y > eps && facePos.y < Ly - eps)
					meshFull_old.setFaceFlag(i, 0);
			}

			//solve the problem using cubical forms of chosen order
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//Refine the full mesh into small cubes and mark cells with appropriate flags.
			BuilderMesh meshFull(3);
			SmallCubePartition3D scpFull(order);
			scpFull.refineMesh(meshFull_old, meshFull);
			meshFull.fillBoundaryFlags(1);
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1) {
					Vector3 nodePos = meshFull.getNodePosition3(i);
					if (nodePos.z < eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps) {
						meshFull.setNodeFlag(i, 2);
					}
					else if (nodePos.z > T_REPEATS * T - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
						meshFull.setNodeFlag(i, 3);
				}
			}

			//prepare memory and variables to store values in the full solution on each time step and set boundary values using boundary condition
			uint nodesUpdated = 0;
			uint edgesUpdated = 0;
			uint facesUpdated = 0;
			uint bodiesUpdated = 0;
			const uint smallNodesInFaces = (order - 1) * (order - 1);
			const uint smallNodesInBodies = (order - 1) * smallNodesInFaces;
			Buffer<double> solFull(meshFull.getNodeSize());
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1 || meshFull.getNodeFlag(i) == 2) {
					solFull[i] = testFunction(meshFull.getNodePosition3(i));
				}
			}

			//Create the mesh covering one time step only and its refinement and mark cells with appropriate flags.
			BuilderMesh mesh_old(3);
			createCartesianMesh(mesh_old, Lx, Ly, dtime, refs, refs, 1);
			mesh_old.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
				Vector3 nodePos = mesh_old.getNodePosition3(i);
				if (nodePos.z > dtime - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
					mesh_old.setNodeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
				Vector3 edgePos = mesh_old.getEdgeAverage3(i);
				if (edgePos.z > dtime - eps && edgePos.x > eps && edgePos.x < Lx - eps && edgePos.y > eps && edgePos.y < Ly - eps)
					mesh_old.setEdgeFlag(i, 0);
			}
			for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
				Vector3 facePos = mesh_old.getFaceAverage3(i);
				if (facePos.z > dtime - eps && facePos.x > eps && facePos.x < Lx - eps && facePos.y > eps && facePos.y < Ly - eps)
					mesh_old.setFaceFlag(i, 0);
			}
			BuilderMesh mesh(3);
			SmallCubePartition3D scp(order);
			scp.refineMesh(mesh_old, mesh);
			mesh.fillBoundaryFlags(1);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1) {
					Vector3 nodePos = mesh.getNodePosition3(i);
					if (nodePos.z < eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps) {
						mesh.setNodeFlag(i, 2);
					}
					else if (nodePos.z > dtime - eps && nodePos.x > eps && nodePos.x < Lx - eps && nodePos.y > eps && nodePos.y < Ly - eps)
						mesh.setNodeFlag(i, 3);
				}
			}

			//form the system matrix
			Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize());
			{
				Buffer<std::unordered_map<uint, double>> star;
				scp.formHodgeMatrix1Forms(star, true);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//restrict to active nodes
			uint boundary2Nodes = (order * refs - 1) * (order * refs - 1);
			uint activeNodes = order * boundary2Nodes;
			Buffer<Buffer<double>> system_matrix_interior(activeNodes);
			for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
					continue;
				std::unordered_map<uint, double>& row = system_matrix[i];
				system_matrix_interior[ii] = Buffer<double>(activeNodes, 0.0);
				for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
					if (mesh.getNodeFlag(j) == 1 || mesh.getNodeFlag(j) == 2)
						continue;
					auto it = row.find(j);
					if (it != row.end() && it->second != 0.0) {
						system_matrix_interior[ii][jj] = it->second;
					}
					++jj;
				}
				++ii;
			}
			Buffer<uint> p;
			decomposeLUP(system_matrix_interior, p, 1e-14);

			//initial condition for update step
			Buffer<double> rhs_prepared(boundary2Nodes);
			Buffer<double> initial_values(boundary2Nodes);
			for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
				if (mesh.getNodeFlag(i) == 2) {
					initial_values[j] = testFunction(mesh.getNodePosition3(i));
					rhs_prepared[j] = -integrateNodeDualBoundaryCartesianMesh(mesh, i, testFunctionGradient, true);
					++j;
				}
			}

			uint stepNumber = 1;
			while (true) {
				//form the right-hand side
				Buffer<double> rhs(activeNodes);
				Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						boundaryValues[i] = testFunction(mesh.getNodePosition3(i));
					else if (mesh.getNodeFlag(i) == 2) {
						boundaryValues[i] = initial_values[j++];
					}
				}
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
						continue;
					rhs[j] = integrateDual3CellCartesianMesh(mesh, i, sourceTerm);
					std::unordered_map<uint, double>& row = system_matrix[i];
					for (auto it = row.begin(); it != row.end(); ++it) {
						rhs[j] -= it->second * boundaryValues[it->first];
					}
					if (mesh.getNodeFlag(i) == 2)
						rhs[j] += rhs_prepared[k++];
					++j;
				}

				//get the solution on active nodes by solving the system and on other nodes from the boundary condition
				Buffer<double> sol_interior;
				solveLUP(system_matrix_interior, p, rhs, sol_interior);
				Buffer<double> sol;
				sol.resize(mesh.getNodeSize());
				for (uint i = 0, j = 0, k = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						sol[i] = testFunction(mesh.getNodePosition3(i));
					else if (mesh.getNodeFlag(i) == 2)
						sol[i] = initial_values[k++];
					else
						sol[i] = sol_interior[j++];
				}

				//store the values in full solution
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					if (mesh_old.getNodeFlag(i) == 1)
						continue;
					while (meshFull_old.getNodeFlag(nodesUpdated) == 1)
						++nodesUpdated;
					solFull[nodesUpdated++] = sol[i];
				}
				for (uint i = 0; i < mesh_old.getEdgeSize(); ++i) {
					if (mesh_old.getEdgeFlag(i) == 1)
						continue;
					while (meshFull_old.getEdgeFlag(edgesUpdated) == 1)
						++edgesUpdated;
					for (uint j = 0; j < order - 1; ++j)
						solFull[scpFull.smallNodesEdgeList(edgesUpdated, j)] = sol[scp.smallNodesEdgeList(i, j)];
					++edgesUpdated;
				}
				for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
					if (mesh_old.getFaceFlag(i) == 1)
						continue;
					while (meshFull_old.getFaceFlag(facesUpdated) == 1)
						++facesUpdated;
					for (uint j = 0; j < smallNodesInFaces; ++j)
						solFull[scpFull.smallNodesFaceList(facesUpdated, j)] = sol[scp.smallNodesFaceList(i, j)];
					++facesUpdated;
				}
				for (uint i = 0; i < mesh_old.getBodySize(); ++i) {
					for (uint j = 0; j < smallNodesInBodies; ++j)
						solFull[scpFull.smallNodesBodyList(bodiesUpdated, j)] = sol[scp.smallNodesBodyList(i, j)];
					++bodiesUpdated;
				}

				//stepwise error analysis
				double maxError = 0.0;
				double averageError = 0.0;
				uint interiorNodes = 0;
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
						continue;
					++interiorNodes;
					double diff = std::abs(sol[i] - testFunction(mesh.getNodePosition3(i)));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageError /= interiorNodes;
				if (stepNumber % 10 == 0) { //choose how often to print the result
					std::cout << "Time elapsed: " << stepNumber * dtime << '\n';
					output << "Time elapsed: " << stepNumber * dtime << '\n';
					std::cout << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					output << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
				}

				//stop if done
				if (++stepNumber > T_REPEATS * refs)
					break;

				//prepare initial condition for the next step
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 3) {
						rhs_prepared[j] = integrateDual3CellCartesianMesh(mesh, i, sourceTerm);
						std::unordered_map<uint, double>& row = system_matrix[i];
						for (auto it = row.begin(); it != row.end(); ++it) {
							rhs_prepared[j] -= it->second * sol[it->first];
						}
						initial_values[j] = sol[i];
						++j;
					}
				}

				//move the mesh
				for (uint i = 0; i < mesh_old.getNodeSize(); ++i) {
					mesh_old.setNodePosition(i, mesh_old.getNodePosition(i) + Vector4(0.0, 0.0, dtime, 0.0));
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(0.0, 0.0, dtime, 0.0));
				}
			}

			//error analysis for the full solution
			std::cout << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";
			output << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";

			//maximum and average error
			double maxError = 0.0;
			double averageError = 0.0;
			uint interiorNodes = 0;
			for (uint i = 0; i < meshFull.getNodeSize(); ++i) {
				if (meshFull.getNodeFlag(i) == 1 || meshFull.getNodeFlag(i) == 2)
					continue;
				++interiorNodes;
				double diff = std::abs(solFull[i] - testFunction(meshFull.getNodePosition3(i)));
				if (maxError < diff)
					maxError = diff;
				averageError += diff;
			}
			averageError /= interiorNodes;
			std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			output << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			std::cout << '\n';
			output << '\n';
			output.save(outputFile.str());
		}
	}

	/*
	Solves inhomogenous wave equation in domain [0, Lx] x [0, Ly] x [0, T * T_REPEATS], using higher order cubical forms and discrete Hodge operators in spacetime.
	Unlike in waveExampleCubical3D_TimeStepping, the values are updated one small step (that corresponds to small cubes) at a time.
	Hence changing the source term in the future does not affect the solution on nodes in the past.
	The exact solution is chosen to be exp(x/4)sin(x+y+t+1)cos(y)(2-0.003t-0.00002t^2).
	This implies boundary conditions and the source term.
	*/
	void waveExampleSmallSteps3D(uint order, const double T, const uint T_REPEATS) {
		if (order == 1) { //method coincides with normal time stepping
			std::cout << "use waveExampleCubical3D_TimeStepping for order=1\n";
			return;
		}
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double eps = 1.0e-12;
		std::function<double(Vector3)> testFunction{ [&](Vector3 p) -> double { return std::exp(0.25 * p.x) * std::sin(p.x + p.y + p.z + 1) * std::cos(p.y) * (2 - 0.003 * p.z - 0.00002 * p.z * p.z); } };
		std::function<double(Vector3)> sourceTerm{ [&](Vector3 p) -> double { return -0.00001 * std::exp(0.25 * p.x) * (((-4 * p.z - 600) * p.z + 400000) * std::sin(p.y) * std::cos(p.z + p.x + p.y + 1)
			+ std::cos(p.y) * (((-3.875 * p.z - 581.25) * p.z + 387496) * std::sin(p.z + p.x + p.y + 1) + (p.z * (p.z + 142) - 100600) * std::cos(p.z + p.x + p.y + 1))); } };
		std::function<Vector3(Vector3)> testFunctionGradient{ [&](Vector3 p) -> Vector3 { return std::exp(0.25 * p.x) * Vector3(0.25 * (2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.y) * (4 * std::cos(p.x + p.y + p.z + 1) + std::sin(p.x + p.y + p.z + 1)),
			(2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + 2 * p.y + p.z + 1), std::cos(p.y) * ((2 - 0.003 * p.z - 0.00002 * p.z * p.z) * std::cos(p.x + p.y + p.z + 1) - (0.003 + 0.00004 * p.z) * std::sin(p.x + p.y + p.z + 1))); } };

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleSmallSteps3D_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		uint maxMeshNumber = 7;
		if (order == 2)
			maxMeshNumber = 6;
		else if (order == 3 || order == 4)
			maxMeshNumber = 5;
		else if (order >= 5)
			maxMeshNumber = 4;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) { //solve the problem on consecutively refined meshes
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			//solve the problem using cubical forms of chosen order
			std::cout << "DEC solution of order " << order << '\n';
			output << "DEC solution of order " << order << '\n';

			//Create one time step of the mesh in [0, Lx] x [0, Ly] and its refinement and mark cells with appropriate flags.
			uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;
			const uint totalSteps = order * refs * T_REPEATS - order + 1;
			BuilderMesh mesh_old(3);
			createCartesianMesh(mesh_old, Lx, Ly, dtime, refs, refs, 1);
			BuilderMesh mesh(3);
			SmallCubePartition3D scp(order);
			scp.refineMesh(mesh_old, mesh);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				Vector3 p = mesh.getNodePosition3(i);
				if (p.x < eps || p.x > Lx - eps || p.y < eps || p.y > Ly - eps)
					mesh.setNodeFlag(i, 1);
				else if (p.z < eps)
					mesh.setNodeFlag(i, 2);
				else if (p.z > dtime - eps)
					mesh.setNodeFlag(i, 3);
				else
					mesh.setNodeFlag(i, 0);
			}

			//find indices of the nodes in the refined mesh in a natural order
			uint Nx_ref = refs * order;
			uint Ny_ref = refs * order;
			uint Nt_ref = order;
			uint timeStepNodes = (Nx_ref + 1) * (Ny_ref + 1);
			uint activeNodes = (Nx_ref - 1) * (Ny_ref - 1);
			Buffer<Buffer<uint>> indices(order + 1);
			for (uint i = 0; i < indices.size(); ++i)
				indices[i].resize(timeStepNodes);
			double edgeLen_X = Lx / Nx_ref;
			double edgeLen_Y = Ly / Ny_ref;
			double edgeLen_T = dtime / Nt_ref;
			Vector4 pos(0.0, 0.0, 0.0, 0.0);
			uint prev = 0;
			for (uint i = 0; i <= order; ++i) {
				for (uint j = 0; j < Ny_ref + 1; ++j) {
					for (uint k = 0; k < Nx_ref + 1; ++k) {
						indices[i][j * (Ny_ref + 1) + k] = prev = mesh.findNode(pos, eps, prev);
						pos.x += edgeLen_X;
					}
					pos.x = 0.0;
					pos.y += edgeLen_Y;
				}
				pos.y = 0.0;
				pos.z += edgeLen_T;
			}

			Buffer<Buffer<double>> system_matrix(mesh.getNodeSize()); //form the system matrix
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				system_matrix[i] = Buffer<double>(mesh.getNodeSize(), 0.0);
			}
			{
				Buffer<std::unordered_map<uint, double>> star;
				scp.formHodgeMatrix1Forms(star, true);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//prepare memory for the full solution and for indices and values needed in updating the solution
			Buffer<Buffer<double>> solFull(totalSteps);
			for (uint i = 0; i < solFull.size(); ++i)
				solFull[i].resize(activeNodes);
			Buffer<uint> activeNodeIndices(activeNodes);
			Buffer<uint> boundaryNodeIndices(timeStepNodes - activeNodes);
			for (uint i = 0, j = 0, k = 0; i < timeStepNodes; ++i) {
				if (mesh.getNodeFlag(indices[order][i]) == 1)
					boundaryNodeIndices[j++] = indices[order][i];
				else
					activeNodeIndices[k++] = indices[order][i];
			}
			Buffer<uint> activeDualIndices(activeNodes);
			for (uint i = 0, j = 0; i < timeStepNodes; ++i) {
				if (mesh.getNodeFlag(indices[order - 1][i]) != 1)
					activeDualIndices[j++] = indices[order - 1][i];
			}
			Buffer<Buffer<double>> previousValues(order);
			for (uint i = 0; i < order; ++i)
				previousValues[i].resize(timeStepNodes);
			Buffer<double> boundaryValues(boundaryNodeIndices.size());
			Buffer<Buffer<Buffer<double>>> updateValues(order);
			for (uint i = 0; i < order; ++i) {
				updateValues[i].resize(activeNodes);
				for (uint j = 0; j < activeNodes; ++j) {
					updateValues[i][j].resize(timeStepNodes);
					uint index_j = activeDualIndices[j];
					for (uint k = 0; k < timeStepNodes; ++k) {
						updateValues[i][j][k] = system_matrix[index_j][indices[i][k]];
					}
				}
			}
			Buffer<Buffer<double>> boundaryUpdateValues(activeNodes);
			for (uint j = 0; j < activeNodes; ++j) {
				boundaryUpdateValues[j].resize(boundaryNodeIndices.size());
				uint index_j = activeDualIndices[j];
				for (uint k = 0; k < boundaryNodeIndices.size(); ++k) {
					boundaryUpdateValues[j][k] = system_matrix[index_j][boundaryNodeIndices[k]];
				}
			}
			Buffer<Buffer<double>> system_matrix_restricted(activeNodes);
			for (uint i = 0; i < activeDualIndices.size(); ++i) {
				system_matrix_restricted[i].resize(activeNodes);
				uint index_i = activeDualIndices[i];
				for (uint j = 0; j < activeNodes; ++j) {
					uint index_j = activeNodeIndices[j];
					system_matrix_restricted[i][j] = system_matrix[index_i][index_j];

				}
			}
			Buffer<uint> p;
			decomposeLUP(system_matrix_restricted, p, 1e-14);

			//initialise values
			for (uint i = 0; i < order; ++i) {
				for (uint j = 0; j < timeStepNodes; ++j) {
					previousValues[i][j] = testFunction(mesh.getNodePosition3(indices[i][j])); //assume known values to initialise the scheme (cheating)
				}
			}
			for (uint i = 0; i < boundaryNodeIndices.size(); ++i)
				boundaryValues[i] = testFunction(mesh.getNodePosition3(boundaryNodeIndices[i]));
			Buffer<double> b(activeNodes, 0.0);
			Buffer<double> x(activeNodes, 0.0);
			Buffer<double> sol(timeStepNodes);

			for (uint smallstep = 0; smallstep < totalSteps; ++smallstep) {
				//1) Integrate source term over active dual cells
				for (uint i = 0; i < activeDualIndices.size(); ++i)
					b[i] = integrateDual3CellCartesianMesh(mesh, activeDualIndices[i], sourceTerm);

				//2) Subtract the contribution of previous values from the left-hand side
				for (uint i = 0; i < order; ++i)
					for (uint j = 0; j < activeDualIndices.size(); ++j) {
						for (uint k = 0; k < updateValues[i][j].size(); ++k) {
							b[j] -= updateValues[i][j][k] * previousValues[i][k];
						}
					}

				//3) Update the boundary values and subtract the contribution of boundary values from the left-hand side
				for (uint k = 0; k < boundaryNodeIndices.size(); ++k) {
					boundaryValues[k] = testFunction(mesh.getNodePosition3(boundaryNodeIndices[k]));
				}
				for (uint j = 0; j < activeDualIndices.size(); ++j) {
					for (uint k = 0; k < boundaryNodeIndices.size(); ++k) {
						b[j] -= boundaryUpdateValues[j][k] * boundaryValues[k];
					}
				}

				//4) Solve the values on active nodes
				solveLUP(system_matrix_restricted, p, b, x);

				//5) Get the solution from these and boundary values
				for (uint i = 0, j = 0, k = 0; i < timeStepNodes; ++i) {
					if (mesh.getNodeFlag(indices[order][i]) == 1)
						sol[i] = boundaryValues[j++];
					else
						sol[i] = x[k++];
				}

				//6) Stepwise error analysis
				double maxError = 0.0;
				double averageError = 0.0;
				for (uint i = 0; i < activeNodes; ++i) {
					double diff = std::abs(x[i] - testFunction(mesh.getNodePosition3(activeNodeIndices[i])));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
				averageError /= activeNodes;
				if ((smallstep + 1) % 10 == 0) { //choose how often to print the result
					std::cout << "Time elapsed: " << dtime + smallstep * dtime / order << '\n';
					output << "Time elapsed: " << dtime + smallstep * dtime / order << '\n';
					std::cout << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
					output << "Max error at current time step: " << maxError << " (average " << averageError << ')' << '\n';
				}

				//7) Store the values in the full solution
				for (uint i = 0; i < activeNodes; ++i) {
					solFull[smallstep][i] = x[i];
				}

				//8) Update the values
				for (uint i = 0; i < order - 1; ++i) {
					previousValues[i] = previousValues[i + 1];
				}
				for (uint i = 0; i < timeStepNodes; ++i) {
					previousValues[order - 1] = sol;
				}

				//9) Move the mesh
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					mesh.setNodePosition(i, mesh.getNodePosition(i) + Vector4(0, 0, dtime / order, 0));
				}
			}

			//error analysis for the full solution
			std::cout << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";
			output << "Full solution error analysis in time interval [0, " << T * T_REPEATS << "]:\n";

			//maximum and average error
			double maxError = 0.0;
			double averageError = 0.0;
			for (uint i = 0; i < totalSteps; ++i) {
				for (uint j = 0; j < activeNodes; ++j) {
					Vector3 pos = mesh.getNodePosition3(activeNodeIndices[j]);
					double diff = std::abs(solFull[i][j] - testFunction(Vector3(pos.x, pos.y, (order + i) * (dtime / order))));
					if (maxError < diff)
						maxError = diff;
					averageError += diff;
				}
			}
			averageError /= (totalSteps * activeNodes);
			std::cout << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			output << "Max error " << maxError << " (average " << averageError << ')' << '\n';
			std::cout << '\n';
			output << '\n';
			output.save(outputFile.str());
		}
	}

	//Computes the norms of the inverses of the system matrices of poissonExampleWhitney2D.
	void poissonExampleWhitney2D_Stability(bool loadFromFile, bool circumcentric) {
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text filepath;
		filepath << "Files/poissonExampleWhitney2D_Stability" << (circumcentric ? "_circumcentric" : "_barycentric") << ".txt";

		for (uint meshNumber = 1; meshNumber <= 7; ++meshNumber) {
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			BuilderMesh mesh_old(2);
			const double r = 2.0;
			mesh_old.createTriangleGrid(Vector2(-r, -r), Vector2(r, r), 4.0);

			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(2);
				SmallSimplexPartition2D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			for (uint order = 1; order <= 8 && (meshNumber <= 4 || order <= 4) && (meshNumber <= 5 || order <= 2) && (meshNumber <= 6 || order <= 1); ++order) {
				BuilderMesh mesh(2);
				SmallSimplexPartition2D ssp(order);
				ssp.refineMesh(mesh_old, mesh);
				uint interiorNodes = 0;
				mesh.fillBoundaryFlags(1);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) != 1)
						++interiorNodes;
				}
				MatrixN system_matrix_interior_inverse;
				if (loadFromFile && meshNumber + order > 2) {
					Text sysmat_text;
					sysmat_text << "Files/2D/order" << order << (circumcentric ? "/circumcentric" : "/barycentric") << "/system_matrix_interior_inverse" << meshNumber << ".dat";
					loadMatrix(system_matrix_interior_inverse, sysmat_text.str());
				}
				else {
					MatrixN system_matrix_interior;
					system_matrix_interior.toMatrixN(interiorNodes);
					{
						Buffer<std::unordered_map<uint, double>> star;
						ssp.formHodgeMatrix1FormsCustom(star, "DefaultTriangle", circumcentric, false);
						for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
							if (mesh.getNodeFlag(j) == 1)
								continue;
							Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
							const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
							for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
								std::unordered_map<uint, double>& row = star[i];
								for (uint e = 0; e < edges_j.size(); ++e) {
									auto it = row.find(edges_j[e]);
									if (it != row.end())
										col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
								}
							}
							for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
								if (mesh.getNodeFlag(i) == 1)
									continue;
								const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
								for (uint e = 0; e < edges_i.size(); ++e) {
									system_matrix_interior[ii][jj] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
								}
								++ii;
							}
							++jj;
						}
					}
					system_matrix_interior_inverse = system_matrix_interior.inverse();
					if (interiorNodes != 0) {
						Text sysmat_text;
						sysmat_text << "Files/2D/order" << order << (circumcentric ? "/circumcentric" : "/barycentric") << "/system_matrix_interior_inverse" << meshNumber << ".dat";
						saveMatrix(system_matrix_interior_inverse, sysmat_text.str());
					}
				}

				Buffer<double> weights(interiorNodes);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						continue;
					weights[j] = dualFaceArea(i, mesh, circumcentric);
					++j;
				}
				double norm = computeWeightedInfinityNorm(system_matrix_interior_inverse, weights);

				std::cout << "Norm of the system matrix inverse of order " << order << ": " << norm << '\n' << '\n';
				output << "Norm of the system matrix inverse of order " << order << ": " << norm << '\n' << '\n';
				output.save(filepath.str());
			}
		}
	}

	//Computes the norms of the inverses of the system matrices of poissonExampleWhitney3D.
	void poissonExampleWhitney3D_Stability(bool loadFromFile, bool circumcentric) {
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text filepath;
		filepath << "Files/poissonExampleWhitney3D_Stability" << (circumcentric ? "_circumcentric" : "_barycentric") << ".txt";

		for (uint meshNumber = 1; meshNumber <= 6; ++meshNumber) {
			std::cout << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";
			output << "### ### ### MESH " << meshNumber << " ### ### ###\n\n";

			BuilderMesh mesh_old(3);
			const double r = 2.0125;
			createRhombicDodecahedronMesh(mesh_old, r, 2.0);
			for (uint i = mesh_old.getNodeSize(); i-- > 0; )
			{
				const Vector4 p = mesh_old.getNodePosition(i);
				if (p.z < 0 - 1e-5) mesh_old.removeNode(i);
				else if (p.x > p.z + 1e-5) mesh_old.removeNode(i);
				else if (-p.x > p.z + 1e-5) mesh_old.removeNode(i);
				else if (p.y > p.z + 1e-5) mesh_old.removeNode(i);
				else if (-p.y > p.z + 1e-5) mesh_old.removeNode(i);
				else mesh_old.setNodePosition(i, Vector4(p.x, p.y, p.z - 1.0, 0.0));
			}

			for (uint refs = 1; refs < meshNumber; ++refs) {
				BuilderMesh helpmesh(3);
				SmallSimplexPartition3D refiner(2);
				refiner.refineMesh(mesh_old, helpmesh);
				mesh_old.createCopy(helpmesh);
			}

			for (uint order = 1; order <= 8 && (meshNumber <= 3 || order <= 4) && (meshNumber <= 4 || order <= 2) && (meshNumber <= 5 || order <= 1); ++order) {
				BuilderMesh mesh(3);
				SmallSimplexPartition3D ssp(order);
				ssp.refineMesh(mesh_old, mesh);
				uint interiorNodes = 0;
				mesh.fillBoundaryFlags(1);
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) != 1)
						++interiorNodes;
				}
				MatrixN system_matrix_interior_inverse;
				if (loadFromFile && meshNumber + order > 2) {
					Text sysmat_text;
					sysmat_text << "Files/3D/order" << order << (circumcentric ? "/circumcentric" : "/barycentric") << "/system_matrix_interior_inverse" << meshNumber << ".dat";
					loadMatrix(system_matrix_interior_inverse, sysmat_text.str());
				}
				else {
					MatrixN system_matrix_interior;
					system_matrix_interior.toMatrixN(interiorNodes);
					{
						Buffer<std::unordered_map<uint, double>> star;
						ssp.formHodgeMatrix1FormsCustom(star, "BccTetrahedron", circumcentric, false);
						for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
							if (mesh.getNodeFlag(j) == 1)
								continue;
							Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
							const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
							for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
								std::unordered_map<uint, double>& row = star[i];
								for (uint e = 0; e < edges_j.size(); ++e) {
									auto it = row.find(edges_j[e]);
									if (it != row.end())
										col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
								}
							}
							for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
								if (mesh.getNodeFlag(i) == 1)
									continue;
								const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
								for (uint e = 0; e < edges_i.size(); ++e) {
									system_matrix_interior[ii][jj] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
								}
								++ii;
							}
							++jj;
						}
					}
					system_matrix_interior_inverse = system_matrix_interior.inverse();
					if (interiorNodes != 0) {
						Text sysmat_text;
						sysmat_text << "Files/3D/order" << order << (circumcentric ? "/circumcentric" : "/barycentric") << "/system_matrix_interior_inverse" << meshNumber << ".dat";
						saveMatrix(system_matrix_interior_inverse, sysmat_text.str());
					}
				}

				Buffer<double> weights(interiorNodes);
				for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
					if (mesh.getNodeFlag(i) == 1)
						continue;
					weights[j] = dualCellVolume(i, mesh, circumcentric);
					++j;
				}
				double norm = computeWeightedInfinityNorm(system_matrix_interior_inverse, weights);

				std::cout << "Norm of the system matrix inverse of order " << order << ": " << norm << '\n' << '\n';
				output << "Norm of the system matrix inverse of order " << order << ": " << norm << '\n' << '\n';
				output.save(filepath.str());
			}
		}
	}

	//Computes the norm of the inverse of the system matrix in waveExampleWhitney2D and waveExampleWhitney2D_TimeStepping.
	void waveExampleWhitney2D_Stability(uint order, const double T, const uint T_REPEATS) {
		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		const double L = 2.0;
		Text outputFile;
		outputFile << "Files/waveExampleWhitney2D_Stability_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		std::cout << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';
		output << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';

		uint maxMeshNumber = 8;
		if (order == 2)
			maxMeshNumber = 7;
		else if (order == 3 || order == 4 || order == 5)
			maxMeshNumber = 6;
		else if (order >= 6)
			maxMeshNumber = 5;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) {
			//build three time steps of the mesh that has been refined meshNumber - 1 times
			const uint refs = std::pow(2, meshNumber - 1);
			const double dtime = 3 * T / refs;
			BuilderMesh mesh_old(2);
			createSpacetimeMesh(mesh_old, L, dtime, refs, 3);
			BuilderMesh mesh(2);
			SmallSimplexPartition2D ssp(order);
			ssp.refineMesh(mesh_old, mesh);

			//find indices of the nodes in the refined mesh in a natural order
			uint Nx_ref = refs * order;
			uint Nt_ref = 3 * order;
			uint activeNodes = (Nx_ref - 1) * Nt_ref;
			uint boundary2Nodes = Nx_ref - 1;
			Buffer<uint> indices((Nx_ref - 1) * (Nt_ref + 1));
			double edgeLen_X = L / Nx_ref;
			double edgeLen_T = dtime / Nt_ref;
			Vector4 pos(0.0, 0.0, 0.0, 0.0);
			uint prev = 0;
			for (uint i = 0; i <= Nt_ref; ++i) {
				for (uint j = 0; j < Nx_ref - 1; ++j) {
					pos.x += edgeLen_X;
					indices[i * (Nx_ref - 1) + j] = prev = mesh.findNode(pos, 1.0e-12, prev);
				}
				pos.x = 0.0;
				pos.y += edgeLen_T;
			}

			//form the blocks A, B, and C of the permuted system matrix
			MatrixN blockA, blockB, blockC;
			uint blockSize = boundary2Nodes * order;
			blockA.toMatrixN(blockSize);
			blockB.toMatrixN(blockSize);
			blockC.toMatrixN(blockSize);
			{
				Buffer<std::unordered_map<uint, double>> star;
				Buffer<uint> refElementIndices(1, 0);
				Buffer <Buffer<uint>> refElements(1);
				refElements[0].resize(mesh_old.getFaceSize());
				for (uint i = 0; i < mesh_old.getFaceSize(); ++i)
					refElements[0][i] = i;
				ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "RightTriangle", false, true, false);
				for (uint j = 0; j < activeNodes; ++j) {
					uint index_j = indices[j + boundary2Nodes];
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(index_j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], index_j);
						}
					}
					for (uint i = 0; i < blockSize; ++i) {
						uint index_i = indices[i + 2 * blockSize];
						const Buffer<uint>& edges_i = mesh.getNodeEdges(index_i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								if (j < blockSize)
									blockC[i][j] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
								else if (j < 2 * blockSize)
									blockB[i][j - blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
								else
									blockA[i][j - 2 * blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//compute the blocks A_i of the inverse and the norm of the inverse
			MatrixN blockA_inv = blockA.inverse();
			MatrixN A0, A1, A2;
			A0 = blockA_inv;
			A1 = -blockA_inv * blockB * blockA_inv;
			A2 = -blockA_inv * (blockB * A1 + blockC * A0);
			Buffer<double> weights0(blockSize);
			Buffer<double> weights1(blockSize);
			for (uint i = 0; i < blockSize; ++i) {
				weights0[i] = dualFaceArea(indices[i], mesh, false);
				weights1[i] = dualFaceArea(indices[i + blockSize], mesh, false);
			}
			Buffer<double> weightedRowSums(blockSize, 0.0);
			Buffer<double> rowSums(blockSize, 0.0);
			for (uint i = 0; i < blockSize; ++i) {
				for (uint j = 0; j < blockSize; ++j) {
					rowSums[i] += std::abs(A0[i][j]);
					weightedRowSums[i] += std::abs(A0[i][j]) * weights1[j];
					if (T_REPEATS * refs > 1) {
						rowSums[i] += std::abs(A1[i][j]);
						weightedRowSums[i] += std::abs(A1[i][j]) * weights1[j];
					}
					if (T_REPEATS * refs > 2) {
						rowSums[i] += std::abs(A2[i][j]);
						weightedRowSums[i] += std::abs(A2[i][j]) * weights1[j];
					}					
				}
			}
			for (uint k = 4; k <= T_REPEATS * refs; ++k) {
				MatrixN temp = -blockA_inv * (blockB * A2 + blockC * A1);
				A0 = A1;
				A1 = A2;
				A2 = temp;
				for (uint i = 0; i < blockSize; ++i) {
					for (uint j = 0; j < blockSize; ++j) {
						rowSums[i] += std::abs(A2[i][j]);
						if (k == T_REPEATS * refs)
							weightedRowSums[i] += std::abs(A2[i][j]) * weights0[j];
						else
							weightedRowSums[i] += std::abs(A2[i][j]) * weights1[j];
					}
				}
				if (k % 10 == 0) { //choose how often to print the result
					double norm = 0.0;
					double weightedNorm = 0.0;

					for (uint i = 0; i < blockSize; ++i) {
						if (rowSums[i] > norm)
							norm = rowSums[i];
						if (weightedRowSums[i] > weightedNorm)
							weightedNorm = weightedRowSums[i];
					}
					std::cout << "Time elapsed: " << k * T / refs << " --- System matrix inverse norm: " << weightedNorm << '\n';
					output << "Time elapsed: " << k * T / refs << " --- System matrix inverse norm: " << weightedNorm << '\n';
				}
			}
			double norm = 0.0;
			double weightedNorm = 0.0;

			for (uint i = 0; i < blockSize; ++i) {
				if (rowSums[i] > norm)
					norm = rowSums[i];
				if (weightedRowSums[i] > weightedNorm)
					weightedNorm = weightedRowSums[i];
			}
			std::cout << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output.save(outputFile.str());
		}
	}

	//Computes the norm of the inverse of the system matrix in waveExampleWhitney3D and waveExampleWhitney3D_TimeStepping.
	void waveExampleWhitney3D_Stability(uint order, const double T, const uint T_REPEATS) {
		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		const double Lx = 2.0;
		const double Ly = 2.0;
		Text outputFile;
		outputFile << "Files/waveExampleWhitney3D_Stability_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		std::cout << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';
		output << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';

		uint maxMeshNumber = 6;
		if (order == 2)
			maxMeshNumber = 5;
		else if (order == 3 || order == 4)
			maxMeshNumber = 4;
		else if (order >= 5)
			maxMeshNumber = 3;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) {
			//build three time steps of the mesh that has been refined meshNumber times
			const uint refs = std::pow(2, meshNumber - 1);
			const double dtime = 3 * T / refs;
			BuilderMesh mesh_old(3);
			createSpacetimeMesh(mesh_old, Lx, Ly, dtime, refs, refs, 3);
			BuilderMesh mesh(3);
			SmallSimplexPartition3D ssp(order);
			ssp.refineMesh(mesh_old, mesh);

			//find indices of the nodes in the refined mesh in a natural order
			uint Nx_ref = refs * order;
			uint Ny_ref = refs * order;
			uint Nt_ref = 3 * order;
			uint activeNodes = (Nx_ref - 1) * (Ny_ref - 1) * Nt_ref;
			uint boundary2Nodes = (Nx_ref - 1) * (Ny_ref - 1);
			Buffer<uint> indices((Nx_ref - 1) * (Ny_ref - 1) * (Nt_ref + 1));
			double edgeLen_X = Lx / Nx_ref;
			double edgeLen_Y = Ly / Ny_ref;
			double edgeLen_T = dtime / Nt_ref;
			Vector4 pos(0.0, 0.0, 0.0, 0.0);
			uint prev = 0;
			for (uint i = 0; i <= Nt_ref; ++i) {
				for (uint j = 0; j < Ny_ref - 1; ++j) {
					pos.y += edgeLen_Y;
					for (uint k = 0; k < Nx_ref - 1; ++k) {
						pos.x += edgeLen_X;
						indices[i * (Nx_ref - 1) * (Ny_ref - 1) + j * (Ny_ref - 1) + k] = prev = mesh.findNode(pos, 1.0e-12, prev);
					}
					pos.x = 0.0;
				}
				pos.y = 0.0;
				pos.z += edgeLen_T;
			}

			//form the blocks A, B, and C of the permuted system matrix
			MatrixN blockA, blockB, blockC;
			uint blockSize = boundary2Nodes * order;
			blockA.toMatrixN(blockSize);
			blockB.toMatrixN(blockSize);
			blockC.toMatrixN(blockSize);
			{
				Buffer<std::unordered_map<uint, double>> star;
				Buffer<uint> refElementIndices(6);
				Buffer<Buffer<uint>> refElements(6);
				for (uint i = 0; i < 6; ++i) {
					refElementIndices[i] = i;
					refElements[i].resize(mesh_old.getBodySize() / 6);
				}

				for (uint i = 0; i < mesh_old.getBodySize(); ++i)
					refElements[i % 6][i / 6] = i;
				ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "CubeTetrahedron", false, true, false);
				for (uint j = 0; j < activeNodes; ++j) {
					uint index_j = indices[j + boundary2Nodes];
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(index_j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], index_j);
						}
					}
					for (uint i = 0; i < blockSize; ++i) {
						uint index_i = indices[i + 2 * blockSize];
						const Buffer<uint>& edges_i = mesh.getNodeEdges(index_i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								if (j < blockSize)
									blockC[i][j] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
								else if (j < 2 * blockSize)
									blockB[i][j - blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
								else
									blockA[i][j - 2 * blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//compute the blocks A_i of the inverse and the norm of the inverse
			MatrixN blockA_inv = blockA.inverse();
			MatrixN A0, A1, A2;
			A0 = blockA_inv;
			A1 = -blockA_inv * blockB * blockA_inv;
			A2 = -blockA_inv * (blockB * A1 + blockC * A0);
			Buffer<double> weights0(blockSize);
			Buffer<double> weights1(blockSize);
			for (uint i = 0; i < blockSize; ++i) {
				weights0[i] = dualCellVolume(indices[i], mesh, false);
				weights1[i] = dualCellVolume(indices[i + blockSize], mesh, false);
			}
			Buffer<double> weightedRowSums(blockSize, 0.0);
			Buffer<double> rowSums(blockSize, 0.0);
			for (uint i = 0; i < blockSize; ++i) {
				for (uint j = 0; j < blockSize; ++j) {
					rowSums[i] += std::abs(A0[i][j]);
					weightedRowSums[i] += std::abs(A0[i][j]) * weights1[j];
					if (T_REPEATS * refs > 1) {
						rowSums[i] += std::abs(A1[i][j]);
						weightedRowSums[i] += std::abs(A1[i][j]) * weights1[j];
					}
					if (T_REPEATS * refs > 2) {
						rowSums[i] += std::abs(A2[i][j]);
						weightedRowSums[i] += std::abs(A2[i][j]) * weights1[j];
					}
				}
			}
			for (uint k = 4; k <= T_REPEATS * refs; ++k) {
				MatrixN temp = -blockA_inv * (blockB * A2 + blockC * A1);
				A0 = A1;
				A1 = A2;
				A2 = temp;
				for (uint i = 0; i < blockSize; ++i) {
					for (uint j = 0; j < blockSize; ++j) {
						rowSums[i] += std::abs(A2[i][j]);
						if (k == T_REPEATS * refs)
							weightedRowSums[i] += std::abs(A2[i][j]) * weights0[j];
						else
							weightedRowSums[i] += std::abs(A2[i][j]) * weights1[j];
					}
				}
				if (k % 10 == 0) { //choose how often to print the result
					double norm = 0.0;
					double weightedNorm = 0.0;

					for (uint i = 0; i < blockSize; ++i) {
						if (rowSums[i] > norm)
							norm = rowSums[i];
						if (weightedRowSums[i] > weightedNorm)
							weightedNorm = weightedRowSums[i];
					}
					std::cout << "Time elapsed: " << k * T / refs << " --- System matrix inverse norm: " << weightedNorm << '\n';
					output << "Time elapsed: " << k * T / refs << " --- System matrix inverse norm: " << weightedNorm << '\n';
				}
			}
			double norm = 0.0;
			double weightedNorm = 0.0;

			for (uint i = 0; i < blockSize; ++i) {
				if (rowSums[i] > norm)
					norm = rowSums[i];
				if (weightedRowSums[i] > weightedNorm)
					weightedNorm = weightedRowSums[i];
			}
			std::cout << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output.save(outputFile.str());
		}
	}

	//same as waveExampleWhitney2D_Stability but with cubical forms
	void waveExampleCubical2D_Stability(uint order, const double T, const uint T_REPEATS) {
		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		const double L = 2.0;
		Text outputFile;
		outputFile << "Files/waveExampleCubical2D_Stability_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		std::cout << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';
		output << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';

		uint maxMeshNumber = 8;
		if (order == 2)
			maxMeshNumber = 7;
		else if (order == 3 || order == 4)
			maxMeshNumber = 6;
		else if (order >= 5)
			maxMeshNumber = 5;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) {
			//build three time steps of the mesh that has been refined meshNumber times
			const uint refs = std::pow(2, meshNumber - 1);
			const double dtime = 3 * T / refs;
			BuilderMesh mesh_old(2);
			createCartesianMesh(mesh_old, L, dtime, refs, 3);
			BuilderMesh mesh(2);
			SmallCubePartition2D scp(order);
			scp.refineMesh(mesh_old, mesh);

			//find indices of the nodes in the refined mesh in a natural order
			uint Nx_ref = refs * order;
			uint Nt_ref = 3 * order;
			uint activeNodes = (Nx_ref - 1) * Nt_ref;
			uint boundary2Nodes = Nx_ref - 1;
			Buffer<uint> indices((Nx_ref - 1) * (Nt_ref + 1));
			double edgeLen_X = L / Nx_ref;
			double edgeLen_T = dtime / Nt_ref;
			Vector4 pos(0.0, 0.0, 0.0, 0.0);
			uint prev = 0;
			for (uint i = 0; i <= Nt_ref; ++i) {
				for (uint j = 0; j < Nx_ref - 1; ++j) {
					pos.x += edgeLen_X;
					indices[i * (Nx_ref - 1) + j] = prev = mesh.findNode(pos, 1.0e-12, prev);
				}
				pos.x = 0.0;
				pos.y += edgeLen_T;
			}

			//form the blocks A, B, and C of the permuted system matrix
			MatrixN blockA, blockB, blockC;
			uint blockSize = boundary2Nodes * order;
			blockA.toMatrixN(blockSize);
			blockB.toMatrixN(blockSize);
			blockC.toMatrixN(blockSize);
			{
				Buffer<std::unordered_map<uint, double>> star;
				scp.formHodgeMatrix1Forms(star, true);
				for (uint j = 0; j < activeNodes; ++j) {
					uint index_j = indices[j + boundary2Nodes];
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(index_j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], index_j);
						}
					}
					for (uint i = 0; i < blockSize; ++i) {
						uint index_i = indices[i + 2 * blockSize];
						const Buffer<uint>& edges_i = mesh.getNodeEdges(index_i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								if (j < blockSize)
									blockC[i][j] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
								else if (j < 2 * blockSize)
									blockB[i][j - blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
								else
									blockA[i][j - 2 * blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//compute the blocks A_i of the inverse and the norm of the inverse
			MatrixN blockA_inv = blockA.inverse();
			MatrixN A0, A1, A2;
			A0 = blockA_inv;
			A1 = -blockA_inv * blockB * blockA_inv;
			A2 = -blockA_inv * (blockB * A1 + blockC * A0);
			Buffer<double> weights0(blockSize);
			Buffer<double> weights1(blockSize);
			for (uint i = 0; i < blockSize; ++i) {
				weights0[i] = dualFaceArea(indices[i], mesh, false);
				weights1[i] = dualFaceArea(indices[i + blockSize], mesh, false);
			}
			Buffer<double> weightedRowSums(blockSize, 0.0);
			Buffer<double> rowSums(blockSize, 0.0);
			for (uint i = 0; i < blockSize; ++i) {
				for (uint j = 0; j < blockSize; ++j) {
					rowSums[i] += std::abs(A0[i][j]);
					weightedRowSums[i] += std::abs(A0[i][j]) * weights1[j];
					if (T_REPEATS * refs > 1) {
						rowSums[i] += std::abs(A1[i][j]);
						weightedRowSums[i] += std::abs(A1[i][j]) * weights1[j];
					}
					if (T_REPEATS * refs > 2) {
						rowSums[i] += std::abs(A2[i][j]);
						weightedRowSums[i] += std::abs(A2[i][j]) * weights1[j];
					}
				}
			}
			for (uint k = 4; k <= T_REPEATS * refs; ++k) {
				MatrixN temp = -blockA_inv * (blockB * A2 + blockC * A1);
				A0 = A1;
				A1 = A2;
				A2 = temp;
				for (uint i = 0; i < blockSize; ++i) {
					for (uint j = 0; j < blockSize; ++j) {
						rowSums[i] += std::abs(A2[i][j]);
						if (k == T_REPEATS * refs)
							weightedRowSums[i] += std::abs(A2[i][j]) * weights0[j];
						else
							weightedRowSums[i] += std::abs(A2[i][j]) * weights1[j];
					}
				}
				if (k % 10 == 0) { //choose how often to print the result
					double norm = 0.0;
					double weightedNorm = 0.0;

					for (uint i = 0; i < blockSize; ++i) {
						if (rowSums[i] > norm)
							norm = rowSums[i];
						if (weightedRowSums[i] > weightedNorm)
							weightedNorm = weightedRowSums[i];
					}
					std::cout << "Time elapsed: " << k * T / refs << " --- System matrix inverse norm: " << weightedNorm << '\n';
					output << "Time elapsed: " << k * T / refs << " --- System matrix inverse norm: " << weightedNorm << '\n';
				}
			}
			double norm = 0.0;
			double weightedNorm = 0.0;

			for (uint i = 0; i < blockSize; ++i) {
				if (rowSums[i] > norm)
					norm = rowSums[i];
				if (weightedRowSums[i] > weightedNorm)
					weightedNorm = weightedRowSums[i];
			}
			std::cout << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output.save(outputFile.str());
		}
	}

	//same as waveExampleWhitney3D_Stability but with cubical forms
	void waveExampleCubical3D_Stability(uint order, const double T, const uint T_REPEATS) {
		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;

		const double Lx = 2.0;
		const double Ly = 2.0;
		Text outputFile;
		outputFile << "Files/waveExampleCubical3D_Stability_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		std::cout << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';
		output << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';

		uint maxMeshNumber = 6;
		if (order == 2)
			maxMeshNumber = 5;
		else if (order == 3 || order == 4)
			maxMeshNumber = 4;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) {
			//build three time steps of the mesh that has been refined meshNumber times
			const uint refs = std::pow(2, meshNumber - 1);
			const double dtime = 3 * T / refs;
			BuilderMesh mesh_old(3);
			createCartesianMesh(mesh_old, Lx, Ly, dtime, refs, refs, 3);
			BuilderMesh mesh(3);
			SmallCubePartition3D scp(order);
			scp.refineMesh(mesh_old, mesh);

			//find indices of the nodes in the refined mesh in a natural order
			uint Nx_ref = refs * order;
			uint Ny_ref = refs * order;
			uint Nt_ref = 3 * order;
			uint activeNodes = (Nx_ref - 1) * (Ny_ref - 1) * Nt_ref;
			uint boundary2Nodes = (Nx_ref - 1) * (Ny_ref - 1);
			Buffer<uint> indices((Nx_ref - 1) * (Ny_ref - 1) * (Nt_ref + 1));
			double edgeLen_X = Lx / Nx_ref;
			double edgeLen_Y = Ly / Ny_ref;
			double edgeLen_T = dtime / Nt_ref;
			Vector4 pos(0.0, 0.0, 0.0, 0.0);
			uint prev = 0;
			for (uint i = 0; i <= Nt_ref; ++i) {
				for (uint j = 0; j < Ny_ref - 1; ++j) {
					pos.y += edgeLen_Y;
					for (uint k = 0; k < Nx_ref - 1; ++k) {
						pos.x += edgeLen_X;
						indices[i * (Nx_ref - 1) * (Ny_ref - 1) + j * (Ny_ref - 1) + k] = prev = mesh.findNode(pos, 1.0e-12, prev);
					}
					pos.x = 0.0;
				}
				pos.y = 0.0;
				pos.z += edgeLen_T;
			}

			//form the blocks A, B, and C of the permuted system matrix
			MatrixN blockA, blockB, blockC;
			uint blockSize = boundary2Nodes * order;
			blockA.toMatrixN(blockSize);
			blockB.toMatrixN(blockSize);
			blockC.toMatrixN(blockSize);
			{
				Buffer<std::unordered_map<uint, double>> star;
				scp.formHodgeMatrix1Forms(star, true);
				for (uint j = 0; j < activeNodes; ++j) {
					uint index_j = indices[j + boundary2Nodes];
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(index_j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], index_j);
						}
					}
					for (uint i = 0; i < blockSize; ++i) {
						uint index_i = indices[i + 2 * blockSize];
						const Buffer<uint>& edges_i = mesh.getNodeEdges(index_i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								if (j < blockSize)
									blockC[i][j] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
								else if (j < 2 * blockSize)
									blockB[i][j - blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
								else
									blockA[i][j - 2 * blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

			//compute the blocks A_i of the inverse
			MatrixN blockA_inv = blockA.inverse();
			MatrixN A0, A1, A2;
			A0 = blockA_inv;
			A1 = -blockA_inv * blockB * blockA_inv;
			A2 = -blockA_inv * (blockB * A1 + blockC * A0);
			Buffer<double> weights0(blockSize);
			Buffer<double> weights1(blockSize);
			for (uint i = 0; i < blockSize; ++i) {
				weights0[i] = dualCellVolumeCartesianMesh(indices[i], mesh);
				weights1[i] = dualCellVolumeCartesianMesh(indices[i + blockSize], mesh);
			}
			Buffer<double> weightedRowSums(blockSize, 0.0);
			Buffer<double> rowSums(blockSize, 0.0);
			for (uint i = 0; i < blockSize; ++i) {
				for (uint j = 0; j < blockSize; ++j) {
					rowSums[i] += std::abs(A0[i][j]);
					weightedRowSums[i] += std::abs(A0[i][j]) * weights1[j];
					if (T_REPEATS * refs > 1) {
						rowSums[i] += std::abs(A1[i][j]);
						weightedRowSums[i] += std::abs(A1[i][j]) * weights1[j];
					}
					if (T_REPEATS * refs > 2) {
						rowSums[i] += std::abs(A2[i][j]);
						weightedRowSums[i] += std::abs(A2[i][j]) * weights1[j];
					}
				}
			}
			for (uint k = 4; k <= T_REPEATS * refs; ++k) {
				MatrixN temp = -blockA_inv * (blockB * A2 + blockC * A1);
				A0 = A1;
				A1 = A2;
				A2 = temp;
				for (uint i = 0; i < blockSize; ++i) {
					for (uint j = 0; j < blockSize; ++j) {
						rowSums[i] += std::abs(A2[i][j]);
						if (k == T_REPEATS * refs)
							weightedRowSums[i] += std::abs(A2[i][j]) * weights0[j];
						else
							weightedRowSums[i] += std::abs(A2[i][j]) * weights1[j];
					}
				}
				if (k % 10 == 0) { //choose how often to print the result
					double norm = 0.0;
					double weightedNorm = 0.0;

					for (uint i = 0; i < blockSize; ++i) {
						if (rowSums[i] > norm)
							norm = rowSums[i];
						if (weightedRowSums[i] > weightedNorm)
							weightedNorm = weightedRowSums[i];
					}
					std::cout << "Time elapsed: " << k * T / refs << " --- System matrix inverse norm: " << weightedNorm << '\n';
					output << "Time elapsed: " << k * T / refs << " --- System matrix inverse norm: " << weightedNorm << '\n';
				}
			}

			//compute the norm of the inverse
			double norm = 0.0;
			double weightedNorm = 0.0;

			for (uint i = 0; i < blockSize; ++i) {
				if (rowSums[i] > norm)
					norm = rowSums[i];
				if (weightedRowSums[i] > weightedNorm)
					weightedNorm = weightedRowSums[i];
			}
			std::cout << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output.save(outputFile.str());
		}
	}

	//Compute norm for the inverse of the matrix that can be interpreted as the system matrix in waveExampleSmallSteps3D.
	void waveExampleSmallSteps3D_Stability(uint order, const double T, const uint T_REPEATS) {
		if (order == 1) { //method coincides with normal time stepping
			std::cout << "use waveExampleCubical3D_Stability for order=1\n";
			return;
		}

		//choose format for displaying results
		std::cout.precision(5);
		std::cout << std::scientific;
		Text output;
		output.precision(5);
		output << std::scientific;
		Text outputFile;
		outputFile << "Files/waveExampleSmallSteps3D_Stability_order" << order << "_time" << T << "_repeats" << T_REPEATS << ".txt";

		const double Lx = 2.0;
		const double Ly = 2.0;
		const double eps = 1.0e-12;

		std::cout << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';
		output << "Order=" << order << ", time=" << T << ", T_REPEATS=" << T_REPEATS << '\n';

		uint maxMeshNumber = 7;
		if (order == 2)
			maxMeshNumber = 6;
		else if (order == 3 || order == 4)
			maxMeshNumber = 5;
		else if (order >= 5)
			maxMeshNumber = 4;
		for (uint meshNumber = 1; meshNumber <= maxMeshNumber; ++meshNumber) {
			//Create one time step of the mesh in [0, Lx] x [0, Ly] and its refinement and mark cells with appropriate flags.
			const uint refs = std::pow(2, meshNumber - 1);
			const double dtime = T / refs;
			const uint totalSteps = order * refs * T_REPEATS - order + 1;
			const double smallstepsize = dtime / order;
			BuilderMesh mesh_old(3);
			createCartesianMesh(mesh_old, Lx, Ly, dtime, refs, refs, 1);
			BuilderMesh mesh(3);
			SmallCubePartition3D ssp(order);
			ssp.refineMesh(mesh_old, mesh);
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				Vector3 p = mesh.getNodePosition3(i);
				if (p.x < eps || p.x > Lx - eps || p.y < eps || p.y > Ly - eps)
					mesh.setNodeFlag(i, 1);
				else if (p.z < eps)
					mesh.setNodeFlag(i, 2);
				else if (p.z > dtime - eps)
					mesh.setNodeFlag(i, 3);
				else
					mesh.setNodeFlag(i, 0);
			}

			//find indices of the nodes in the refined mesh in a natural order
			uint Nx_ref = refs * order;
			uint Ny_ref = refs * order;
			uint Nt_ref = order;
			uint timeStepNodes = (Nx_ref + 1) * (Ny_ref + 1);
			uint activeNodes = (Nx_ref - 1) * (Ny_ref - 1);
			Buffer<Buffer<uint>> indices(order + 1);
			for (uint i = 0; i < indices.size(); ++i)
				indices[i].resize(timeStepNodes);
			double edgeLen_X = Lx / Nx_ref;
			double edgeLen_Y = Ly / Ny_ref;
			double edgeLen_T = dtime / Nt_ref;
			Vector4 pos(0.0, 0.0, 0.0, 0.0);
			uint prev = 0;
			for (uint i = 0; i <= order; ++i) {
				for (uint j = 0; j < Ny_ref + 1; ++j) {
					for (uint k = 0; k < Nx_ref + 1; ++k) {
						indices[i][j * (Ny_ref + 1) + k] = prev = mesh.findNode(pos, 1.0e-13, prev);
						pos.x += edgeLen_X;
					}
					pos.x = 0.0;
					pos.y += edgeLen_Y;
				}
				pos.y = 0.0;
				pos.z += edgeLen_T;
			}

			Buffer<Buffer<double>> system_matrix(mesh.getNodeSize()); //form the system matrix
			for (uint i = 0; i < mesh.getNodeSize(); ++i) {
				system_matrix[i] = Buffer<double>(mesh.getNodeSize(), 0.0);
			}
			{
				Buffer<std::unordered_map<uint, double>> star;
				ssp.formHodgeMatrix1Forms(star, true);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}
			Buffer<uint> activeNodeIndices(activeNodes);
			for (uint i = 0, j = 0; i < timeStepNodes; ++i) {
				if (mesh.getNodeFlag(indices[order][i]) != 1)
					activeNodeIndices[j++] = indices[order][i];
			}
			Buffer<uint> activeDualIndices(activeNodes);
			for (uint i = 0, j = 0; i < timeStepNodes; ++i) {
				if (mesh.getNodeFlag(indices[order - 1][i]) != 1)
					activeDualIndices[j++] = indices[order - 1][i];
			}
			Buffer<Buffer<Buffer<double>>> updateValues(order);
			for (uint i = 0; i < order; ++i) {
				updateValues[i].resize(activeNodes);
				for (uint j = 0; j < activeNodes; ++j) {
					updateValues[i][j].resize(activeNodes);
					uint index_j = activeDualIndices[j];
					for (uint k = 0, l = 0; k < timeStepNodes; ++k) {
						if (mesh.getNodeFlag(indices[i][k]) == 1)
							continue;
						updateValues[i][j][l++] = system_matrix[index_j][indices[i][k]];
					}
				}
			}
			Buffer<Buffer<double>> system_matrix_restricted(activeNodes);
			for (uint i = 0; i < activeDualIndices.size(); ++i) {
				system_matrix_restricted[i].resize(activeNodes);
				uint index_i = activeDualIndices[i];
				for (uint j = 0; j < activeNodes; ++j) {
					uint index_j = activeNodeIndices[j];
					system_matrix_restricted[i][j] = system_matrix[index_i][index_j];

				}
			}

			//compute the blocks A_i of the inverse
			MatrixN A(system_matrix_restricted);
			Buffer<MatrixN> B(order);
			for (uint i = 0; i < order; ++i)
				B[i] = MatrixN(updateValues[order - i - 1]);
			MatrixN A_inv = A.inverse();
			Buffer<MatrixN> A_i(order);
			A_i[0] = A_inv;
			for (uint i = 1; i < order; ++i) {
				MatrixN sum = B[0] * A_i[i - 1];
				for (uint k = 1; k < i && k < order; ++k) {
					sum += B[k] * A_i[i - 1 - k];
				}
				A_i[i] = -A_inv * sum;
			}
			Buffer<double> rowSums(activeNodes, 0.0);
			for (uint i = 0; i < activeNodes; ++i) {
				for (uint j = 0; j < order; ++j) {
					for (uint k = 0; i < activeNodes; ++i) {
						rowSums[i] += std::abs(A_i[j][i][k]);
					}
				}
			}
			double weight = edgeLen_X * edgeLen_Y * edgeLen_T;
			for (uint smallstep = 0; smallstep < totalSteps; ++smallstep) {
				MatrixN sum = B[0] * A_i[order - 1];
				for (uint i = 1; i < order; ++i) {
					sum += B[i] * A_i[order - 1 - i];
				}
				MatrixN A_new = -A_inv * sum;
				for (uint i = 1; i < order; ++i) {
					A_i[i - 1] = A_i[i];
				}
				A_i[order - 1] = A_new;
				for (uint i = 0; i < activeNodes; ++i) {
					for (uint j = 0; j < activeNodes; ++j) {
						rowSums[i] += std::abs(A_new[i][j]);
					}
				}
				if ((smallstep + 1) % 10 == 0) { //choose how often to print the result
					double norm = 0.0;
					for (uint i = 0; i < activeNodes; ++i) {
						if (rowSums[i] > norm)
							norm = rowSums[i];
					}
					double weightedNorm = weight * norm;
					std::cout << "Time elapsed: " << dtime + smallstep * dtime / order << " --- System matrix inverse norm: " << weightedNorm << '\n';
					output << "Time elapsed: " << dtime + smallstep * dtime / order << " --- System matrix inverse norm: " << weightedNorm << '\n';
				}
			}

			//compute the norm of the inverse
			double norm = 0.0;
			for (uint i = 0; i < activeNodes; ++i) {
				if (rowSums[i] > norm)
					norm = rowSums[i];
			}
			double weightedNorm = weight * norm;
			std::cout << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output << "Mesh " << meshNumber << " --- System matrix inverse norm in time interval [0, " << T * T_REPEATS << "]: " << weightedNorm << '\n';
			output.save(outputFile.str());
		}
	}

	/*
	Permute the system matrix of waveExampleWhitney2D to a natural order, illustrating the nonzero structure.
	*/
	void waveExampleWhitney2D_PermutationTest(uint order, const uint meshNumber) {
		const double L = 2.0;
		const double T = 2.0;
		const uint Nx = std::pow(2, meshNumber - 1);
		const uint Nt = std::pow(2, meshNumber - 1);

		//Create triangle mesh in [0, L] x [0, T].
		BuilderMesh mesh_old(2);
		createSpacetimeMesh(mesh_old, L, T, Nx, Nt);

		//refine the mesh into small simplices
		BuilderMesh mesh(2);
		SmallSimplexPartition2D ssp(order);
		ssp.refineMesh(mesh_old, mesh);

		//mark boundary nodes of refined mesh
		mesh.fillBoundaryFlags(1);
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1) {
				if (mesh.getNodePosition2(i).y < 1e-14 && mesh.getNodePosition2(i).x > 1e-14 && mesh.getNodePosition2(i).x < L - 1e-14) {
					mesh.setNodeFlag(i, 2);
				}
				else if (mesh.getNodePosition2(i).y > T - 1e-14 && mesh.getNodePosition2(i).x > 1e-14 && mesh.getNodePosition2(i).x < L - 1e-14)
					mesh.setNodeFlag(i, 3);
			}
		}

		//find indices of the nodes in the refined mesh in a natural order
		uint Nx_ref = Nx * order;
		uint Nt_ref = Nt * order;
		uint activeNodes = (Nx_ref - 1) * Nt_ref;
		uint boundary2Nodes = Nx_ref - 1;
		Buffer<uint> indices((Nx_ref - 1) * (Nt_ref + 1));
		double edgeLen_X = L / Nx_ref;
		double edgeLen_T = T / Nt_ref;
		Vector4 pos(0.0, 0.0, 0.0, 0.0);
		uint prev = 0;
		for (uint i = 0; i <= Nt_ref; ++i) {
			for (uint j = 0; j < Nx_ref - 1; ++j) {
				pos.x += edgeLen_X;
				indices[i * (Nx_ref - 1) + j] = prev = mesh.findNode(pos, 1.0e-13, prev);
			}
			pos.x = 0.0;
			pos.y += edgeLen_T;
		}

		Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize()); //form the system matrix
		{
			Buffer<std::unordered_map<uint, double>> star;
			Buffer<uint> refElementIndices(1, 0);
			Buffer <Buffer<uint>> refElements(1);
			refElements[0].resize(mesh_old.getFaceSize());
			for (uint i = 0; i < mesh_old.getFaceSize(); ++i)
				refElements[0][i] = i;
			ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "RightTriangle", false, true, false);
			for (uint j = 0; j < mesh.getNodeSize(); ++j) {
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
					}
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0) {
							system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
						}
					}
				}
			}
		}

		//permute the system matrix to the natural order
		Buffer<Buffer<double>> system_matrix_permuted(activeNodes);
		for (uint i = 0; i < activeNodes; ++i) {
			system_matrix_permuted[i] = Buffer<double>(activeNodes, 0.0);
			std::unordered_map<uint, double>& row = system_matrix[indices[i]];
			for (uint j = 0; j < activeNodes; ++j) {
				auto it = row.find(indices[j + boundary2Nodes]);
				if (it != row.end() && std::abs(it->second) > 1.0e-14) {
					system_matrix_permuted[i][j] = it->second;
				}
			}

		}

		printMatrix(system_matrix_permuted);

		//print the system matrix to file
		Text output;
		output.precision(3);
		output << std::scientific;
		for (uint i = 0; i < system_matrix_permuted.size(); ++i) {
			for (uint j = 0; j < system_matrix_permuted[i].size(); ++j) {
				if (std::abs(system_matrix_permuted[i][j]) > 0.0)
					output << system_matrix_permuted[i][j] << ' ';
				else
					output << '0' << ' ';
			}
			output << '\n';
		}
		Text filepath;
		filepath << "Files/waveExampleWhitney2D_sysmat_order" << order << ".txt";
		output.save(filepath.str());
	}

	/*
	Permute the system matrix of waveExampleWhitney3D to a natural order, illustrating the nonzero structure.
	*/
	void waveExampleWhitney3D_PermutationTest(uint order, const uint meshNumber) {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double T = 1.41;
		const uint Nx = std::pow(2, meshNumber - 1);;
		const uint Ny = std::pow(2, meshNumber - 1);;
		const uint Nt = std::pow(2, meshNumber - 1);;

		//Create tetrahedral mesh in [0, Lx] x [0, Ly] x [0, T].
		BuilderMesh mesh_old(3);
		createSpacetimeMesh(mesh_old, Lx, Ly, T, Nx, Ny, Nt);

		//refine the mesh into small simplices
		BuilderMesh mesh(3);
		SmallSimplexPartition3D ssp(order);
		ssp.refineMesh(mesh_old, mesh);

		//find indices of the nodes in the refined mesh in a natural order
		uint Nx_ref = Nx * order;
		uint Ny_ref = Ny * order;
		uint Nt_ref = Nt * order;
		uint activeNodes = (Nx_ref - 1) * (Ny_ref - 1) * Nt_ref;
		uint boundary2Nodes = (Nx_ref - 1) * (Ny_ref - 1);
		Buffer<uint> indices((Nx_ref - 1) * (Ny_ref - 1) * (Nt_ref + 1));
		double edgeLen_X = Lx / Nx_ref;
		double edgeLen_Y = Ly / Ny_ref;
		double edgeLen_T = T / Nt_ref;
		Vector4 pos(0.0, 0.0, 0.0, 0.0);
		uint prev = 0;
		for (uint i = 0; i <= Nt_ref; ++i) {
			for (uint j = 0; j < Ny_ref - 1; ++j) {
				pos.y += edgeLen_Y;
				for (uint k = 0; k < Nx_ref - 1; ++k) {
					pos.x += edgeLen_X;
					indices[i * (Nx_ref - 1) * (Ny_ref - 1) + j * (Ny_ref - 1) + k] = prev = mesh.findNode(pos, 1.0e-13, prev);
				}
				pos.x = 0.0;
			}
			pos.y = 0.0;
			pos.z += edgeLen_T;
		}

		Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize()); //form the system matrix
			{
				Buffer<std::unordered_map<uint, double>> star;
				Buffer<uint> refElementIndices(6);
				Buffer<Buffer<uint>> refElements(6);
				for (uint i = 0; i < 6; ++i) {
					refElementIndices[i] = i;
					refElements[i].resize(mesh_old.getBodySize() / 6);
				}

				for (uint i = 0; i < mesh_old.getBodySize(); ++i)
					refElements[i % 6][i / 6] = i;
				ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "CubeTetrahedron", false, true, false);
				for (uint j = 0; j < mesh.getNodeSize(); ++j) {
					Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
					const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
					for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
						std::unordered_map<uint, double>& row = star[i];
						for (uint e = 0; e < edges_j.size(); ++e) {
							auto it = row.find(edges_j[e]);
							if (it != row.end())
								col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
						}
					}
					for (uint i = 0; i < mesh.getNodeSize(); ++i) {
						const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
						for (uint e = 0; e < edges_i.size(); ++e) {
							if (col_j_of_star_d[edges_i[e]] != 0.0) {
								system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
							}
						}
					}
				}
			}

		//permute the system matrix to the natural order
		Buffer<Buffer<double>> system_matrix_permuted(activeNodes);
		for (uint i = 0; i < activeNodes; ++i) {
			system_matrix_permuted[i] = Buffer<double>(activeNodes, 0.0);
			std::unordered_map<uint, double>& row = system_matrix[indices[i]];
			for (uint j = 0; j < activeNodes; ++j) {
				auto it = row.find(indices[j + boundary2Nodes]);
				if (it != row.end() && std::abs(it->second) > 1.0e-14) {
					system_matrix_permuted[i][j] = it->second;
				}
			}

		}

		//print the system matrix to file
		Text output;
		output.precision(3);
		output << std::scientific;
		for (uint i = 0; i < system_matrix_permuted.size(); ++i) {
			for (uint j = 0; j < system_matrix_permuted[i].size(); ++j) {
				if (std::abs(system_matrix_permuted[i][j]) > 0.0)
					output << system_matrix_permuted[i][j] << ' ';
				else
					output << '0' << ' ';
			}
			output << '\n';
		}
		Text filepath;
		filepath << "Files/waveExampleWhitney3D_sysmat_order" << order << ".txt";
		output.save(filepath.str());
	}

	/*
	Permute the system matrix of waveExampleCubical3D to a natural order, illustrating the nonzero structure.
	*/
	void waveExampleCubical3D_PermutationTest(uint order, const uint meshNumber) {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const double T = std::sqrt(2);
		const uint Nx = std::pow(2, meshNumber - 1);;
		const uint Ny = std::pow(2, meshNumber - 1);;
		const uint Nt = std::pow(2, meshNumber - 1);;

		//Create cubical mesh in [0, Lx] x [0, Ly] x [0, T].
		BuilderMesh mesh_old(3);
		createCartesianMesh(mesh_old, Lx, Ly, T, Nx, Ny, Nt);

		//refine the mesh into small cubes
		BuilderMesh mesh(3);
		SmallCubePartition3D scp(order);
		scp.refineMesh(mesh_old, mesh);

		//find indices of the nodes in the refined mesh in a natural order
		uint Nx_ref = Nx * order;
		uint Ny_ref = Ny * order;
		uint Nt_ref = Nt * order;
		uint activeNodes = (Nx_ref - 1) * (Ny_ref - 1) * Nt_ref;
		uint boundary2Nodes = (Nx_ref - 1) * (Ny_ref - 1);
		Buffer<uint> indices((Nx_ref - 1) * (Ny_ref - 1) * (Nt_ref + 1));
		double edgeLen_X = Lx / Nx_ref;
		double edgeLen_Y = Ly / Ny_ref;
		double edgeLen_T = T / Nt_ref;
		Vector4 pos(0.0, 0.0, 0.0, 0.0);
		uint prev = 0;
		for (uint i = 0; i <= Nt_ref; ++i) {
			for (uint j = 0; j < Ny_ref - 1; ++j) {
				pos.y += edgeLen_Y;
				for (uint k = 0; k < Nx_ref - 1; ++k) {
					pos.x += edgeLen_X;
					indices[i * (Nx_ref - 1) * (Ny_ref - 1) + j * (Ny_ref - 1) + k] = prev = mesh.findNode(pos, 1.0e-13, prev);
				}
				pos.x = 0.0;
			}
			pos.y = 0.0;
			pos.z += edgeLen_T;
		}

		Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize()); //form the system matrix
		{
			Buffer<std::unordered_map<uint, double>> star;
			scp.formHodgeMatrix1Forms(star, true);
			for (uint j = 0; j < mesh.getNodeSize(); ++j) {
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
					}
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0) {
							system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
						}
					}
				}
			}
		}

		//permute the system matrix to the natural order
		Buffer<Buffer<double>> system_matrix_permuted(activeNodes);
		for (uint i = 0; i < activeNodes; ++i) {
			system_matrix_permuted[i] = Buffer<double>(activeNodes, 0.0);
			std::unordered_map<uint, double>& row = system_matrix[indices[i]];
			for (uint j = 0; j < activeNodes; ++j) {
				auto it = row.find(indices[j + boundary2Nodes]);
				if (it != row.end() && std::abs(it->second) > 1.0e-13) {
					system_matrix_permuted[i][j] = it->second;
				}
			}

		}

		//print the system matrix to file
		Text output;
		output.precision(3);
		output << std::scientific;
		for (uint i = 0; i < system_matrix_permuted.size(); ++i) {
			for (uint j = 0; j < system_matrix_permuted[i].size(); ++j) {
				if (std::abs(system_matrix_permuted[i][j]) > 0.0)
					output << system_matrix_permuted[i][j] << ' ';
				else
					output << '0' << ' ';
			}
			output << '\n';
		}
		Text filepath;
		filepath << "Files/waveExampleCubical3D_sysmat_order" << order << ".txt";
		output.save(filepath.str());
	}

	/*
	Print the blocks in the system matrix of waveExampleCubical2D to files.
	*/
	void waveExampleCubical2D_MatrixBlocks(uint order, const double T, const uint meshNumber) {
		const double L = 2.0;
		const uint spacesteps = std::pow(2, meshNumber - 1);
		const double dtime = 3 * T / spacesteps;
		const double eps = 1e-12;

		//build three time steps of the mesh

		BuilderMesh mesh_old(3);
		createCartesianMesh(mesh_old, L, dtime, spacesteps, 3);
		BuilderMesh mesh(3);
		SmallCubePartition2D scp(order);
		scp.refineMesh(mesh_old, mesh);

		//find indices of the nodes in the refined mesh in a natural order
		uint Nx_ref = spacesteps * order;
		uint Nt_ref = 3 * order;
		uint activeNodes = (Nx_ref - 1) * Nt_ref;
		uint boundary2Nodes = Nx_ref - 1;
		Buffer<uint> indices((Nx_ref - 1) * (Nt_ref + 1));
		double edgeLen_X = L / Nx_ref;
		double edgeLen_T = dtime / Nt_ref;
		Vector4 pos(0.0, 0.0, 0.0, 0.0);
		uint prev = 0;
		for (uint i = 0; i <= Nt_ref; ++i) {
			for (uint j = 0; j < Nx_ref - 1; ++j) {
				pos.x += edgeLen_X;
				indices[i * (Nx_ref - 1) + j] = prev = mesh.findNode(pos, 1.0e-12, prev);
			}
			pos.x = 0.0;
			pos.y += edgeLen_T;
		}

		//form the blocks A, B, and C of the permuted system matrix
		MatrixN blockA, blockB, blockC;
		uint blockSize = boundary2Nodes * order;
		blockA.toMatrixN(blockSize);
		blockB.toMatrixN(blockSize);
		blockC.toMatrixN(blockSize);
		{
			Buffer<std::unordered_map<uint, double>> star;
			scp.formHodgeMatrix1Forms(star, true);
			for (uint j = 0; j < activeNodes; ++j) {
				uint index_j = indices[j + boundary2Nodes];
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(index_j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], index_j);
					}
				}
				for (uint i = 0; i < blockSize; ++i) {
					uint index_i = indices[i + 2 * blockSize];
					const Buffer<uint>& edges_i = mesh.getNodeEdges(index_i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0) {
							if (j < blockSize)
								blockC[i][j] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
							else if (j < 2 * blockSize)
								blockB[i][j - blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
							else
								blockA[i][j - 2 * blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
						}
					}
				}
			}
		}

		//print the blocks A, B, and C into files
		Text filepathA;
		filepathA << "Files/waveExampleCubical2D_order" << order << "_blockA.txt";
		Text filepathB;
		filepathB << "Files/waveExampleCubical2D_order" << order << "_blockB.txt";
		Text filepathC;
		filepathC << "Files/waveExampleCubical2D_order" << order << "_blockC.txt";
		Text outputA;
		outputA.precision(15);
		outputA << std::scientific;
		Text outputB;
		outputB.precision(15);
		outputB << std::scientific;
		Text outputC;
		outputC.precision(15);
		outputC << std::scientific;
		for (uint i = 0; i < blockSize; ++i) {
			for (uint j = 0; j < blockSize; ++j) {
				if (std::abs(blockA[i][j]) > 0)
					outputA << blockA[i][j] << ' ';
				else
					outputA << '0' << ' ';
			}
			outputA << '\n';
		}
		outputA.save(filepathA.str());
		for (uint i = 0; i < blockSize; ++i) {
			for (uint j = 0; j < blockSize; ++j) {
				if (std::abs(blockB[i][j]) > 0)
					outputB << blockB[i][j] << ' ';
				else
					outputB << '0' << ' ';
			}
			outputB << '\n';
		}
		outputB.save(filepathB.str());
		for (uint i = 0; i < blockSize; ++i) {
			for (uint j = 0; j < blockSize; ++j) {
				if (std::abs(blockC[i][j]) > 0)
					outputC << blockC[i][j] << ' ';
				else
					outputC << '0' << ' ';
			}
			outputC << '\n';
		}
		outputC.save(filepathC.str());
	}

	/*
	Print the blocks in the system matrix of waveExampleCubical3D to files.
	*/
	void waveExampleCubical3D_MatrixBlocks(uint order, const double T, const uint meshNumber) {
		const double Lx = 2.0;
		const double Ly = 2.0;
		const uint spacesteps = std::pow(2, meshNumber - 1);
		const double dtime = 3 * T / spacesteps;
		const double eps = 1e-12;

		//build three time steps of the mesh

		BuilderMesh mesh_old(3);
		createCartesianMesh(mesh_old, Lx, Ly, dtime, spacesteps, spacesteps, 3);
		BuilderMesh mesh(3);
		SmallCubePartition3D scp(order);
		scp.refineMesh(mesh_old, mesh);

		//find indices of the nodes in the refined mesh in a natural order
		uint Nx_ref = spacesteps * order;
		uint Ny_ref = spacesteps * order;
		uint Nt_ref = 3 * order;
		uint activeNodes = (Nx_ref - 1) * (Ny_ref - 1) * Nt_ref;
		uint boundary2Nodes = (Nx_ref - 1) * (Ny_ref - 1);
		Buffer<uint> indices((Nx_ref - 1) * (Ny_ref - 1) * (Nt_ref + 1));
		double edgeLen_X = Lx / Nx_ref;
		double edgeLen_Y = Ly / Ny_ref;
		double edgeLen_T = dtime / Nt_ref;
		Vector4 pos(0.0, 0.0, 0.0, 0.0);
		uint prev = 0;
		for (uint i = 0; i <= Nt_ref; ++i) {
			for (uint j = 0; j < Ny_ref - 1; ++j) {
				pos.y += edgeLen_Y;
				for (uint k = 0; k < Nx_ref - 1; ++k) {
					pos.x += edgeLen_X;
					indices[i * (Nx_ref - 1) * (Ny_ref - 1) + j * (Ny_ref - 1) + k] = prev = mesh.findNode(pos, 1.0e-12, prev);
				}
				pos.x = 0.0;
			}
			pos.y = 0.0;
			pos.z += edgeLen_T;
		}

		//form the blocks A, B, and C of the permuted system matrix
		MatrixN blockA, blockB, blockC;
		uint blockSize = boundary2Nodes * order;
		blockA.toMatrixN(blockSize);
		blockB.toMatrixN(blockSize);
		blockC.toMatrixN(blockSize);
		{
			Buffer<std::unordered_map<uint, double>> star;
			scp.formHodgeMatrix1Forms(star, true);
			for (uint j = 0; j < activeNodes; ++j) {
				uint index_j = indices[j + boundary2Nodes];
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(index_j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], index_j);
					}
				}
				for (uint i = 0; i < blockSize; ++i) {
					uint index_i = indices[i + 2 * blockSize];
					const Buffer<uint>& edges_i = mesh.getNodeEdges(index_i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0) {
							if (j < blockSize)
								blockC[i][j] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
							else if (j < 2 * blockSize)
								blockB[i][j - blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
							else
								blockA[i][j - 2 * blockSize] -= mesh.getEdgeIncidence(edges_i[e], index_i) * col_j_of_star_d[edges_i[e]];
						}
					}
				}
			}
		}

		//print the blocks A, B, and C into files
		Text filepathA;
		filepathA << "Files/waveExampleCubical3D_order" << order << "_blockA.txt";
		Text filepathB;
		filepathB << "Files/waveExampleCubical3D_order" << order << "_blockB.txt";
		Text filepathC;
		filepathC << "Files/waveExampleCubical3D_order" << order << "_blockC.txt";
		Text outputA;
		outputA.precision(15);
		outputA << std::scientific;
		Text outputB;
		outputB.precision(15);
		outputB << std::scientific;
		Text outputC;
		outputC.precision(15);
		outputC << std::scientific;
		for (uint i = 0; i < blockSize; ++i) {
			for (uint j = 0; j < blockSize; ++j) {
				if (std::abs(blockA[i][j]) > 0)
					outputA << blockA[i][j] << ' ';
				else
					outputA << '0' << ' ';
			}
			outputA << '\n';
		}
		outputA.save(filepathA.str());
		for (uint i = 0; i < blockSize; ++i) {
			for (uint j = 0; j < blockSize; ++j) {
				if (std::abs(blockB[i][j]) > 0)
					outputB << blockB[i][j] << ' ';
				else
					outputB << '0' << ' ';
			}
			outputB << '\n';
		}
		outputB.save(filepathB.str());
		for (uint i = 0; i < blockSize; ++i) {
			for (uint j = 0; j < blockSize; ++j) {
				if (std::abs(blockC[i][j]) > 0)
					outputC << blockC[i][j] << ' ';
				else
					outputC << '0' << ' ';
			}
			outputC << '\n';
		}
		outputC.save(filepathC.str());
	}

	/*
	Print the blocks in the system matrix of waveExampleSmallSteps3D to files.
	*/
	void waveExampleSmallSteps3D_MatrixBlocks(uint order, const double T, const uint meshNumber) {
		const double Lx = 2;
		const double Ly = 2;
		const uint spacesteps = std::pow(2, meshNumber - 1);
		const double dtime = T / spacesteps;
		const double eps = 1e-12;
		const double smallstepsize = T / order;
		BuilderMesh mesh_old(3);
		createCartesianMesh(mesh_old, Lx, Ly, dtime, spacesteps, spacesteps, 1);
		BuilderMesh mesh(3);
		SmallCubePartition3D scp(order);
		scp.refineMesh(mesh_old, mesh);
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			Vector3 p = mesh.getNodePosition3(i);
			if (p.x < eps || p.x > Lx - eps || p.y < eps || p.y > Ly - eps)
				mesh.setNodeFlag(i, 1);
			else if (p.z < eps)
				mesh.setNodeFlag(i, 2);
			else if (p.z > dtime - eps)
				mesh.setNodeFlag(i, 3);
			else
				mesh.setNodeFlag(i, 0);
		}

		//find indices of the nodes in the refined mesh in a natural order
		uint Nx_ref = spacesteps * order;
		uint Ny_ref = spacesteps * order;
		uint Nt_ref = order;
		uint timeStepNodes = (Nx_ref + 1) * (Ny_ref + 1);
		uint activeNodes = (Nx_ref - 1) * (Ny_ref - 1);
		Buffer<Buffer<uint>> indices(order + 1);
		for (uint i = 0; i < indices.size(); ++i)
			indices[i].resize(timeStepNodes);
		double edgeLen_X = Lx / Nx_ref;
		double edgeLen_Y = Ly / Ny_ref;
		double edgeLen_T = dtime / Nt_ref;
		Vector4 pos(0.0, 0.0, 0.0, 0.0);
		uint prev = 0;
		for (uint i = 0; i <= order; ++i) {
			for (uint j = 0; j < Ny_ref + 1; ++j) {
				for (uint k = 0; k < Nx_ref + 1; ++k) {
					indices[i][j * (Ny_ref + 1) + k] = prev = mesh.findNode(pos, 1.0e-13, prev);
					pos.x += edgeLen_X;
				}
				pos.x = 0.0;
				pos.y += edgeLen_Y;
			}
			pos.y = 0.0;
			pos.z += edgeLen_T;
		}

		Buffer<Buffer<double>> system_matrix(mesh.getNodeSize()); //form the system matrix
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			system_matrix[i] = Buffer<double>(mesh.getNodeSize(), 0.0);
		}
		{
			Buffer<std::unordered_map<uint, double>> star;
			scp.formHodgeMatrix1Forms(star, true);
			for (uint j = 0; j < mesh.getNodeSize(); ++j) {
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
					}
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0) {
							system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
						}
					}
				}
			}
		}
		Buffer<uint> activeNodeIndices(activeNodes);
		for (uint i = 0, j = 0; i < timeStepNodes; ++i) {
			if (mesh.getNodeFlag(indices[order][i]) != 1)
				activeNodeIndices[j++] = indices[order][i];
		}
		Buffer<uint> activeDualIndices(activeNodes);
		for (uint i = 0, j = 0; i < timeStepNodes; ++i) {
			if (mesh.getNodeFlag(indices[order - 1][i]) != 1)
				activeDualIndices[j++] = indices[order - 1][i];
		}
		Buffer<Buffer<Buffer<double>>> updateValues(order);
		for (uint i = 0; i < order; ++i) {
			updateValues[i].resize(activeNodes);
			for (uint j = 0; j < activeNodes; ++j) {
				updateValues[i][j].resize(activeNodes);
				uint index_j = activeDualIndices[j];
				for (uint k = 0, l = 0; k < timeStepNodes; ++k) {
					if (mesh.getNodeFlag(indices[i][k]) == 1)
						continue;
					updateValues[i][j][l++] = system_matrix[index_j][indices[i][k]];
				}
			}
		}
		Buffer<Buffer<double>> system_matrix_restricted(activeNodes);
		for (uint i = 0; i < activeDualIndices.size(); ++i) {
			system_matrix_restricted[i].resize(activeNodes);
			uint index_i = activeDualIndices[i];
			for (uint j = 0; j < activeNodes; ++j) {
				uint index_j = activeNodeIndices[j];
				system_matrix_restricted[i][j] = system_matrix[index_i][index_j];

			}
		}

		Buffer<MatrixN> B(order + 1);
		B[0] = MatrixN(system_matrix_restricted);
		for (uint i = 0; i < order; ++i)
			B[i + 1] = MatrixN(updateValues[order - i - 1]);

		//print the matrices to file
		for (uint k = 0; k <= order; ++k) {
			Text filepath;
			filepath << "Files/waveExampleSmallSteps3D_order" << order << "_block" << k << ".txt";
			Text output;
			output.precision(15);
			output << std::scientific;
			for (uint i = 0; i < B[k].size(); ++i) {
				for (uint j = 0; j < B[k][i].size(); ++j) {
					if (std::abs(B[k][i][j]) > 0)
						output << B[k][i][j] << ' ';
					else
						output << '0' << ' ';
				}
				output << '\n';
			}
			output.save(filepath.str());
		}
	}

	/*
	Compute error L2 norm when approximating fn with the interpolant of discreteForm. If mesh_old has at least minElements elements, we use it to compute the error.
	Otherwise, mesh_old is refined and the error is computed using a refined mesh with at least minElements elements.
	*/
	double computeDifferenceL2Norm0Form(const Mesh& mesh_old, const SmallSimplexPartition2D& ssp, const Buffer<double>& discreteForm, const std::function<double(Vector2)>& fn, uint minElements) {
		Buffer<double> nodeCoefficients;
		Buffer<VectorN> edgeCoefficients;
		Buffer<VectorN> faceCoefficients;
		ssp.solve0FormCoefficients(discreteForm, nodeCoefficients, edgeCoefficients, faceCoefficients);
		double l2NormSquared = 0.0;
		if (minElements <= mesh_old.getFaceSize()) {
			for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
				const Buffer<uint> nodes = mesh_old.getFaceNodes(i);
				double A = std::abs(mesh_old.getFaceVector2(i).determinant());
				l2NormSquared += A * integralAverage([&ssp, &nodeCoefficients, &edgeCoefficients, &faceCoefficients, &fn, i](Vector2 p) -> double {
					return std::pow(ssp.evaluate0FormWithCoefficients(p, nodeCoefficients, edgeCoefficients, faceCoefficients, i) - fn(p), 2);
				}, mesh_old.getNodePosition2(nodes[0]), mesh_old.getNodePosition2(nodes[1]), mesh_old.getNodePosition2(nodes[2]));
			}
		}
		else {
			BuilderMesh refinement(2);
			refinement.createCopy(mesh_old);
			while (refinement.getFaceSize() < minElements) {
				BuilderMesh helpmesh(2);
				SmallSimplexPartition2D refiner(2);
				refiner.refineMesh(refinement, helpmesh);
				refinement.createCopy(helpmesh);
			}
			uint partsPerElement = refinement.getFaceSize() / mesh_old.getFaceSize();
			for (uint i = 0; i < refinement.getFaceSize(); ++i) {
				uint i_big = i / partsPerElement;
				const Buffer<uint> nodes = refinement.getFaceNodes(i);
				double A = std::abs(refinement.getFaceVector2(i).determinant());
				l2NormSquared += A * integralAverage([&ssp, &nodeCoefficients, &edgeCoefficients, &faceCoefficients, &fn, i_big](Vector2 p) -> double {
					return std::pow(ssp.evaluate0FormWithCoefficients(p, nodeCoefficients, edgeCoefficients, faceCoefficients, i_big) - fn(p), 2);
				}, refinement.getNodePosition2(nodes[0]), refinement.getNodePosition2(nodes[1]), refinement.getNodePosition2(nodes[2]));
			}
		}
		return std::sqrt(l2NormSquared);
	}

	/*
	Compute error L2 norm when approximating fn with the interpolant of discreteForm. If mesh_old has at least minElements elements, we use it to compute the error.
	Otherwise, mesh_old is refined and the error is computed using a refined mesh with at least minElements elements.
	*/
	double computeDifferenceL2Norm0Form(const Mesh& mesh_old, const SmallSimplexPartition3D& ssp, const Buffer<double>& discreteForm, const std::function<double(Vector3)>& fn, uint minElements) {
		Buffer<double> nodeCoefficients;
		Buffer<VectorN> edgeCoefficients;
		Buffer<VectorN> faceCoefficients;
		Buffer<VectorN> bodyCoefficients;
		ssp.solve0FormCoefficients(discreteForm, nodeCoefficients, edgeCoefficients, faceCoefficients, bodyCoefficients);
		double l2NormSquared = 0.0;
		if (minElements <= mesh_old.getBodySize()) {
			for (uint i = 0; i < mesh_old.getBodySize(); ++i) {
				const Buffer<uint> nodes = mesh_old.getBodyNodes(i);
				double V = std::abs(mesh_old.getBodyVector3(i).determinant());
				l2NormSquared += V * integralAverage([&ssp, &nodeCoefficients, &edgeCoefficients, &faceCoefficients, &bodyCoefficients, &fn, i](Vector3 p) -> double {
					return std::pow(ssp.evaluate0FormWithCoefficients(p, nodeCoefficients, edgeCoefficients, faceCoefficients, bodyCoefficients, i) - fn(p), 2);
				}, mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[1]), mesh_old.getNodePosition3(nodes[2]), mesh_old.getNodePosition3(nodes[3]));
			}
		}
		else {
			BuilderMesh refinement(3);
			refinement.createCopy(mesh_old);
			while (refinement.getBodySize() < minElements) {
				BuilderMesh helpmesh(3);
				SmallSimplexPartition3D refiner(2);
				refiner.refineMesh(refinement, helpmesh);
				refinement.createCopy(helpmesh);
			}
			uint partsPerElement = refinement.getBodySize() / mesh_old.getBodySize();
			for (uint i = 0; i < refinement.getBodySize(); ++i) {
				uint i_big = i / partsPerElement;
				const Buffer<uint> nodes = refinement.getBodyNodes(i);
				double V = std::abs(refinement.getBodyVector3(i).determinant());
				l2NormSquared += V * integralAverage([&ssp, &nodeCoefficients, &edgeCoefficients, &faceCoefficients, &bodyCoefficients, &fn, i_big](Vector3 p) -> double {
					return std::pow(ssp.evaluate0FormWithCoefficients(p, nodeCoefficients, edgeCoefficients, faceCoefficients, bodyCoefficients, i_big) - fn(p), 2);
				}, refinement.getNodePosition3(nodes[0]), refinement.getNodePosition3(nodes[1]), refinement.getNodePosition3(nodes[2]), refinement.getNodePosition3(nodes[3]));
			}
		}
		return std::sqrt(l2NormSquared);
	}

	/*
	Compute error L2 norm when approximating fn with the interpolant of discreteForm. If mesh_old has at least minElements elements, we use it to compute the error.
	Otherwise, mesh_old is refined and the error is computed using a refined mesh with at least minElements elements.
	*/
	double computeDifferenceL2Norm1Form(const Mesh& mesh_old, const SmallSimplexPartition2D& ssp, const Buffer<double>& discreteForm, const std::function<Vector2(Vector2)>& fn, uint minElements) {
		Buffer<VectorN> edgeCoefficients;
		Buffer<VectorN> faceCoefficients;
		ssp.solve1FormCoefficients(discreteForm, edgeCoefficients, faceCoefficients);
		double l2NormSquared = 0.0;
		if (minElements <= mesh_old.getFaceSize()) {
			for (uint i = 0; i < mesh_old.getFaceSize(); i++) {
				const Buffer<uint> nodes = mesh_old.getFaceNodes(i);
				double A = std::abs(mesh_old.getFaceVector2(i).determinant());
				l2NormSquared += A * integralAverage([&ssp, &edgeCoefficients, &faceCoefficients, &fn, i](Vector2 p) -> double {
					return (ssp.evaluate1FormWithCoefficients(p, edgeCoefficients, faceCoefficients, i) - fn(p)).lensq();
				}, mesh_old.getNodePosition2(nodes[0]), mesh_old.getNodePosition2(nodes[1]), mesh_old.getNodePosition2(nodes[2]));
			}
		}
		else {
			BuilderMesh refinement(2);
			refinement.createCopy(mesh_old);
			while (refinement.getFaceSize() < minElements) {
				BuilderMesh helpmesh(2);
				SmallSimplexPartition2D refiner(2);
				refiner.refineMesh(refinement, helpmesh);
				refinement.createCopy(helpmesh);
			}
			uint partsPerElement = refinement.getFaceSize() / mesh_old.getFaceSize();
			for (uint i = 0; i < refinement.getFaceSize(); ++i) {
				uint i_big = i / partsPerElement;
				const Buffer<uint> nodes = refinement.getFaceNodes(i);
				double A = std::abs(refinement.getFaceVector2(i).determinant());
				l2NormSquared += A * integralAverage([&ssp, &edgeCoefficients, &faceCoefficients, &fn, i_big](Vector2 p) -> double {
					return (ssp.evaluate1FormWithCoefficients(p, edgeCoefficients, faceCoefficients, i_big) - fn(p)).lensq();
				}, refinement.getNodePosition2(nodes[0]), refinement.getNodePosition2(nodes[1]), refinement.getNodePosition2(nodes[2]));
			}
		}
		return std::sqrt(l2NormSquared);
	}

	/*
	Compute error L2 norm when approximating fn with the interpolant of discreteForm. If mesh_old has at least minElements elements, we use it to compute the error.
	Otherwise, mesh_old is refined and the error is computed using a refined mesh with at least minElements elements.
	*/
	double computeDifferenceL2Norm1Form(const Mesh& mesh_old, const SmallSimplexPartition3D& ssp, const Buffer<double>& discreteForm, const std::function<Vector3(Vector3)>& fn, uint minElements) {
		Buffer<VectorN> edgeCoefficients;
		Buffer<VectorN> faceCoefficients;
		Buffer<VectorN> bodyCoefficients;
		ssp.solve1FormCoefficients(discreteForm, edgeCoefficients, faceCoefficients, bodyCoefficients);
		double l2NormSquared = 0.0;
		if (minElements <= mesh_old.getBodySize()) {
			for (uint i = 0; i < mesh_old.getBodySize(); i++) {
				const Buffer<uint> nodes = mesh_old.getBodyNodes(i);
				double V = std::abs(mesh_old.getBodyVector3(i).determinant());
				l2NormSquared += V * integralAverage([&ssp, &edgeCoefficients, &faceCoefficients, &bodyCoefficients, &fn, i](Vector3 p) -> double {
					return (ssp.evaluate1FormWithCoefficients(p, edgeCoefficients, faceCoefficients, bodyCoefficients, i) - fn(p)).lensq();
				}, mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[1]), mesh_old.getNodePosition3(nodes[2]), mesh_old.getNodePosition3(nodes[3]));
			}
		}
		else {
			BuilderMesh refinement(3);
			refinement.createCopy(mesh_old);
			while (refinement.getBodySize() < minElements) {
				BuilderMesh helpmesh(3);
				SmallSimplexPartition3D refiner(2);
				refiner.refineMesh(refinement, helpmesh);
				refinement.createCopy(helpmesh);
			}
			uint partsPerElement = refinement.getBodySize() / mesh_old.getBodySize();
			for (uint i = 0; i < refinement.getBodySize(); ++i) {
				uint i_big = i / partsPerElement;
				const Buffer<uint> nodes = refinement.getBodyNodes(i);
				double V = std::abs(refinement.getBodyVector3(i).determinant());
				l2NormSquared += V * integralAverage([&ssp, &edgeCoefficients, &faceCoefficients, &bodyCoefficients, &fn, i_big](Vector3 p) -> double {
					return (ssp.evaluate1FormWithCoefficients(p, edgeCoefficients, faceCoefficients, bodyCoefficients, i_big) - fn(p)).lensq();
				}, refinement.getNodePosition3(nodes[0]), refinement.getNodePosition3(nodes[1]), refinement.getNodePosition3(nodes[2]), refinement.getNodePosition3(nodes[3]));
			}
		}
		return std::sqrt(l2NormSquared);
	}

	/*
	Compute error L2 norm when approximating fn with the interpolant of discreteForm. If mesh_old has at least minElements elements, we use it to compute the error.
	Otherwise, mesh_old is refined and the error is computed using a refined mesh with at least minElements elements.
	*/
	double computeDifferenceL2Norm2Form(const Mesh& mesh_old, const SmallSimplexPartition2D& ssp, const Buffer<double>& discreteForm, const std::function<double(Vector2)>& fn, uint minElements) {
		Buffer<VectorN> faceCoefficients;
		ssp.solve2FormCoefficients(discreteForm, faceCoefficients);
		double l2NormSquared = 0.0;
		if (minElements <= mesh_old.getFaceSize()) {
			for (uint i = 0; i < mesh_old.getFaceSize(); i++) {
				const Buffer<uint> nodes = mesh_old.getFaceNodes(i);
				double A = std::abs(mesh_old.getFaceVector2(i).determinant());
				l2NormSquared += A * integralAverage([&ssp, &faceCoefficients, &fn, i](Vector2 p) -> double {
					return std::pow(ssp.evaluate2FormWithCoefficients(p, faceCoefficients, i) - fn(p), 2);
				}, mesh_old.getNodePosition2(nodes[0]), mesh_old.getNodePosition2(nodes[1]), mesh_old.getNodePosition2(nodes[2]));
			}
		}
		else {
			BuilderMesh refinement(2);
			refinement.createCopy(mesh_old);
			while (refinement.getFaceSize() < minElements) {
				BuilderMesh helpmesh(2);
				SmallSimplexPartition2D refiner(2);
				refiner.refineMesh(refinement, helpmesh);
				refinement.createCopy(helpmesh);
			}
			uint partsPerElement = refinement.getFaceSize() / mesh_old.getFaceSize();
			for (uint i = 0; i < refinement.getFaceSize(); ++i) {
				uint i_big = i / partsPerElement;
				const Buffer<uint> nodes = refinement.getFaceNodes(i);
				double A = std::abs(refinement.getFaceVector2(i).determinant());
				l2NormSquared += A * integralAverage([&ssp, &faceCoefficients, &fn, i_big](Vector2 p) -> double {
					return std::pow(ssp.evaluate2FormWithCoefficients(p, faceCoefficients, i_big) - fn(p), 2);
				}, refinement.getNodePosition2(nodes[0]), refinement.getNodePosition2(nodes[1]), refinement.getNodePosition2(nodes[2]));
			}
		}
		return std::sqrt(l2NormSquared);
	}

	/*
	Compute error L2 norm when approximating fn with the interpolant of discreteForm. If mesh_old has at least minElements elements, we use it to compute the error.
	Otherwise, mesh_old is refined and the error is computed using a refined mesh with at least minElements elements.
	*/
	double computeDifferenceL2Norm2Form(const Mesh& mesh_old, const SmallSimplexPartition3D& ssp, const Buffer<double>& discreteForm, const std::function<Vector3(Vector3)>& fn, uint minElements) {
		Buffer<VectorN> faceCoefficients;
		Buffer<VectorN> bodyCoefficients;
		ssp.solve2FormCoefficients(discreteForm, faceCoefficients, bodyCoefficients);
		double l2NormSquared = 0.0;
		if (minElements <= mesh_old.getBodySize()) {
			for (uint i = 0; i < mesh_old.getBodySize(); i++) {
				const Buffer<uint> nodes = mesh_old.getBodyNodes(i);
				double V = std::abs(mesh_old.getBodyVector3(i).determinant());
				l2NormSquared += V * integralAverage([&ssp, &faceCoefficients, &bodyCoefficients, &fn, i](Vector3 p) -> double {
					return (ssp.evaluate2FormWithCoefficients(p, faceCoefficients, bodyCoefficients, i) - fn(p)).lensq();
				}, mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[1]), mesh_old.getNodePosition3(nodes[2]), mesh_old.getNodePosition3(nodes[3]));
			}
		}
		else {
			BuilderMesh refinement(3);
			refinement.createCopy(mesh_old);
			while (refinement.getBodySize() < minElements) {
				BuilderMesh helpmesh(3);
				SmallSimplexPartition3D refiner(2);
				refiner.refineMesh(refinement, helpmesh);
				refinement.createCopy(helpmesh);
			}
			uint partsPerElement = refinement.getBodySize() / mesh_old.getBodySize();
			for (uint i = 0; i < refinement.getBodySize(); ++i) {
				uint i_big = i / partsPerElement;
				const Buffer<uint> nodes = refinement.getBodyNodes(i);
				double V = std::abs(refinement.getBodyVector3(i).determinant());
				l2NormSquared += V * integralAverage([&ssp, &faceCoefficients, &bodyCoefficients, &fn, i_big](Vector3 p) -> double {
					return (ssp.evaluate2FormWithCoefficients(p, faceCoefficients, bodyCoefficients, i_big) - fn(p)).lensq();
				}, refinement.getNodePosition3(nodes[0]), refinement.getNodePosition3(nodes[1]), refinement.getNodePosition3(nodes[2]), refinement.getNodePosition3(nodes[3]));
			}
		}
		return std::sqrt(l2NormSquared);
	}

	/*
	Compute error L2 norm when approximating fn with the interpolant of discreteForm. If mesh_old has at least minElements elements, we use it to compute the error.
	Otherwise, mesh_old is refined and the error is computed using a refined mesh with at least minElements elements.
	*/
	double computeDifferenceL2Norm3Form(const Mesh& mesh_old, const SmallSimplexPartition3D& ssp, const Buffer<double>& discreteForm, const std::function<double(Vector3)>& fn, uint minElements) {
		Buffer<VectorN> bodyCoefficients;
		ssp.solve3FormCoefficients(discreteForm, bodyCoefficients);
		double l2NormSquared = 0.0;
		if (minElements <= mesh_old.getBodySize()) {
			for (uint i = 0; i < mesh_old.getBodySize(); i++) {
				const Buffer<uint> nodes = mesh_old.getBodyNodes(i);
				double V = std::abs(mesh_old.getBodyVector3(i).determinant());
				l2NormSquared += V * integralAverage([&ssp, &bodyCoefficients, &fn, i](Vector3 p) -> double {
					return std::pow(ssp.evaluate3FormWithCoefficients(p, bodyCoefficients, i) - fn(p), 2);
				}, mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[1]), mesh_old.getNodePosition3(nodes[2]), mesh_old.getNodePosition3(nodes[3]));
			}
		}
		else {
			BuilderMesh refinement(3);
			refinement.createCopy(mesh_old);
			while (refinement.getBodySize() < minElements) {
				BuilderMesh helpmesh(3);
				SmallSimplexPartition3D refiner(2);
				refiner.refineMesh(refinement, helpmesh);
				refinement.createCopy(helpmesh);
			}
			uint partsPerElement = refinement.getBodySize() / mesh_old.getBodySize();
			for (uint i = 0; i < refinement.getBodySize(); ++i) {
				uint i_big = i / partsPerElement;
				const Buffer<uint> nodes = refinement.getBodyNodes(i);
				double V = std::abs(refinement.getBodyVector3(i).determinant());
				l2NormSquared += V * integralAverage([&ssp, &bodyCoefficients, &fn, i_big](Vector3 p) -> double {
					return std::pow(ssp.evaluate3FormWithCoefficients(p, bodyCoefficients, i_big) - fn(p), 2);
				}, refinement.getNodePosition3(nodes[0]), refinement.getNodePosition3(nodes[1]), refinement.getNodePosition3(nodes[2]), refinement.getNodePosition3(nodes[3]));
			}
		}
		return std::sqrt(l2NormSquared);
	}

	//Compute error L2 norm when approximating fn with the interpolant of discreteForm.
	double computeDifferenceL2Norm1Form(const Mesh& mesh_old, const SmallCubePartition2D& scp, const Buffer<double>& discreteForm, const std::function<Vector2(Vector2)>& fn,
		double maxLen) {
		Buffer<VectorN> edgeCoefs;
		Buffer<VectorN> faceCoefs;
		scp.solve1FormCoefficients(discreteForm, edgeCoefs, faceCoefs);
		double l2NormSquared = 0.0;
		for (uint i = 0; i < mesh_old.getFaceSize(); ++i) {
			const Buffer<uint> nodes = getQuadrilateralNodes(i, mesh_old);
			l2NormSquared += parallelogramIntegral([&scp, &edgeCoefs, &faceCoefs, &fn, i](Vector2 p) -> double {
				return (scp.evaluate1FormWithCoefficients(p, edgeCoefs, faceCoefs, i) - fn(p)).lensq();
			}, mesh_old.getNodePosition2(nodes[0]), mesh_old.getNodePosition2(nodes[1]), mesh_old.getNodePosition2(nodes[2]), mesh_old.getNodePosition2(nodes[3]), maxLen);
		}
		return std::sqrt(l2NormSquared);
	}

	//Compute error L2 norm when approximating fn with the interpolant of discreteForm.
	double computeDifferenceL2Norm1Form(const Mesh& mesh_old, const SmallCubePartition3D& scp, const Buffer<double>& discreteForm, const std::function<Vector3(Vector3)>& fn,
		double maxLen) {
		Buffer<VectorN> edgeCoefs;
		Buffer<VectorN> faceCoefs;
		Buffer<VectorN> bodyCoefs;
		scp.solve1FormCoefficients(discreteForm, edgeCoefs, faceCoefs, bodyCoefs);
		double l2NormSquared = 0.0;
		for (uint i = 0; i < mesh_old.getBodySize(); ++i) {
			const Buffer<uint> nodes = getCubeNodes(i, mesh_old);
			l2NormSquared += parallelepipedIntegral([&scp, &edgeCoefs, &faceCoefs, &bodyCoefs, &fn, i](Vector3 p) -> double {
				return (scp.evaluate1FormWithCoefficients(p, edgeCoefs, faceCoefs, bodyCoefs, i) - fn(p)).lensq();
			}, mesh_old.getNodePosition3(nodes[0]), mesh_old.getNodePosition3(nodes[1]), mesh_old.getNodePosition3(nodes[2]), mesh_old.getNodePosition3(nodes[3]),
				mesh_old.getNodePosition3(nodes[4]), mesh_old.getNodePosition3(nodes[5]), mesh_old.getNodePosition3(nodes[6]), mesh_old.getNodePosition3(nodes[7]), maxLen);
		}
		return std::sqrt(l2NormSquared);
	}

	/*
	Get the solution to the linear system using the LU decomposition.
	Matrix mat (of size (nodes x nodes)) is restricted to rows correponding to those nodes for which mesh.getNodeFlag(i) != 1 and mesh.getNodeFlag(i) != rowFlag and
	columns correponding to those nodes for which mesh.getNodeFlag(i) != 1 and mesh.getNodeFlag(i) != colFlag.
	The right-hand side b should have size equal to the number of such nodes.
	*/
	void solveSystem(const Buffer<std::unordered_map<uint, double>>& mat, Buffer<double>& x, const Buffer<double>& b, const Mesh& mesh, int rowFlag, int colFlag) {
		Buffer<Buffer<double>> system_matrix_interior(b.size());
		for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == rowFlag)
				continue;
			system_matrix_interior[ii].resize(b.size());
			for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
				if (mesh.getNodeFlag(j) == 1 || mesh.getNodeFlag(j) == colFlag)
					continue;
				system_matrix_interior[ii][jj] = mat[i][j];
				++jj;
			}
			++ii;
		}
		Buffer<uint> p;
		decomposeLUP(system_matrix_interior, p, 1e-15);
		solveLUP(system_matrix_interior, p, b, x);
	}

	//Computes the solution using the refined mesh in polynomialPoissonExampleWhitney2D and poissonExampleWhitney2D.
	void getDECSolutionPoissonEq(const BuilderMesh& mesh, const SmallSimplexPartition2D& ssp, const std::function<double(Vector2)>& sourceTerm, const std::function<double(Vector2)>& dirichletBC,
		const std::function<Vector2(Vector2)>& neumannBC, Buffer<double>& sol, bool circumcentric, bool oneElement) {
		Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize()); //form the system matrix
		{
			Buffer<std::unordered_map<uint, double>> star;
			if (oneElement) //if all elements have same shape, use precomputed integrals
				ssp.formHodgeMatrix1FormsCustom(star, "DefaultTriangle", circumcentric, false);
			else
				ssp.formHodgeMatrix1Forms(star, circumcentric);
			for (uint j = 0; j < mesh.getNodeSize(); ++j) {
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
					}
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0)
							system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
					}
				}
			}
		}
		//count the number of active nodes
		uint activeNodes = 0;
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) != 1)
				++activeNodes;
		}
		//form the right-hand side
		Buffer<double> rhs(activeNodes, 0.0);
		Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1)
				boundaryValues[i] = dirichletBC(mesh.getNodePosition2(i));
		}
		for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1)
				continue;
			rhs[j] = integrateDual2Cell(mesh, i, sourceTerm, circumcentric);
			std::unordered_map<uint, double>& row = system_matrix[i];
			for (auto it = row.begin(); it != row.end(); ++it) {
				rhs[j] -= it->second * boundaryValues[it->first];
			}
			if (mesh.getNodeFlag(i) == 2)
				rhs[j] -= integrateNodeDualBoundary(mesh, i, neumannBC, circumcentric);
			++j;
		}
		//get the solution on active nodes by solving the system and on other nodes from the boundary condition
		Buffer<double> sol_interior;
		Buffer<Buffer<double>> system_matrix_interior(activeNodes);
		for (uint i = 0, ii = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1)
				continue;
			system_matrix_interior[ii].resize(activeNodes);
			for (uint j = 0, jj = 0; j < mesh.getNodeSize(); ++j) {
				if (mesh.getNodeFlag(j) == 1)
					continue;
				system_matrix_interior[ii][jj] = system_matrix[i][j];
				++jj;

			}
			++ii;
		}
		Buffer<uint> p;
		decomposeLUP(system_matrix_interior, p, 1e-15);
		solveLUP(system_matrix_interior, p, rhs, sol_interior);
		sol.resize(mesh.getNodeSize());
		for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1)
				sol[i] = dirichletBC(mesh.getNodePosition2(i));
			else
				sol[i] = sol_interior[j++];
		}
	}

	//Computes the solution using the refined mesh in polynomialPoissonExampleWhitney3D and poissonExampleWhitney3D.
	void getDECSolutionPoissonEq(const BuilderMesh& mesh, const SmallSimplexPartition3D& ssp, const std::function<double(Vector3)>& sourceTerm, const std::function<double(Vector3)>& dirichletBC,
		const std::function<Vector3(Vector3)>& neumannBC, Buffer<double>& sol, bool circumcentric, bool oneElement) {
		Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize()); //form the system matrix
		{
			Buffer<std::unordered_map<uint, double>> star;
			if (oneElement) //if all elements have same shape, use precomputed integrals
				ssp.formHodgeMatrix1FormsCustom(star, "BccTetrahedron", circumcentric, false);
			else
				ssp.formHodgeMatrix1Forms(star, circumcentric);
			for (uint j = 0; j < mesh.getNodeSize(); ++j) {
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
					}
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0)
							system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
					}
				}
			}
		}
		//count the number of active nodes
		uint activeNodes = 0;
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) != 1)
				++activeNodes;
		}
		//form the right-hand side
		Buffer<double> rhs(activeNodes, 0.0);
		Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1)
				boundaryValues[i] = dirichletBC(mesh.getNodePosition3(i));
		}
		for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1)
				continue;
			rhs[j] = integrateDual3Cell(mesh, i, sourceTerm, circumcentric);
			std::unordered_map<uint, double>& row = system_matrix[i];
			for (auto it = row.begin(); it != row.end(); ++it) {
				rhs[j] -= it->second * boundaryValues[it->first];
			}
			if (mesh.getNodeFlag(i) == 2)
				rhs[j] -= integrateNodeDualBoundary(mesh, i, neumannBC, circumcentric);
			++j;
		}
		//get the solution on active nodes by solving the system and on other nodes from the boundary condition
		Buffer<double> sol_interior;
		solveSystem(system_matrix, sol_interior, rhs, mesh);
		sol.resize(mesh.getNodeSize());
		for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1)
				sol[i] = dirichletBC(mesh.getNodePosition3(i));
			else
				sol[i] = sol_interior[j++];
		}
	}

	//Computes the solution using the refined mesh in polynomialWaveExampleWhitney2D and waveExampleWhitney2D.
	void getDECSolutionWaveEq(const BuilderMesh& mesh_old, const BuilderMesh& mesh, const SmallSimplexPartition2D& ssp, const std::function<double(Vector2)>& sourceTerm,
		const std::function<double(Vector2)>& dirichletBC, const std::function<Vector2(Vector2)>& neumannBC, Buffer<double>& sol) {
		Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize()); //form the system matrix
		{
			Buffer<std::unordered_map<uint, double>> star;
			Buffer<uint> refElementIndices(1, 0);
			Buffer <Buffer<uint>> refElements(1);
			refElements[0].resize(mesh_old.getFaceSize());
			for (uint i = 0; i < mesh_old.getFaceSize(); ++i)
				refElements[0][i] = i;
			ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "RightTriangle", false, true, false);
			for (uint j = 0; j < mesh.getNodeSize(); ++j) {
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
					}
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0) {
							system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
						}
					}
				}
			}
		}
		//count the number of active nodes
		uint activeNodes = 0;
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 0 || mesh.getNodeFlag(i) == 3)
				++activeNodes;
		}
		//form the right-hand side
		Buffer<double> rhs(activeNodes, 0.0);
		Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
				boundaryValues[i] = dirichletBC(mesh.getNodePosition2(i));
		}
		for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
				continue;
			rhs[j] = integrateDual2Cell(mesh, i, sourceTerm);
			std::unordered_map<uint, double>& row = system_matrix[i];
			for (auto it = row.begin(); it != row.end(); ++it) {
				rhs[j] -= it->second * boundaryValues[it->first];
			}
			if (mesh.getNodeFlag(i) == 2)
				rhs[j] -= integrateNodeDualBoundary(mesh, i, neumannBC, false, true);
			++j;
		}
		//get the solution on active nodes by solving the system and on other nodes from the boundary condition
		Buffer<double> sol_interior;

		solveSystem(system_matrix, sol_interior, rhs, mesh, 3, 2);

		sol.resize(mesh.getNodeSize());
		for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
				sol[i] = dirichletBC(mesh.getNodePosition2(i));
			else
				sol[i] = sol_interior[j++];
		}
	}

	//Computes the solution using the refined mesh in polynomialWaveExampleWhitney3D and waveExampleWhitney3D.
	void getDECSolutionWaveEq(const BuilderMesh& mesh_old, const BuilderMesh& mesh, const SmallSimplexPartition3D& ssp, const std::function<double(Vector3)>& sourceTerm,
		const std::function<double(Vector3)>& dirichletBC, const std::function<Vector3(Vector3)>& neumannBC, Buffer<double>& sol) {
		Buffer<std::unordered_map<uint, double>> system_matrix(mesh.getNodeSize()); //form the system matrix
		{
			Buffer<std::unordered_map<uint, double>> star;
			Buffer<uint> refElementIndices(6);
			Buffer<Buffer<uint>> refElements(6);
			for (uint i = 0; i < 6; ++i) {
				refElementIndices[i] = i;
				refElements[i].resize(mesh_old.getBodySize() / 6);
			}

			for (uint i = 0; i < mesh_old.getBodySize(); ++i)
				refElements[i % 6][i / 6] = i;
			ssp.formHodgeMatrix1FormsCustom(star, refElementIndices, refElements, "CubeTetrahedron", false, true, false);
			for (uint j = 0; j < mesh.getNodeSize(); ++j) {
				Buffer<double> col_j_of_star_d(mesh.getEdgeSize(), 0.0);
				const Buffer<uint>& edges_j = mesh.getNodeEdges(j);
				for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
					std::unordered_map<uint, double>& row = star[i];
					for (uint e = 0; e < edges_j.size(); ++e) {
						auto it = row.find(edges_j[e]);
						if (it != row.end())
							col_j_of_star_d[i] += it->second * mesh.getEdgeIncidence(edges_j[e], j);
					}
				}
				for (uint i = 0; i < mesh.getNodeSize(); ++i) {
					const Buffer<uint>& edges_i = mesh.getNodeEdges(i);
					for (uint e = 0; e < edges_i.size(); ++e) {
						if (col_j_of_star_d[edges_i[e]] != 0.0) {
							system_matrix[i][j] -= mesh.getEdgeIncidence(edges_i[e], i) * col_j_of_star_d[edges_i[e]];
						}
					}
				}
			}
		}
		//count the number of active nodes
		uint activeNodes = 0;
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 0 || mesh.getNodeFlag(i) == 3)
				++activeNodes;
		}
		//form the right-hand side
		Buffer<double> rhs(activeNodes, 0.0);
		Buffer<double> boundaryValues(mesh.getNodeSize(), 0.0);
		for (uint i = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
				boundaryValues[i] = dirichletBC(mesh.getNodePosition3(i));
		}
		for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 3)
				continue;
			rhs[j] = integrateDual3Cell(mesh, i, sourceTerm);
			std::unordered_map<uint, double>& row = system_matrix[i];
			for (auto it = row.begin(); it != row.end(); ++it) {
				rhs[j] -= it->second * boundaryValues[it->first];
			}
			if (mesh.getNodeFlag(i) == 2)
				rhs[j] -= integrateNodeDualBoundary(mesh, i, neumannBC, false, true);
			++j;
		}
		//get the solution on active nodes by solving the system and on other nodes from the boundary condition
		Buffer<double> sol_interior;

		solveSystem(system_matrix, sol_interior, rhs, mesh, 3, 2);

		sol.resize(mesh.getNodeSize());
		for (uint i = 0, j = 0; i < mesh.getNodeSize(); ++i) {
			if (mesh.getNodeFlag(i) == 1 || mesh.getNodeFlag(i) == 2)
				sol[i] = dirichletBC(mesh.getNodePosition3(i));
			else
				sol[i] = sol_interior[j++];
		}
	}
}