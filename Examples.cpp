#include "ErrorAnalysisFunctions.hpp"

using namespace gfd;

void example1();
void example2();
void example3();
void example4();
void example5();
void example6();
void example7();
void example8();
void example9();
void example10();
void example11();
void example12();
void example13();

int main(int argc, const char* argv[])
{
	//int exampleNumber = atoi(argv[1]);
	int exampleNumber = 1;
	switch (exampleNumber) {
	case 1:
		example1();
		break;
	case 2:
		example2();
		break;
	case 3:
		example3();
		break;
	case 4:
		example4();
		break;
	case 5:
		example5();
		break;
	case 6:
		example6();
		break;
	case 7:
		example7();
		break;
	case 8:
		example8();
		break;
	case 9:
		example9();
		break;
	case 10:
		example10();
		break;
	case 11:
		example11();
		break;
	case 12:
		example12();
		break;
	case 13:
		example13();
		break;
	}
}

void example1() {
	polynomialPoissonExampleWhitney2D('a', false);
	polynomialPoissonExampleWhitney2D('b', false);
	polynomialPoissonExampleWhitney2D('c', false);
}

void example2() {
	poissonExampleWhitney2D();
}

void example3() {
	polynomialPoissonExampleWhitney3D('a', false);
	polynomialPoissonExampleWhitney3D('b', false);
	polynomialPoissonExampleWhitney3D('c', false);
}

void example4() {
	poissonExampleWhitney3D();
}

void example5() {
	poissonExampleWhitney2D_Stability();
	poissonExampleWhitney3D_Stability();
}

void example6() {
	for (uint order = 1; order <= 6; ++order) {
		waveExampleWhitney2D_TimeStepping(order, 2.0, 100);
	}
	for (uint order = 1; order <= 6; ++order) {
		waveExampleWhitney2D_TimeStepping(order, 2.001, 20);
	}
}

void example7() {
	for (uint order = 1; order <= 6; ++order) {
		waveExampleWhitney2D_Stability(order, 2.0, 100);
	}
	for (uint order = 1; order <= 6; ++order) {
		waveExampleWhitney2D_Stability(order, 2.001, 20);
	}
	for (uint order = 1; order <= 6; ++order) {
		waveExampleWhitney2D_Stability(order, 2.0, 10000);
	}
	waveExampleWhitney2D_Stability(2, 5.0/3.0, 120);
}

void example8() {
	waveExampleCubical2D_TimeStepping(1, 2.0, 100);
	waveExampleCubical2D_TimeStepping(2, 1.0, 200);
	waveExampleCubical2D_TimeStepping(3, 1.0, 200);
	waveExampleCubical2D_TimeStepping(4, 200.0 / 211.0, 211);
	waveExampleCubical2D_TimeStepping(5, 200.0 / 211.0, 211);
}

void example9() {
	waveExampleWhitney3D_Stability(1, 1.41, 100);
	waveExampleWhitney3D_Stability(1, 1.42, 20);
	waveExampleWhitney3D_Stability(2, 0.25, 100);
	waveExampleWhitney3D_Stability(3, 0.25, 40);
	waveExampleWhitney3D_Stability(4, 0.25, 40);
}

void example10() {
	waveExampleCubical3D_Stability(1, 1.41, 100);
	waveExampleCubical3D_Stability(1, 1.42, 20);
	waveExampleCubical3D_Stability(2, 0.71, 200);
	waveExampleCubical3D_Stability(2, 0.75, 50);
	waveExampleCubical3D_Stability(3, 0.5, 1000);
	waveExampleCubical3D_Stability(3, 0.53, 200);
	waveExampleCubical3D_Stability(4, 0.4, 100);
	waveExampleCubical3D_Stability(4, 0.42, 100);
}

void example11() {
	waveExampleCubical3D_TimeStepping(1, 100.0 / 71.0, 71);
	waveExampleCubical3D_TimeStepping(2, 100.0 / 142.0, 142);
	waveExampleCubical3D_TimeStepping(3, 0.5, 200);
	waveExampleCubical3D_TimeStepping(4, 0.4, 250);
	waveExampleCubical3D_TimeStepping(5, 0.4, 250);
}

void example12() {
	polynomialTestWhitney0Forms3D();
	polynomialTestWhitney1Forms3D();
	polynomialTestWhitney2Forms3D();
	polynomialTestWhitney3Forms3D();
}

void example13() {
	testWhitney1Forms3D();
}