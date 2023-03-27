#include "Vector.hpp"
#include "Matrix.hpp"

using namespace gfd;

Matrix2 Vector2::outerProduct(const Vector2 &v) const { return Matrix2(x*v, y*v); }
SymMatrix2 Vector2::outerProduct() const { return SymMatrix2(x*x, x*y, y*y); }
SymTwoMatrix2 TwoVector2::outerProduct() const { return SymTwoMatrix2(xy*xy); }

Matrix3 Vector3::outerProduct(const Vector3 &v) const { return Matrix3(x*v, y*v, z*v); }
SymMatrix3 Vector3::outerProduct() const { return SymMatrix3(x*x, x*y, y*y, x*z, y*z, z*z); }
SymTwoMatrix3 TwoVector3::outerProduct() const { return SymTwoMatrix3(xy*xy, xy*xz, xz*xz, xy*yz, xz*yz, yz*yz); }
SymThreeMatrix3 ThreeVector3::outerProduct() const { return SymThreeMatrix3(xyz*xyz); }

Matrix4 Vector4::outerProduct(const Vector4 &v) const { return Matrix4(x*v, y*v, z*v, t*v); }
SymMatrix4 Vector4::outerProduct() const { return SymMatrix4(x*x, x*y, y*y, x*z, y*z, z*z, x*t, y*t, z*t, t*t); }
SymTwoMatrix4 TwoVector4::outerProduct() const
{
	return SymTwoMatrix4(xy*xy, xy*xz, xz*xz, xy*yz, xz*yz, yz*yz,
		xy*xt, xz*xt, yz*xt, xt*xt, xy*yt, xz*yt, yz*yt, xt*yt, yt*yt,
		xy*zt, xz*zt, yz*zt, xt*zt, yt*zt, zt*zt);
}
SymThreeMatrix4 ThreeVector4::outerProduct() const
{
	return SymThreeMatrix4(xyz*xyz, xyz*xyt, xyt*xyt, xyz*xzt, xyt*xzt, xzt*xzt, xyz*yzt, xyt*yzt, xzt*yzt, yzt*yzt);
}
SymFourMatrix4 FourVector4::outerProduct() const { return SymFourMatrix4(xyzt*xyzt); }

TwoVector3 Vector3::dual() const
{
	return TwoVector3(z, -y, x);
}
TwoVector3 Vector3::dualof() const
{
	return TwoVector3(z, -y, x);
}

ThreeVector4 Vector4::dual() const
{
	return ThreeVector4(-t, z, -y, x);
}
ThreeVector4 Vector4::dualof() const
{
	return ThreeVector4(t, -z, y, -x);
}

MatrixN VectorN::outerProduct(const VectorN &v) const
{
	MatrixN m;
	m.val.resize(size() > v.size() ? size() : v.size());
	for(uint i=0; i<size(); i++) m[i] = val[i] * v;
	for(uint i=0; i<m.size(); i++) m[i].toVectorN(m.size());
	return m;
}

double SymTwoMatrix4::determinant() const
{
	MatrixN m;
	m.toMatrixN(6);
	m[0][0] = xyxy; m[0][1] = xyxz; m[0][2] = xyyz; m[0][3] = xyxt; m[0][4] = xyyt; m[0][5] = xyzt;
	m[1][0] = xyxz; m[1][1] = xzxz; m[1][2] = xzyz; m[1][3] = xzxt; m[1][4] = xzyt; m[1][5] = xzzt;
	m[2][0] = xyyz; m[2][1] = xzyz; m[2][2] = yzyz; m[2][3] = yzxt; m[2][4] = yzyt; m[2][5] = yzzt;
	m[3][0] = xyxt; m[3][1] = xzxt; m[3][2] = yzxt; m[3][3] = xtxt; m[3][4] = xtyt; m[3][5] = xtzt;
	m[4][0] = xyyt; m[4][1] = xzyt; m[4][2] = yzyt; m[4][3] = xtyt; m[4][4] = ytyt; m[4][5] = ytzt;
	m[5][0] = xyzt; m[5][1] = xzzt; m[5][2] = yzzt; m[5][3] = xtzt; m[5][4] = ytzt; m[5][5] = ztzt;
	return m.determinant();
}

SymTwoMatrix4 SymTwoMatrix4::inverse(const double d) const
{
	MatrixN m;
	m.toMatrixN(6);
	m[0][0] = xyxy; m[0][1] = xyxz; m[0][2] = xyyz; m[0][3] = xyxt; m[0][4] = xyyt; m[0][5] = xyzt;
	m[1][0] = xyxz; m[1][1] = xzxz; m[1][2] = xzyz; m[1][3] = xzxt; m[1][4] = xzyt; m[1][5] = xzzt;
	m[2][0] = xyyz; m[2][1] = xzyz; m[2][2] = yzyz; m[2][3] = yzxt; m[2][4] = yzyt; m[2][5] = yzzt;
	m[3][0] = xyxt; m[3][1] = xzxt; m[3][2] = yzxt; m[3][3] = xtxt; m[3][4] = xtyt; m[3][5] = xtzt;
	m[4][0] = xyyt; m[4][1] = xzyt; m[4][2] = yzyt; m[4][3] = xtyt; m[4][4] = ytyt; m[4][5] = ytzt;
	m[5][0] = xyzt; m[5][1] = xzzt; m[5][2] = yzzt; m[5][3] = xtzt; m[5][4] = ytzt; m[5][5] = ztzt;
	const MatrixN n(m.inverse(d));
	return SymTwoMatrix4(n[0][0], n[1][0], n[1][1], n[2][0], n[2][1], n[2][2], n[3][0], n[3][1], n[3][2], n[3][3],
		n[4][0], n[4][1], n[4][2], n[4][3], n[4][4], n[5][0], n[5][1], n[5][2], n[5][3], n[5][4], n[5][5]);
}
