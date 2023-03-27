/*
Matrix.hpp implements 2, 3, and 4 dimensional square matrices.
The square matrix classes are called Matrix2, Matrix3, and Matrix4.

Matrix.hpp also implements 2, 3, and 4 dimensional symmetric matrices, two-matrices, three-matrices, and four-matrices,
The symmetric matrix classes are called SymMatrix2, SymMatrix3, and SymMatrix4.
The symmetric two-matrix classes are called SymmTwoMatrix2, SymmTwoMatrix3, and SymmTwoMatrix4.
The symmetric three-matrix classes are called SymmThreeMatrix3, and SymmThreeMatrix4.
The symmetric four-matrix is called SymmFourMatrix4.

Each n-matrices operate with n-vectors of corresponding dimension.
The classes implement basic matrix operators (+,-,*,/), transpose, normalization, determinant, and inverse.
*/

#ifndef _MATRIX_HPP_INCLUDED_
#define _MATRIX_HPP_INCLUDED_

#include "Types.hpp"
#include "Complex.hpp"

namespace gfd
{

template<class M, class V> class Matrix
{
public:
	// virtual functions
	virtual M copy() const { return M(*(M*)this); }
	virtual M transpose() const { return M(*(M*)this); }
	virtual M negation() const = 0;
	virtual M &increase(const M &m) = 0;
	virtual M &decrease(const M &m) = 0;
	virtual M &scale(const double d) = 0;
	virtual M product(const M &m) const = 0;
	virtual V vectorproduct(const V &v) const = 0;
	virtual M inverse(const double d = 1.0) const = 0;
	virtual bool equals(const M &m) const = 0;
	virtual double dot(const M &m) const = 0;
	virtual double determinant() const = 0;

	// shared operators and functions
	M operator+() const { return copy(); }
	M &operator+=(const M &m) { return increase(m); }
	M operator+(const M &m) const { return copy().increase(m); }

	M operator-() const { return negation(); }
	M &operator-=(const M &m) { return decrease(m); }
	M operator-(const M &m) const { return copy().decrease(m); }

	M &operator*=(const double d) { return scale(d); }
	M operator*(const double d) const { return copy().scale(d); }
	M &operator*=(const M &m) { return (*(M*)this) = product(m); }
	M operator*(const M &m) const { return product(m); }
	V operator*(const V &v) const { return vectorproduct(v); }
	friend M operator*(const double d, const M &m) { return m.copy().scale(d); }
	friend V operator*(const V &v, const M &m) { return m.transpose().vectorproduct(v); }

	M &operator/=(const double d) { return scale(1.0 / d); }
	M operator/(const double d) const { return copy().scale(1.0 / d); }
	M &operator/=(const M &m) { return (*(M*)this) = product(m.inverse()); }
	M operator/(const M &m) const { return product(m.inverse()); }
	friend M operator/(const double d, const M &m) { return m.inverse(d); }
	friend V operator/(const V &v, const M &m) { return m.inverse().transpose().vectorproduct(v); }

	bool operator==(const M &m) const { return equals(m); }
	bool operator!=(const M &m) const { return !equals(m); }

	double lensq() const { return dot(*(M*)this); }
	double len() const { return std::sqrt(lensq()); }
	M unit() const { return copy().scale(1.0 / len()); }
	M &normalize() { return scale(1.0 / len()); }
};

class SymMatrix2 : public Matrix<SymMatrix2,Vector2>
{
public:
	SymMatrix2() {}
	SymMatrix2(const double mxx, const double mxy, const double myy) { xx = mxx; xy = mxy; yy = myy; }
	SymMatrix2(const SymMatrix2 &m) { xx = m.xx; xy = m.xy; yy = m.yy; }

	SymMatrix2 negation() const { return SymMatrix2(-xx,-xy,-yy); }
	SymMatrix2 &increase(const SymMatrix2 &m) { xx += m.xx; xy += m.xy; yy += m.yy; return (*this); }
	SymMatrix2 &decrease(const SymMatrix2 &m) { xx -= m.xx; xy -= m.xy; yy -= m.yy; return (*this); }
	SymMatrix2 &scale(const double d) { xx *= d; xy *= d; yy *= d; return (*this); }
	Vector2 vectorproduct(const Vector2 &v) const { return Vector2(xx * v.x + xy * v.y, xy * v.x + yy * v.y); }
	SymMatrix2 product(const SymMatrix2 &m) const { return SymMatrix2(xx * m.xx + xy * m.xy, xx * m.xy + xy * m.yy, xy * m.xy + yy * m.yy); }
	bool equals(const SymMatrix2 &m) const { return (xx == m.xx) && (xy == m.xy) && (yy == m.yy); }
	double dot(const SymMatrix2 &m) const { return xx * m.xx + yy * m.yy + 2.0 * xy * m.xy; }
	double determinant() const { return xx * yy - xy * xy; }
	SymMatrix2 inverse(const double d = 1.0) const { return SymMatrix2(yy, -xy, xx) *= d / determinant(); }

	double xx;
	double xy;
	double yy;
};

class Matrix2 : public Matrix<Matrix2,Vector2>
{
public:
	Matrix2() {}
	Matrix2(const double mxx, const double mxy, const double myx, const double myy) : x(mxx, mxy), y(myx, myy) {}
	Matrix2(const Vector2 &mx, const Vector2 &my) : x(mx), y(my) {}
	Matrix2(const Matrix2 &m) : x(m.x), y(m.y) {}
	Matrix2(const SymMatrix2 &m) : x(m.xx, m.xy), y(m.xy, m.yy) {}

	Matrix2 transpose() const { return Matrix2(x.x,y.x, x.y,y.y); }
	Matrix2 negation() const { return Matrix2(-x,-y); }
	Matrix2 &increase(const Matrix2 &m) { x += m.x; y += m.y; return (*this); }
	Matrix2 &decrease(const Matrix2 &m) { x -= m.x; y -= m.y; return (*this); }
	Matrix2 &scale(const double d) { x *= d; y *= d; return (*this); }
	Vector2 vectorproduct(const Vector2 &v) const { return Vector2(x.dot(v), y.dot(v)); }
	Matrix2 product(const Matrix2 &m) const { return Matrix2(x.x * m.x + x.y * m.y, y.x * m.x + y.y * m.y); }
	bool equals(const Matrix2 &m) const { return (x == m.x) && (y == m.y); }
	double dot(const Matrix2 &m) const { return x.dot(m.x) + y.dot(m.y); }
	double determinant() const { return TwoVector2(x,y).xy; }
	Matrix2 inverse(const double d = 1.0) const
	{
		Matrix2 m(y.y, -x.y, -y.x, x.x);
		return m *= d / (x.x * m.x.x + x.y * m.y.x);
	}

	SymMatrix2 toSymMatrix2() const { return SymMatrix2(x.x, x.y, y.y); }

	Vector2 x;
	Vector2 y;
};

class SymTwoMatrix2 : public Matrix<SymTwoMatrix2,TwoVector2>
{
public:
	SymTwoMatrix2() {}
	SymTwoMatrix2(const double mxyxy) { xyxy = mxyxy; }
	SymTwoMatrix2(const SymTwoMatrix2 &m) { xyxy = m.xyxy; }
	SymTwoMatrix2(const SymMatrix2 &m) { xyxy = m.determinant(); }

	SymTwoMatrix2 negation() const { return SymTwoMatrix2(-xyxy); }
	SymTwoMatrix2 &increase(const SymTwoMatrix2 &m) { xyxy += m.xyxy; return (*this); }
	SymTwoMatrix2 &decrease(const SymTwoMatrix2 &m) { xyxy -= m.xyxy; return (*this); }
	SymTwoMatrix2 &scale(const double d) { xyxy *= d; return (*this); }
	TwoVector2 vectorproduct(const TwoVector2 &v) const { return TwoVector2(xyxy * v.xy); }
	SymTwoMatrix2 product(const SymTwoMatrix2 &m) const { return SymTwoMatrix2(xyxy * m.xyxy); }
	bool equals(const SymTwoMatrix2 &m) const { return (xyxy == m.xyxy); }
	double dot(const SymTwoMatrix2 &m) const { return xyxy * m.xyxy; }
	double determinant() const { return xyxy; }
	SymTwoMatrix2 inverse(const double d = 1.0) const { return SymTwoMatrix2(d / xyxy); }

	double xyxy;
};

class SymMatrix3 : public Matrix<SymMatrix3,Vector3>
{
public:
	SymMatrix3() {}
	SymMatrix3(const double mxx, const double mxy, const double myy, const double mxz, const double myz, const double mzz)
	{
		xx = mxx; xy = mxy; yy = myy; xz = mxz; yz = myz; zz = mzz;
	}
	SymMatrix3(const SymMatrix2 &m, const double mxz, const double myz, const double mzz)
	{
		xx = m.xx; xy = m.xy; yy = m.yy; xz = mxz; yz = myz; zz = mzz;
	}
	SymMatrix3(const SymMatrix3 &m)
	{
		xx = m.xx; xy = m.xy; yy = m.yy; xz = m.xz; yz = m.yz; zz = m.zz;
	}

	SymMatrix3 negation() const { return SymMatrix3(-xx,-xy,-yy,-xz,-yz,-zz); }
	SymMatrix3 &increase(const SymMatrix3 &m)
	{
		xx += m.xx; xy += m.xy; yy += m.yy; xz += m.xz; yz += m.yz; zz += m.zz;
		return (*this);
	}
	SymMatrix3 &decrease(const SymMatrix3 &m)
	{
		xx -= m.xx; xy -= m.xy; yy -= m.yy; xz -= m.xz; yz -= m.yz; zz -= m.zz;
		return (*this);
	}
	SymMatrix3 &scale(const double d)
	{
		xx *= d; xy *= d; yy *= d; xz *= d; yz *= d; zz *= d;
		return (*this);
	}
	Vector3 vectorproduct(const Vector3 &v) const
	{
		return Vector3(xx * v.x + xy * v.y + xz * v.z, xy * v.x + yy * v.y + yz * v.z, xz * v.x + yz * v.y + zz * v.z);
	}
	SymMatrix3 product(const SymMatrix3 &m) const
	{
		return SymMatrix3(xx * m.xx + xy * m.xy + xz * m.xz, xx * m.xy + xy * m.yy + xz * m.yz, xy * m.xy + yy * m.yy + yz * m.yz,
			xx * m.xz + xy * m.yz + xz * m.zz, xy * m.xz + yy * m.yz + yz * m.zz, xz * m.xz + yz * m.yz + zz * m.zz);
	}
	bool equals(const SymMatrix3 &m) const
	{
		return (xx == m.xx) && (xy == m.xy) && (yy == m.yy) && (xz == m.xz) && (yz == m.yz) && (zz == m.zz);
	}
	double dot(const SymMatrix3 &m) const
	{
		return xx * m.xx + yy * m.yy + zz * m.zz + 2.0 * (xy * m.xy + xz * m.xz + yz * m.yz);
	}
	double determinant() const
	{
		return xx * yy * zz + 2.0 * xy * yz * xz - yz * yz * xx - xz * xz * yy - xy * xy * zz;
	}
	SymMatrix3 inverse(const double d = 1.0) const
	{
		const Vector3 cross(xy * yz - xz * yy, xz * xy - xx * yz, xx * yy - xy * xy);
		return SymMatrix3(yy * zz - yz * yz, yz * xz - zz * xy, zz * xx - xz * xz, cross.x, cross.y, cross.z) *=
			d / (cross.x * xz + cross.y * yz + cross.z * zz);
	}

	SymMatrix2 toSymMatrix2() const { return SymMatrix2(xx, xy, yy); }

	double xx;
	double xy;
	double yy;
	double xz;
	double yz;
	double zz;
};

class Matrix3 : public Matrix<Matrix3,Vector3>
{
public:
	Matrix3() {}
	Matrix3(const double mxx, const double mxy, const double mxz,
		const double myx, const double myy, const double myz,
		const double mzx, const double mzy, const double mzz)
		: x(mxx, mxy, mxz), y(myx, myy, myz), z(mzx, mzy, mzz) {}
	Matrix3(const Vector3 &mx, const Vector3 &my, const Vector3 &mz) : x(mx), y(my), z(mz) {}
	Matrix3(const Matrix2 &m, const Vector2 &zcolumn, const Vector3 &mz) : x(m.x, zcolumn.x), y(m.y, zcolumn.y), z(mz) {}
	Matrix3(const Matrix3 &m) : x(m.x), y(m.y), z(m.z) {}
	Matrix3(const SymMatrix3 &m) : x(m.xx, m.xy, m.xz), y(m.xy, m.yy, m.yz), z(m.xz, m.yz, m.zz) {}

	Matrix3 transpose() const { return Matrix3(x.x,y.x,z.x, x.y,y.y,z.y, x.z,y.z,z.z); }
	Matrix3 negation() const { return Matrix3(-x,-y,z); }
	Matrix3 &increase(const Matrix3 &m) { x += m.x; y += m.y; z += m.z; return (*this); }
	Matrix3 &decrease(const Matrix3 &m) { x -= m.x; y -= m.y; z -= m.z; return (*this); }
	Matrix3 &scale(const double d) { x *= d; y *= d; z *= d; return (*this); }
	Vector3 vectorproduct(const Vector3 &v) const { return Vector3(x.dot(v), y.dot(v), z.dot(v)); }
	Matrix3 product(const Matrix3 &m) const { return Matrix3(x.x * m.x + x.y * m.y + x.z * m.z, y.x * m.x + y.y * m.y + y.z * m.z, z.x * m.x + z.y * m.y + z.z * m.z); }
	bool equals(const Matrix3 &m) const { return (x == m.x) && (y == m.y) && (z == m.z); }
	double dot(const Matrix3 &m) const { return x.dot(m.x) + y.dot(m.y) + z.dot(m.z); }
	double determinant() const { return ThreeVector3(x,y,z).xyz; }
	Matrix3 inverse(const double d = 1.0) const
	{
		const TwoVector3 vx(y, z);
		const TwoVector3 vy(x, z);
		const TwoVector3 vz(x, y);
		Matrix3 m(vx.yz,-vy.yz,vz.yz, -vx.xz,vy.xz,-vz.xz, vx.xy,-vy.xy,vz.xy);
		return m *= d / (x.x * m.x.x + x.y * m.y.x + x.z * m.z.x);
	}

	Matrix2 toMatrix2() const {	return Matrix2(x.x, x.y, y.x, y.y); }
	SymMatrix3 toSymMatrix3() const { return SymMatrix3(x.x, x.y, y.y, x.z, y.z, z.z); }

	Vector3 x;
	Vector3 y;
	Vector3 z;
};

class SymTwoMatrix3 : public Matrix<SymTwoMatrix3,TwoVector3>
{
public:
	SymTwoMatrix3() {}
	SymTwoMatrix3(const double mxyxy, const double mxyxz, const double mxzxz, const double mxyyz, const double mxzyz, const double myzyz)
	{
		xyxy = mxyxy; xyxz = mxyxz; xzxz = mxzxz; xyyz = mxyyz; xzyz = mxzyz; yzyz = myzyz;
	}
	SymTwoMatrix3(const SymTwoMatrix2 &m, const double mxyxz, const double mxzxz, const double mxyyz, const double mxzyz, const double myzyz)
	{
		xyxy = m.xyxy; xyxz = mxyxz; xzxz = mxzxz; xyyz = mxyyz; xzyz = mxzyz; yzyz = myzyz;
	}
	SymTwoMatrix3(const SymTwoMatrix3 &m)
	{
		xyxy = m.xyxy; xyxz = m.xyxz; xzxz = m.xzxz; xyyz = m.xyyz; xzyz = m.xzyz; yzyz = m.yzyz;
	}
	SymTwoMatrix3(const SymMatrix3 &m)
	{
		xyxy = m.xx * m.yy - m.xy * m.xy;
		xyxz = m.xx * m.yz - m.xz * m.xy;
		xzxz = m.xx * m.zz - m.xz * m.xz;
		xyyz = m.xy * m.yz - m.xz * m.yy;
		xzyz = m.xy * m.zz - m.xz * m.yz;
		yzyz = m.yy * m.zz - m.yz * m.yz;
	}

	SymTwoMatrix3 negation() const { return SymTwoMatrix3(-xyxy,-xyxz,-xzxz,-xyyz,-xzyz,-yzyz); }
	SymTwoMatrix3 &increase(const SymTwoMatrix3 &m)
	{
		xyxy += m.xyxy; xyxz += m.xyxz; xzxz += m.xzxz; xyyz += m.xyyz; xzyz += m.xzyz; yzyz += m.yzyz;
		return (*this);
	}
	SymTwoMatrix3 &decrease(const SymTwoMatrix3 &m)
	{
		xyxy -= m.xyxy; xyxz -= m.xyxz; xzxz -= m.xzxz; xyyz -= m.xyyz; xzyz -= m.xzyz; yzyz -= m.yzyz;
		return (*this);
	}
	SymTwoMatrix3 &scale(const double d)
	{
		xyxy *= d; xyxz *= d; xzxz *= d; xyyz *= d; xzyz *= d; yzyz *= d;
		return (*this);
	}
	TwoVector3 vectorproduct(const TwoVector3 &v) const
	{
		return TwoVector3(xyxy * v.xy + xyxz * v.xz + xyyz * v.yz, xyxz * v.xy + xzxz * v.xz + xzyz * v.yz, xyyz * v.xy + xzyz * v.xz + yzyz * v.yz);
	}
	SymTwoMatrix3 product(const SymTwoMatrix3 &m) const
	{
		return SymTwoMatrix3(xyxy * m.xyxy + xyxz * m.xyxz + xyyz * m.xyyz, xyxy * m.xyxz + xyxz * m.xzxz + xyyz * m.xzyz, xyxz * m.xyxz + xzxz * m.xzxz + xzyz * m.xzyz,
			xyxy * m.xyyz + xyxz * m.xzyz + xyyz * m.yzyz, xyxz * m.xyyz + xzxz * m.xzyz + xzyz * m.yzyz, xyyz * m.xyyz + xzyz * m.xzyz + yzyz * m.yzyz);
	}
	bool equals(const SymTwoMatrix3 &m) const
	{
		return (xyxy == m.xyxy) && (xyxz == m.xyxz) && (xzxz == m.xzxz) && (xyyz == m.xyyz) && (xzyz == m.xzyz) && (yzyz == m.yzyz);
	}
	double dot(const SymTwoMatrix3 &m) const
	{
		return xyxy * m.xyxy + xzxz * m.xzxz + yzyz * m.yzyz + 2.0 * (xyxz * m.xyxz + xyyz * m.xyyz + xzyz * m.xzyz);
	}
	double determinant() const
	{
		return xyxy * xzxz * yzyz + 2.0 * xyxz * xzyz * xyyz - xzyz * xzyz * xyxy - xyyz * xyyz * xzxz - xyxz * xyxz * yzyz;
	}
	SymTwoMatrix3 inverse(const double d = 1.0) const
	{
		const TwoVector3 cross(xyxz * xzyz - xyyz * xzxz, xyyz * xyxz - xyxy * xzyz, xyxy * xzxz - xyxz * xyxz);
		return SymTwoMatrix3(xzxz * yzyz - xzyz * xzyz, xzyz * xyyz - yzyz * xyxz, yzyz * xyxy - xyyz * xyyz, cross.xy, cross.xz, cross.yz) *=
			d / (cross.xy * xyyz + cross.xz * xzyz + cross.yz * yzyz);
	}

	SymTwoMatrix2 toSymTwoMatrix2() const { return SymTwoMatrix2(xyxy); }

	double xyxy;
	double xyxz;
	double xzxz;
	double xyyz;
	double xzyz;
	double yzyz;
};

class SymThreeMatrix3 : public Matrix<SymThreeMatrix3,ThreeVector3>
{
public:
	SymThreeMatrix3() {}
	SymThreeMatrix3(const double mxyzxyz) { xyzxyz = mxyzxyz; }
	SymThreeMatrix3(const SymThreeMatrix3 &m) { xyzxyz = m.xyzxyz; }
	SymThreeMatrix3(const SymMatrix3 &m) { xyzxyz = m.determinant(); }

	SymThreeMatrix3 negation() const { return SymThreeMatrix3(-xyzxyz); }
	SymThreeMatrix3 &increase(const SymThreeMatrix3 &m) { xyzxyz += m.xyzxyz; return (*this); }
	SymThreeMatrix3 &decrease(const SymThreeMatrix3 &m) { xyzxyz -= m.xyzxyz; return (*this); }
	SymThreeMatrix3 &scale(const double d) { xyzxyz *= d; return (*this); }
	ThreeVector3 vectorproduct(const ThreeVector3 &v) const { return ThreeVector3(xyzxyz * v.xyz); }
	SymThreeMatrix3 product(const SymThreeMatrix3 &m) const { return SymThreeMatrix3(xyzxyz * m.xyzxyz); }
	bool equals(const SymThreeMatrix3 &m) const { return (xyzxyz == m.xyzxyz); }
	double dot(const SymThreeMatrix3 &m) const { return xyzxyz * m.xyzxyz; }
	double determinant() const { return xyzxyz; }
	SymThreeMatrix3 inverse(const double d = 1.0) const { return SymThreeMatrix3(d / xyzxyz); }

	double xyzxyz;
};

class SymMatrix4 : public Matrix<SymMatrix4,Vector4>
{
public:
	SymMatrix4() {}
	SymMatrix4(const double mxx, const double mxy, const double myy, const double mxz, const double myz, const double mzz,
		const double mxt, const double myt, const double mzt, const double mtt)
	{
		xx = mxx; xy = mxy; yy = myy; xz = mxz; yz = myz; zz = mzz;
		xt = mxt; yt = myt; zt = mzt; tt = mtt;
	}
	SymMatrix4(const SymMatrix2 &m, const double mxz, const double myz, const double mzz,
		const double mxt, const double myt, const double mzt, const double mtt)
	{
		xx = m.xx; xy = m.xy; yy = m.yy; xz = mxz; yz = myz; zz = mzz;
		xt = mxt; yt = myt; zt = mzt; tt = mtt;
	}
	SymMatrix4(const SymMatrix3 &m, const double mxt, const double myt, const double mzt, const double mtt)
	{
		xx = m.xx; xy = m.xy; yy = m.yy; xz = m.xz; yz = m.yz; zz = m.zz;
		xt = mxt; yt = myt; zt = mzt; tt = mtt;
	}
	SymMatrix4(const SymMatrix4 &m)
	{
		xx = m.xx; xy = m.xy; yy = m.yy; xz = m.xz; yz = m.yz; zz = m.zz;
		xt = m.xt; yt = m.yt; zt = m.zt; tt = m.tt;
	}

	SymMatrix4 negation() const { return SymMatrix4(-xx,-xy,-yy,-xz,-yz,-zz,-xt,-yt,-zt,-tt); }
	SymMatrix4 &increase(const SymMatrix4 &m)
	{
		xx += m.xx; xy += m.xy; yy += m.yy; xz += m.xz; yz += m.yz; zz += m.zz;
		xt += m.xt; yt += m.yt; zt += m.zt; tt += m.tt;
		return (*this);
	}
	SymMatrix4 &decrease(const SymMatrix4 &m)
	{
		xx -= m.xx; xy -= m.xy; yy -= m.yy; xz -= m.xz; yz -= m.yz; zz -= m.zz;
		xt -= m.xt; yt -= m.yt; zt -= m.zt; tt -= m.tt;
		return (*this);
	}
	SymMatrix4 &scale(const double d)
	{
		xx *= d; xy *= d; yy *= d; xz *= d; yz *= d; zz *= d;
		xt *= d; yt *= d; zt *= d; tt *= d;
		return (*this);
	}
	Vector4 vectorproduct(const Vector4 &v) const
	{
		return Vector4(xx * v.x + xy * v.y + xz * v.z + xt * v.t, xy * v.x + yy * v.y + yz * v.z + yt * v.t,
			xz * v.x + yz * v.y + zz * v.z + zt * v.t, xt * v.x + yt * v.y + zt * v.z + tt * v.t);
	}
	SymMatrix4 product(const SymMatrix4 &m) const
	{
		return SymMatrix4(xx * m.xx + xy * m.xy + xz * m.xz + xt * m.xt,
			xx * m.xy + xy * m.yy + xz * m.yz + xt * m.yt,
			xy * m.xy + yy * m.yy + yz * m.yz + yt * m.yt,
			xx * m.xz + xy * m.yz + xz * m.zz + xt * m.zt,
			xy * m.xz + yy * m.yz + yz * m.zz + yt * m.zt,
			xz * m.xz + yz * m.yz + zz * m.zz + zt * m.zt,
			xx * m.xt + xy * m.yt + xz * m.zt + xt * m.tt,
			xy * m.xt + yy * m.yt + yz * m.zt + yt * m.tt,
			xz * m.xt + yz * m.yt + zz * m.zt + zt * m.tt,
			xt * m.xt + yt * m.yt + zt * m.zt + tt * m.tt);
	}
	bool equals(const SymMatrix4 &m) const
	{
		return (xx == m.xx) && (xy == m.xy) && (yy == m.yy) && (xz == m.xz) && (yz == m.yz) && (zz == m.zz) &&
			(xt == m.xt) && (yt == m.yt) && (zt == m.zt) && (tt == m.tt);
	}
	double dot(const SymMatrix4 &m) const
	{
		return xx * m.xx + yy * m.yy + zz * m.zz + tt * m.tt +
			2.0 * (xy * m.xy + xz * m.xz + yz * m.yz + xt * m.xt + yt * m.yt + zt * m.zt);
	}
	double determinant() const
	{
		const double dxy = xz * yt - xt * yz;
		const double dyz = yz * zt - yt * zz;
		const double dzt = zz * tt - zt * zt;
		const double dxz = xz * zt - xt * zz;
		const double dxt = xz * tt - xt * zt;
		const double dyt = yz * tt - yt * zt;
		return xx * (yy * dzt - yz * dyt + yt * dyz) - xy * (xy * dzt - yz * dxt + yt * dxz) +
			xz * (xy * dyt - yy * dxt + yt * dxy) - xt * (xy * dyz - yy * dxz + yz * dxy);
	}
	SymMatrix4 inverse(const double d = 1.0) const
	{
		const double dxy = xz * yt - xt * yz;
		const double dyz = yz * zt - yt * zz;
		const double dzt = zz * tt - zt * zt;
		const double dxt = xz * tt - xt * zt;
		const double dyt = yz * tt - yt * zt;
		const double dxz = xz * zt - xt * zz;
		const double dyzt = yy * dzt - yz * dyt + yt * dyz;
		const double dxzt = xy * dzt - yz * dxt + yt * dxz;
		const double dxyt = xy * dyt - yy * dxt + yt * dxy;
		const double dxyz = xy * dyz - yy * dxz + yz * dxy;
		const double axy = xx * yy - xy * xy;
		const double ayz = xy * yz - yy * xz;
		const double axz = xx * yz - xy * xz;
		const double ayt = xy * yt - yy * xt;
		const double axt = xx * yt - xy * xt;
		const double det = xx * dyzt - xy * dxzt + xz * dxyt - xt * dxyz;
		return SymMatrix4(dyzt, -dxzt, xx * dzt - xz * dxt + xt * dxz, dxyt, -xx * dyt + xy * dxt - xt * dxy, xt * ayt - yt * axt + tt * axy,
			-dxyz, xx * dyz - xy * dxz + xz * dxy, -xt * ayz + yt * axz - zt * axy, xz * ayz - yz * axz + zz * axy) *= d / det;
	}

	SymMatrix2 toSymMatrix2() const { return SymMatrix2(xx, xy, yy); }
	SymMatrix3 toSymMatrix3() const { return SymMatrix3(xx, xy, yy, xz, yz, zz); }

	double xx;
	double xy;
	double yy;
	double xz;
	double yz;
	double zz;
	double xt;
	double yt;
	double zt;
	double tt;
};

class Matrix4 : public Matrix<Matrix4,Vector4>
{
public:
	Matrix4() {}
	Matrix4(const double mxx, const double mxy, const double mxz, const double mxt,
					const double myx, const double myy, const double myz, const double myt,
					const double mzx, const double mzy, const double mzz, const double mzt,
					const double mtx, const double mty, const double mtz, const double mtt)
		: x(mxx, mxy, mxz, mxt), y(myx, myy, myz, myt), z(mzx, mzy, mzz, mzt), t(mtx, mty, mtz, mtt) {}
	Matrix4(const Vector4 &mx, const Vector4 &my, const Vector4 &mz, const Vector4 &mt) : x(mx), y(my), z(mz), t(mt) {}
	Matrix4(const Matrix2 &m, const Vector2 &zcolumn, const Vector2 &tcolumn, const Vector4 &mz, const Vector4 &mt)
		: x(m.x, zcolumn.x, tcolumn.x), y(m.y, zcolumn.y, tcolumn.y), z(mz), t(mt) {}
	Matrix4(const Matrix3 &m, const Vector3 &tcolumn, const Vector4 &mt)
		: x(m.x, tcolumn.x), y(m.y, tcolumn.y), z(m.z, tcolumn.z), t(mt) {}
	Matrix4(const Matrix4 &m) : x(m.x), y(m.y), z(m.z), t(m.t) {}
	Matrix4(const SymMatrix4 &m) : x(m.xx, m.xy, m.xz, m.xt), y(m.xy, m.yy, m.yz, m.yt), z(m.xz, m.yz, m.zz, m.zt), t(m.xt, m.yt, m.zt, m.tt) {}

	Matrix4 transpose() const { return Matrix4(x.x,y.x,z.x,t.x, x.y,y.y,z.y,t.y, x.z,y.z,z.z,t.z, x.t,y.t,z.t,t.t); }
	Matrix4 negation() const { return Matrix4(-x,-y,-z,-t); }
	Matrix4 &increase(const Matrix4 &m) { x += m.x; y += m.y; z += m.z; t += m.t; return (*this); }
	Matrix4 &decrease(const Matrix4 &m) { x -= m.x; y -= m.y; z -= m.z; t -= m.t; return (*this); }
	Matrix4 &scale(const double d) { x *= d; y *= d; z *= d; t *= d; return (*this); }
	Vector4 vectorproduct(const Vector4 &v) const { return Vector4(x.dot(v), y.dot(v), z.dot(v), t.dot(v)); }
	Matrix4 product(const Matrix4 &m) const
	{
		return Matrix4(x.x * m.x + x.y * m.y + x.z * m.z + x.t * m.t,
			y.x * m.x + y.y * m.y + y.z * m.z + y.t * m.t,
			z.x * m.x + z.y * m.y + z.z * m.z + z.t * m.t,
			t.x * m.x + t.y * m.y + t.z * m.z + t.t * m.t);
	}
	bool equals(const Matrix4 &m) const { return (x == m.x) && (y == m.y) && (z == m.z) && (t == m.t); }
	double dot(const Matrix4 &m) const { return x.dot(m.x) + y.dot(m.y) + z.dot(m.z) + t.dot(m.t); }
	double determinant() const { return FourVector4(x,y,z,t).xyzt; }
	Matrix4 inverse(const double d = 1.0) const
	{
		const ThreeVector4 vx(y, z, t);
		const ThreeVector4 vy(x, z, t);
		const ThreeVector4 vz(x, y, t);
		const ThreeVector4 vt(x, y, z);
		Matrix4 m(vx.yzt,-vy.yzt,vz.yzt,-vt.yzt, -vx.xzt,vy.xzt,-vz.xzt,vt.xzt, vx.xyt,-vy.xyt,vz.xyt,-vt.xyt, -vx.xyz,vy.xyz,-vz.xyz,vt.xyz);
		return m *= d / (x.x * m.x.x + x.y * m.y.x + x.z * m.z.x + x.t * m.t.x);
	}

	Matrix2 toMatrix2() const {	return Matrix2(x.x, x.y, y.x, y.y); }
	Matrix3 toMatrix3() const { return Matrix3(x.x, x.y, x.z, y.x, y.y, y.z, z.x, z.y, z.z); }
	SymMatrix4 toSymMatrix4() const { return SymMatrix4(x.x, x.y, y.y, x.z, y.z, z.z, x.t, y.t, z.t, t.t); }

	Vector4 x;
	Vector4 y;
	Vector4 z;
	Vector4 t;
};

class SymTwoMatrix4 : public Matrix<SymTwoMatrix4,TwoVector4>
{
public:
	SymTwoMatrix4() {}
	SymTwoMatrix4(const double mxyxy, const double mxyxz, const double mxzxz, const double mxyyz, const double mxzyz, const double myzyz,
		const double mxyxt, const double mxzxt, const double myzxt, const double mxtxt,
		const double mxyyt, const double mxzyt, const double myzyt, const double mxtyt, const double mytyt,
		const double mxyzt, const double mxzzt, const double myzzt, const double mxtzt, const double mytzt, const double mztzt)
	{
		xyxy = mxyxy; xyxz = mxyxz; xzxz = mxzxz; xyyz = mxyyz; xzyz = mxzyz; yzyz = myzyz;
		xyxt = mxyxt; xzxt = mxzxt; yzxt = myzxt; xtxt = mxtxt;
		xyyt = mxyyt; xzyt = mxzyt; yzyt = myzyt; xtyt = mxtyt; ytyt = mytyt;
		xyzt = mxyzt; xzzt = mxzzt; yzzt = myzzt; xtzt = mxtzt; ytzt = mytzt; ztzt = mztzt;
	}
	SymTwoMatrix4(const SymTwoMatrix2 &m, const double mxyxz, const double mxzxz, const double mxyyz, const double mxzyz, const double myzyz,
		const double mxyxt, const double mxzxt, const double myzxt, const double mxtxt,
		const double mxyyt, const double mxzyt, const double myzyt, const double mxtyt, const double mytyt,
		const double mxyzt, const double mxzzt, const double myzzt, const double mxtzt, const double mytzt, const double mztzt)
	{
		xyxy = m.xyxy; xyxz = mxyxz; xzxz = mxzxz; xyyz = mxyyz; xzyz = mxzyz; yzyz = myzyz;
		xyxt = mxyxt; xzxt = mxzxt; yzxt = myzxt; xtxt = mxtxt;
		xyyt = mxyyt; xzyt = mxzyt; yzyt = myzyt; xtyt = mxtyt; ytyt = mytyt;
		xyzt = mxyzt; xzzt = mxzzt; yzzt = myzzt; xtzt = mxtzt; ytzt = mytzt; ztzt = mztzt;
	}
	SymTwoMatrix4(const SymTwoMatrix3 &m, const double mxyxt, const double mxzxt, const double myzxt, const double mxtxt,
		const double mxyyt, const double mxzyt, const double myzyt, const double mxtyt, const double mytyt,
		const double mxyzt, const double mxzzt, const double myzzt, const double mxtzt, const double mytzt, const double mztzt)
	{
		xyxy = m.xyxy; xyxz = m.xyxz; xzxz = m.xzxz; xyyz = m.xyyz; xzyz = m.xzyz; yzyz = m.yzyz;
		xyxt = mxyxt; xzxt = mxzxt; yzxt = myzxt; xtxt = mxtxt;
		xyyt = mxyyt; xzyt = mxzyt; yzyt = myzyt; xtyt = mxtyt; ytyt = mytyt;
		xyzt = mxyzt; xzzt = mxzzt; yzzt = myzzt; xtzt = mxtzt; ytzt = mytzt; ztzt = mztzt;
	}
	SymTwoMatrix4(const SymTwoMatrix4 &m)
	{
		xyxy = m.xyxy; xyxz = m.xyxz; xzxz = m.xzxz; xyyz = m.xyyz; xzyz = m.xzyz; yzyz = m.yzyz;
		xyxt = m.xyxt; xzxt = m.xzxt; yzxt = m.yzxt; xtxt = m.xtxt;
		xyyt = m.xyyt; xzyt = m.xzyt; yzyt = m.yzyt; xtyt = m.xtyt; ytyt = m.ytyt;
		xyzt = m.xyzt; xzzt = m.xzzt; yzzt = m.yzzt; xtzt = m.xtzt; ytzt = m.ytzt; ztzt = m.ztzt;
	}
	SymTwoMatrix4(const SymMatrix4 &m)
	{
		xyxy = m.xx * m.yy - m.xy * m.xy;
		xyxz = m.xx * m.yz - m.xz * m.xy;
		xzxz = m.xx * m.zz - m.xz * m.xz;
		xyyz = m.xy * m.yz - m.xz * m.yy;
		xzyz = m.xy * m.zz - m.xz * m.yz;
		yzyz = m.yy * m.zz - m.yz * m.yz;
		xyxt = m.xx * m.yt - m.xt * m.xy;
		xzxt = m.xx * m.zt - m.xt * m.xz;
		yzxt = m.xy * m.zt - m.yt * m.xz;
		xtxt = m.xx * m.tt - m.xt * m.xt;
		xyyt = m.xy * m.yt - m.xt * m.yy;
		xzyt = m.xy * m.zt - m.xt * m.yz;
		yzyt = m.yy * m.zt - m.yt * m.yz;
		xtyt = m.xy * m.tt - m.xt * m.yt;
		ytyt = m.yy * m.tt - m.yt * m.yt;
		xyzt = m.xz * m.yt - m.xt * m.yz;
		xzzt = m.xz * m.zt - m.xt * m.zz;
		yzzt = m.yz * m.zt - m.yt * m.zz;
		xtzt = m.xz * m.tt - m.xt * m.zt;
		ytzt = m.yz * m.tt - m.yt * m.zt;
		ztzt = m.zz * m.tt - m.zt * m.zt;
	}

	SymTwoMatrix4 negation() const
	{
		return SymTwoMatrix4(-xyxy,-xyxz,-xzxz,-xyyz,-xzyz,-yzyz,-xyxt,-xzxt,-yzxt,-xtxt,-xyyt,-xzyt,-yzyt,-xtyt,-ytyt,-xyzt,-xzzt,-yzzt,-xtzt,-ytzt,-ztzt);
	}
	SymTwoMatrix4 &increase(const SymTwoMatrix4 &m)
	{
		xyxy += m.xyxy; xyxz += m.xyxz; xzxz += m.xzxz; xyyz += m.xyyz; xzyz += m.xzyz; yzyz += m.yzyz;
		xyxt += m.xyxt; xzxt += m.xzxt; yzxt += m.yzxt; xtxt += m.xtxt;
		xyyt += m.xyyt; xzyt += m.xzyt; yzyt += m.yzyt; xtyt += m.xtyt; ytyt += m.ytyt;
		xyzt += m.xyzt; xzzt += m.xzzt; yzzt += m.yzzt; xtzt += m.xtzt; ytzt += m.ytzt; ztzt += m.ztzt;
		return (*this);
	}
	SymTwoMatrix4 &decrease(const SymTwoMatrix4 &m)
	{
		xyxy -= m.xyxy; xyxz -= m.xyxz; xzxz -= m.xzxz; xyyz -= m.xyyz; xzyz -= m.xzyz; yzyz -= m.yzyz;
		xyxt -= m.xyxt; xzxt -= m.xzxt; yzxt -= m.yzxt; xtxt -= m.xtxt;
		xyyt -= m.xyyt; xzyt -= m.xzyt; yzyt -= m.yzyt; xtyt -= m.xtyt; ytyt -= m.ytyt;
		xyzt -= m.xyzt; xzzt -= m.xzzt; yzzt -= m.yzzt; xtzt -= m.xtzt; ytzt -= m.ytzt; ztzt -= m.ztzt;
		return (*this);
	}
	SymTwoMatrix4 &scale(const double d)
	{
		xyxy *= d; xyxz *= d; xzxz *= d; xyyz *= d; xzyz *= d; yzyz *= d;
		xyxt *= d; xzxt *= d; yzxt *= d; xtxt *= d;
		xyyt *= d; xzyt *= d; yzyt *= d; xtyt *= d; ytyt *= d;
		xyzt *= d; xzzt *= d; yzzt *= d; xtzt *= d; ytzt *= d; ztzt *= d;
		return (*this);
	}
	TwoVector4 vectorproduct(const TwoVector4 &v) const
	{
		return TwoVector4(xyxy * v.xy + xyxz * v.xz + xyyz * v.yz + xyxt * v.xt + xyyt * v.yt + xyzt * v.zt,
			xyxz * v.xy + xzxz * v.xz + xzyz * v.yz + xzxt * v.xt + xzyt * v.yt + xzzt * v.zt,
			xyyz * v.xy + xzyz * v.xz + yzyz * v.yz + yzxt * v.xt + yzyt * v.yt + yzzt * v.zt,
			xyxt * v.xy + xzxt * v.xz + yzxt * v.yz + xtxt * v.xt + xtyt * v.yt + xtzt * v.zt,
			xyyt * v.xy + xzyt * v.xz + yzyt * v.yz + xtyt * v.xt + ytyt * v.yt + ytzt * v.zt,
			xyzt * v.xy + xzzt * v.xz + yzzt * v.yz + xtzt * v.xt + ytzt * v.yt + ztzt * v.zt);
	}
	SymTwoMatrix4 product(const SymTwoMatrix4 &m) const
	{
		return SymTwoMatrix4(xyxy * m.xyxy + xyxz * m.xyxz + xyyz * m.xyyz + xyxt * m.xyxt + xyyt * m.xyyt + xyzt * m.xyzt,
			xyxy * m.xyxz + xyxz * m.xzxz + xyyz * m.xzyz + xyxt * m.xzxt + xyyt * m.xzyt + xyzt * m.xzzt,
			xyxz * m.xyxz + xzxz * m.xzxz + xzyz * m.xzyz + xzxt * m.xzxt + xzyt * m.xzyt + xzzt * m.xzzt,
			xyxy * m.xyyz + xyxz * m.xzyz + xyyz * m.yzyz + xyxt * m.yzxt + xyyt * m.yzyt + xyzt * m.yzzt,
			xyxz * m.xyyz + xzxz * m.xzyz + xzyz * m.yzyz + xzxt * m.yzxt + xzyt * m.yzyt + xzzt * m.yzzt,
			xyyz * m.xyyz + xzyz * m.xzyz + yzyz * m.yzyz + yzxt * m.yzxt + yzyt * m.yzyt + yzzt * m.yzzt,
			xyxy * m.xyxt + xyxz * m.xzxt + xyyz * m.yzxt + xyxt * m.xtxt + xyyt * m.xtyt + xyzt * m.xtzt,
			xyxz * m.xyxt + xzxz * m.xzxt + xzyz * m.yzxt + xzxt * m.xtxt + xzyt * m.xtyt + xzzt * m.xtzt,
			xyyz * m.xyxt + xzyz * m.xzxt + yzyz * m.yzxt + yzxt * m.xtxt + yzyt * m.xtyt + yzzt * m.xtzt,
			xyxt * m.xyxt + xzxt * m.xzxt + yzxt * m.yzxt + xtxt * m.xtxt + xtyt * m.xtyt + xtzt * m.xtzt,
			xyxy * m.xyyt + xyxz * m.xzyt + xyyz * m.yzyt + xyxt * m.xtyt + xyyt * m.ytyt + xyzt * m.ytzt,
			xyxz * m.xyyt + xzxz * m.xzyt + xzyz * m.yzyt + xzxt * m.xtyt + xzyt * m.ytyt + xzzt * m.ytzt,
			xyyz * m.xyyt + xzyz * m.xzyt + yzyz * m.yzyt + yzxt * m.xtyt + yzyt * m.ytyt + yzzt * m.ytzt,
			xyxt * m.xyyt + xzxt * m.xzyt + yzxt * m.yzyt + xtxt * m.xtyt + xtyt * m.ytyt + xtzt * m.ytzt,
			xyyt * m.xyyt + xzyt * m.xzyt + yzyt * m.yzyt + xtyt * m.xtyt + ytyt * m.ytyt + ytzt * m.ytzt,
			xyxy * m.xyzt + xyxz * m.xzzt + xyyz * m.yzzt + xyxt * m.xtzt + xyyt * m.ytzt + xyzt * m.ztzt,
			xyxz * m.xyzt + xzxz * m.xzzt + xzyz * m.yzzt + xzxt * m.xtzt + xzyt * m.ytzt + xzzt * m.ztzt,
			xyyz * m.xyzt + xzyz * m.xzzt + yzyz * m.yzzt + yzxt * m.xtzt + yzyt * m.ytzt + yzzt * m.ztzt,
			xyxt * m.xyzt + xzxt * m.xzzt + yzxt * m.yzzt + xtxt * m.xtzt + xtyt * m.ytzt + xtzt * m.ztzt,
			xyyt * m.xyzt + xzyt * m.xzzt + yzyt * m.yzzt + xtyt * m.xtzt + ytyt * m.ytzt + ytzt * m.ztzt,
			xyzt * m.xyzt + xzzt * m.xzzt + yzzt * m.yzzt + xtzt * m.xtzt + ytzt * m.ytzt + ztzt * m.ztzt);
	}
	bool equals(const SymTwoMatrix4 &m) const
	{
		return (xyxy == m.xyxy) && (xyxz == m.xyxz) && (xzxz == m.xzxz) && (xyyz == m.xyyz) && (xzyz == m.xzyz) && (yzyz == m.yzyz) &&
			(xyxt == m.xyxt) && (xzxt == m.xzxt) && (yzxt == m.yzxt) && (xtxt == m.xtxt) &&
			(xyyt == m.xyyt) && (xzyt == m.xzyt) && (yzyt == m.yzyt) && (xtyt == m.xtyt) && (ytyt == m.ytyt) &&
			(xyzt == m.xyzt) && (xzzt == m.xzzt) && (yzzt == m.yzzt) && (xtzt == m.xtzt) && (ytzt == m.ytzt) && (ztzt == m.ztzt);
	}
	double dot(const SymTwoMatrix4 &m) const
	{
		return xyxy * m.xyxy + xzxz * m.xzxz + yzyz * m.yzyz + xtxt * m.xtxt + ytyt * m.ytyt + ztzt * m.ztzt +
			2.0 * (xyxz * m.xyxz + xyyz * m.xyyz + xzyz * m.xzyz + xyxt * m.xyxt + xzxt * m.xzxt + yzxt * m.yzxt +
			xyyt * m.xyyt + xzyt * m.xzyt + yzyt * m.yzyt + xtyt * m.xtyt +
			xyzt * m.xyzt + xzzt * m.xzzt + yzzt * m.yzzt + xtzt * m.xtzt + ytzt * m.ytzt);
	}
	double determinant() const;
	SymTwoMatrix4 inverse(const double d = 1.0) const;

	SymTwoMatrix2 toSymTwoMatrix2() const { return SymTwoMatrix2(xyxy); }
	SymTwoMatrix3 toSymTwoMatrix3() const { return SymTwoMatrix3(xyxy, xyxz, xzxz, xyyz, xzyz, yzyz); }

	double xyxy;
	double xyxz;
	double xzxz;
	double xyyz;
	double xzyz;
	double yzyz;
	double xyxt;
	double xzxt;
	double yzxt;
	double xtxt;
	double xyyt;
	double xzyt;
	double yzyt;
	double xtyt;
	double ytyt;
	double xyzt;
	double xzzt;
	double yzzt;
	double xtzt;
	double ytzt;
	double ztzt;

};


class SymThreeMatrix4 : public Matrix<SymThreeMatrix4,ThreeVector4>
{
public:
	SymThreeMatrix4() {}
	SymThreeMatrix4(const double mxyzxyz, const double mxyzxyt, const double mxytxyt, const double mxyzxzt, const double mxytxzt, const double mxztxzt,
		const double mxyzyzt, const double mxytyzt, const double mxztyzt, const double myztyzt)
	{
		xyzxyz = mxyzxyz; xyzxyt = mxyzxyt; xytxyt = mxytxyt; xyzxzt = mxyzxzt; xytxzt = mxytxzt; xztxzt = mxztxzt;
		xyzyzt = mxyzyzt; xytyzt = mxytyzt; xztyzt = mxztyzt; yztyzt = myztyzt;
	}
	SymThreeMatrix4(const SymThreeMatrix3 &m, const double mxyzxyt, const double mxytxyt, const double mxyzxzt, const double mxytxzt, const double mxztxzt,
		const double mxyzyzt, const double mxytyzt, const double mxztyzt, const double myztyzt)
	{
		xyzxyz = m.xyzxyz; xyzxyt = mxyzxyt; xytxyt = mxytxyt; xyzxzt = mxyzxzt; xytxzt = mxytxzt; xztxzt = mxztxzt;
		xyzyzt = mxyzyzt; xytyzt = mxytyzt; xztyzt = mxztyzt; yztyzt = myztyzt;
	}
	SymThreeMatrix4(const SymThreeMatrix4 &m)
	{
		xyzxyz = m.xyzxyz; xyzxyt = m.xyzxyt; xytxyt = m.xytxyt; xyzxzt = m.xyzxzt; xytxzt = m.xytxzt; xztxzt = m.xztxzt;
		xyzyzt = m.xyzyzt; xytyzt = m.xytyzt; xztyzt = m.xztyzt; yztyzt = m.yztyzt;
	}
	SymThreeMatrix4(const SymMatrix4 &m)
	{
		xyzxyz = SymMatrix3(m.xx, m.xy, m.yy, m.xz, m.yz, m.zz).determinant();
		xyzxyt = Matrix3(m.xx, m.xy, m.xz, m.xy, m.yy, m.yz, m.xt, m.yt, m.zt).determinant();
		xytxyt = SymMatrix3(m.xx, m.xy, m.yy, m.xt, m.yt, m.tt).determinant();
		xyzxzt = Matrix3(m.xx, m.xy, m.xz, m.xz, m.yz, m.zz, m.xt, m.yt, m.zt).determinant();
		xytxzt = Matrix3(m.xx, m.xy, m.xt, m.xz, m.yz, m.zt, m.xt, m.yt, m.tt).determinant();
		xztxzt = SymMatrix3(m.xx, m.xz, m.zz, m.xt, m.zt, m.tt).determinant();
		xyzyzt = Matrix3(m.xy, m.yy, m.yz, m.xz, m.yz, m.zz, m.xt, m.yt, m.zt).determinant();
		xytyzt = Matrix3(m.xy, m.yy, m.yt, m.xz, m.yz, m.zt, m.xt, m.yt, m.tt).determinant();
		xztyzt = Matrix3(m.xy, m.yz, m.yt, m.xz, m.zz, m.zt, m.xt, m.zt, m.tt).determinant();
		yztyzt = SymMatrix3(m.yy, m.yz, m.zz, m.yt, m.zt, m.tt).determinant();
	}

	SymThreeMatrix4 negation() const { return SymThreeMatrix4(-xyzxyz,-xyzxyt,-xytxyt,-xyzxzt,-xytxzt,-xztxzt,-xyzyzt,-xytyzt,-xztyzt,-yztyzt); }
	SymThreeMatrix4 &increase(const SymThreeMatrix4 &m)
	{
		xyzxyz += m.xyzxyz; xyzxyt += m.xyzxyt; xytxyt += m.xytxyt; xyzxzt += m.xyzxzt; xytxzt += m.xytxzt; xztxzt += m.xztxzt;
		xyzyzt += m.xyzyzt; xytyzt += m.xytyzt; xztyzt += m.xztyzt; yztyzt += m.yztyzt;
		return (*this);
	}
	SymThreeMatrix4 &decrease(const SymThreeMatrix4 &m)
	{
		xyzxyz -= m.xyzxyz; xyzxyt -= m.xyzxyt; xytxyt -= m.xytxyt; xyzxzt -= m.xyzxzt; xytxzt -= m.xytxzt; xztxzt -= m.xztxzt;
		xyzyzt -= m.xyzyzt; xytyzt -= m.xytyzt; xztyzt -= m.xztyzt; yztyzt -= m.yztyzt;
		return (*this);
	}
	SymThreeMatrix4 &scale(const double d)
	{
		xyzxyz *= d; xyzxyt *= d; xytxyt *= d; xyzxzt *= d; xytxzt *= d; xztxzt *= d;
		xyzyzt *= d; xytyzt *= d; xztyzt *= d; yztyzt *= d;
		return (*this);
	}
	ThreeVector4 vectorproduct(const ThreeVector4 &v) const
	{
		return ThreeVector4(xyzxyz * v.xyz + xyzxyt * v.xyt + xyzxzt * v.xzt + xyzyzt * v.yzt, xyzxyt * v.xyz + xytxyt * v.xyt + xytxzt * v.xzt + xytyzt * v.yzt,
			xyzxzt * v.xyz + xytxzt * v.xyt + xztxzt * v.xzt + xztyzt * v.yzt, xyzyzt * v.xyz + xytyzt * v.xyt + xztyzt * v.xzt + yztyzt * v.yzt);
	}
	SymThreeMatrix4 product(const SymThreeMatrix4 &m) const
	{
		return SymThreeMatrix4(xyzxyz * m.xyzxyz + xyzxyt * m.xyzxyt + xyzxzt * m.xyzxzt + xyzyzt * m.xyzyzt,
			xyzxyz * m.xyzxyt + xyzxyt * m.xytxyt + xyzxzt * m.xytxzt + xyzyzt * m.xytyzt,
			xyzxyt * m.xyzxyt + xytxyt * m.xytxyt + xytxzt * m.xytxzt + xytyzt * m.xytyzt,
			xyzxyz * m.xyzxzt + xyzxyt * m.xytxzt + xyzxzt * m.xztxzt + xyzyzt * m.xztyzt,
			xyzxyt * m.xyzxzt + xytxyt * m.xytxzt + xytxzt * m.xztxzt + xytyzt * m.xztyzt,
			xyzxzt * m.xyzxzt + xytxzt * m.xytxzt + xztxzt * m.xztxzt + xztyzt * m.xztyzt,
			xyzxyz * m.xyzyzt + xyzxyt * m.xytyzt + xyzxzt * m.xztyzt + xyzyzt * m.yztyzt,
			xyzxyt * m.xyzyzt + xytxyt * m.xytyzt + xytxzt * m.xztyzt + xytyzt * m.yztyzt,
			xyzxzt * m.xyzyzt + xytxzt * m.xytyzt + xztxzt * m.xztyzt + xztyzt * m.yztyzt,
			xyzyzt * m.xyzyzt + xytyzt * m.xytyzt + xztyzt * m.xztyzt + yztyzt * m.yztyzt);
	}
	bool equals(const SymThreeMatrix4 &m) const
	{
		return (xyzxyz == m.xyzxyz) && (xyzxyt == m.xyzxyt) && (xytxyt == m.xytxyt) && (xyzxzt == m.xyzxzt) && (xytxzt == m.xytxzt) && (xztxzt == m.xztxzt) &&
			(xyzyzt == m.xyzyzt) && (xytyzt == m.xytyzt) && (xztyzt == m.xztyzt) && (yztyzt == m.yztyzt);
	}
	double dot(const SymThreeMatrix4 &m) const
	{
		return xyzxyz * m.xyzxyz + xytxyt * m.xytxyt + xztxzt * m.xztxzt + yztyzt * m.yztyzt +
			2.0 * (xyzxyt * m.xyzxyt + xyzxzt * m.xyzxzt + xytxzt * m.xytxzt + xyzyzt * m.xyzyzt + xytyzt * m.xytyzt + xztyzt * m.xztyzt);
	}
	double determinant() const
	{
		const double dxyzxyt = xyzxzt * xytyzt - xyzyzt * xytxzt;
		const double dxytxzt = xytxzt * xztyzt - xytyzt * xztxzt;
		const double dxztyzt = xztxzt * yztyzt - xztyzt * xztyzt;
		const double dxyzxzt = xyzxzt * xztyzt - xyzyzt * xztxzt;
		const double dxyzyzt = xyzxzt * yztyzt - xyzyzt * xztyzt;
		const double dxytyzt = xytxzt * yztyzt - xytyzt * xztyzt;
		return xyzxyz * (xytxyt * dxztyzt - xytxzt * dxytyzt + xytyzt * dxytxzt) - xyzxyt * (xyzxyt * dxztyzt - xytxzt * dxyzyzt + xytyzt * dxyzxzt) +
			xyzxzt * (xyzxyt * dxytyzt - xytxyt * dxyzyzt + xytyzt * dxyzxyt) - xyzyzt * (xyzxyt * dxytxzt - xytxyt * dxyzxzt + xytxzt * dxyzxyt);
	}
	SymThreeMatrix4 inverse(const double d = 1.0) const
	{
		const double dxyzxyt = xyzxzt * xytyzt - xyzyzt * xytxzt;
		const double dxytxzt = xytxzt * xztyzt - xytyzt * xztxzt;
		const double dxztyzt = xztxzt * yztyzt - xztyzt * xztyzt;
		const double dxyzyzt = xyzxzt * yztyzt - xyzyzt * xztyzt;
		const double dxytyzt = xytxzt * yztyzt - xytyzt * xztyzt;
		const double dxyzxzt = xyzxzt * xztyzt - xyzyzt * xztxzt;
		const double dxytxztyzt = xytxyt * dxztyzt - xytxzt * dxytyzt + xytyzt * dxytxzt;
		const double dxyzxztyzt = xyzxyt * dxztyzt - xytxzt * dxyzyzt + xytyzt * dxyzxzt;
		const double dxyzxytyzt = xyzxyt * dxytyzt - xytxyt * dxyzyzt + xytyzt * dxyzxyt;
		const double dxyzxytxzt = xyzxyt * dxytxzt - xytxyt * dxyzxzt + xytxzt * dxyzxyt;
		const double axyzxyt = xyzxyz * xytxyt - xyzxyt * xyzxyt;
		const double axytxzt = xyzxyt * xytxzt - xytxyt * xyzxzt;
		const double axyzxzt = xyzxyz * xytxzt - xyzxyt * xyzxzt;
		const double axytyzt = xyzxyt * xytyzt - xytxyt * xyzyzt;
		const double axyzyzt = xyzxyz * xytyzt - xyzxyt * xyzyzt;
		const double det = xyzxyz * dxytxztyzt - xyzxyt * dxyzxztyzt + xyzxzt * dxyzxytyzt - xyzyzt * dxyzxytxzt;
		return SymThreeMatrix4(dxytxztyzt, -dxyzxztyzt, xyzxyz * dxztyzt - xyzxzt * dxyzyzt + xyzyzt * dxyzxzt, dxyzxytyzt, -xyzxyz * dxytyzt + xyzxyt * dxyzyzt - xyzyzt * dxyzxyt, xyzyzt * axytyzt - xytyzt * axyzyzt + yztyzt * axyzxyt,
			-dxyzxytxzt, xyzxyz * dxytxzt - xyzxyt * dxyzxzt + xyzxzt * dxyzxyt, -xyzyzt * axytxzt + xytyzt * axyzxzt - xztyzt * axyzxyt, xyzxzt * axytxzt - xytxzt * axyzxzt + xztxzt * axyzxyt) *= d / det;
	}

	SymThreeMatrix3 toSymThreeMatrix3() const { return SymThreeMatrix3(xyzxyz); }

	double xyzxyz;
	double xyzxyt;
	double xytxyt;
	double xyzxzt;
	double xytxzt;
	double xztxzt;
	double xyzyzt;
	double xytyzt;
	double xztyzt;
	double yztyzt;
};

class SymFourMatrix4 : public Matrix<SymFourMatrix4,FourVector4>
{
public:
	SymFourMatrix4() {}
	SymFourMatrix4(const double mxyztxyzt) { xyztxyzt = mxyztxyzt; }
	SymFourMatrix4(const SymFourMatrix4 &m) { xyztxyzt = m.xyztxyzt; }
	SymFourMatrix4(const SymMatrix4 &m) { xyztxyzt = m.determinant(); }

	SymFourMatrix4 negation() const { return SymFourMatrix4(-xyztxyzt); }
	SymFourMatrix4 &increase(const SymFourMatrix4 &m) { xyztxyzt += m.xyztxyzt; return (*this); }
	SymFourMatrix4 &decrease(const SymFourMatrix4 &m) { xyztxyzt -= m.xyztxyzt; return (*this); }
	SymFourMatrix4 &scale(const double d) { xyztxyzt *= d; return (*this); }
	FourVector4 vectorproduct(const FourVector4 &v) const { return FourVector4(xyztxyzt * v.xyzt); }
	SymFourMatrix4 product(const SymFourMatrix4 &m) const { return SymFourMatrix4(xyztxyzt * m.xyztxyzt); }
	bool equals(const SymFourMatrix4 &m) const { return (xyztxyzt == m.xyztxyzt); }
	double dot(const SymFourMatrix4 &m) const { return xyztxyzt * m.xyztxyzt; }
	double determinant() const { return xyztxyzt; }
	SymFourMatrix4 inverse(const double d = 1.0) const { return SymFourMatrix4(d / xyztxyzt); }

	double xyztxyzt;
};

class MatrixN : public Matrix<MatrixN,VectorN>
{
public:
	MatrixN() {}
	MatrixN(const double m) {
		val.resize(1);
		val[0] = VectorN(m);
	}
	MatrixN(const Matrix2 &m) {
		val.resize(2);
		val[0] = VectorN(m.x); val[1] = VectorN(m.y);
	}
	MatrixN(const Matrix3 &m) {
		val.resize(3);
		val[0] = VectorN(m.x); val[1] = VectorN(m.y); val[2] = VectorN(m.z);
	}
	MatrixN(const Matrix4 &m) {
		val.resize(4);
		val[0] = VectorN(m.x); val[1] = VectorN(m.y); val[2] = VectorN(m.z); val[3] = VectorN(m.t);
	}
	MatrixN(const Buffer< Buffer<double> > &mval) {
		toMatrixN(mval.size());
		for(uint i=0; i<mval.size(); i++) {
			if(mval[i].size() > mval.size()) toMatrixN(mval[i].size());
			for(uint j=0; j<mval[i].size(); j++) val[i][j] = mval[i][j];
		}
	}
	MatrixN(const Buffer<VectorN> &mval) {
		toMatrixN(mval.size());
		for(uint i=0; i<mval.size(); i++) {
			if(mval[i].size() > mval.size()) toMatrixN(mval[i].size());
			val[i] = mval[i];
			val[i].toVectorN(size());
		}
	}
	MatrixN(const MatrixN &m) {
		val.resize(m.size());
		for(uint i=0; i<size(); i++) val[i] = m.val[i];
	}
	MatrixN &operator=(const MatrixN &m) {
		val.resize(m.size());
		for(uint i=0; i<size(); i++) val[i] = m.val[i];
		return (*this);
	}

	MatrixN transpose() const {
		MatrixN m;
		m.val.resize(size());
		for(uint i=0; i<size(); i++) {
			m[i].val.resize(size());
			for(uint j=0; j<size(); j++) m[i][j] = val[j][i];
		}
		return m;
	}
	MatrixN negation() const {
		MatrixN m;
		m.val.resize(size());
		for(uint i=0; i<size(); i++) m[i] = val[i].negation();
		return m;
	}
	MatrixN &increase(const MatrixN &m) {
		if(m.size() > size()) toMatrixN(m.size());
		for(uint i=0; i<m.size(); i++) val[i].increase(m[i]);
		return (*this);
	}
	MatrixN &decrease(const MatrixN &m) {
		if(m.size() > size()) toMatrixN(m.size());
		for(uint i=0; i<m.size(); i++) val[i].decrease(m[i]);
		return (*this);
	}
	MatrixN &scale(const double d) {
		for(uint i=0; i<size(); i++) val[i].scale(d);
		return (*this);
	}
	VectorN vectorproduct(const VectorN &v) const {
		VectorN r;
		r.val.resize(size());
		for(uint i=0; i<size(); i++) r[i] = val[i].dot(v);
		return r;
	}
	MatrixN product(const MatrixN &m) const {
		MatrixN r;
		r.toMatrixN(size() > m.size() ? size() : m.size());
		const uint minsize = size() < m.size() ? size() : m.size();
		for(uint i=0; i<size(); i++) {
			for(uint j=0; j<minsize; j++) r[i] += val[i][j] * m[j];
		}
		return r;
	}
	bool equals(const MatrixN &m) const {
		if(size() != m.size()) return false;
		for(uint i=0; i<size(); i++) {
			if(val[i] != m[i]) return false;
		}
		return true;
	}
	double dot(const MatrixN &m) const {
		const uint minsize = (size() < m.size() ? size() : m.size());
		double sum = 0.0;
		for(uint i=0; i<minsize; i++) sum += val[i].dot(m[i]);
		return sum;
	}
	double determinant() const {
		MatrixN l, u;
		Buffer<uint> p;
		if(!decomposeLU(l, u, p)) return 0.0;

		uint i, j;
		double prod = 1.0; // determinant of l equals one
		for(j=1; j<size(); j++) { // determinant of permutation
			for(i=0; i<j; i++) {
				if(p[i] > p[j]) prod *= -1.0;
			}
		}
		for(i=0; i<size(); i++) prod *= u[i][i]; // determinant of u
		return prod;
	}
	MatrixN inverse(const double d = 1.0) const {
		MatrixN l, u;
		Buffer<uint> p;
		if(!decomposeLU(l, u, p)) return MatrixN();

		// invert l, u, and p
		uint i, j, k;
		const uint n = size();
		for(i=0; i<n; i++) {
			for(j=i+1; j<n; j++) {
				double &lji = l[j][i];
				lji *= -1.0;
				for(k=0; k<i; k++) l[j][k] += lji * l[i][k];
			}
		}
		for(i=n; i-->0; ) {
			double div = 1.0 / u[i][i];
			for(j=0; j<i; j++) {
				double &uji = u[j][i];
				uji *= -div;
				for(k=i+1; k<n; k++) u[j][k] += uji * u[i][k];
			}
			div *= d;
			u[i][i] = div;
			for(j=i+1; j<n; j++) u[i][j] *= div;
		}
		MatrixN pu;
		pu.val.resize(size());
		for(i=0; i<n; i++) pu[i].val.swap(u[p[i]].val);
		return pu * l;
	}
	uint size() const { return val.size(); }
	VectorN &operator [](const uint i) const { return val[i]; }
	void toMatrixN(const uint size) {
		Buffer<VectorN> buf(size);
		for(uint i=0; i<size; i++) {
			if(i < val.size()) buf[i] = val[i];
			buf[i].toVectorN(size);
		}
		val.swap(buf);
	}
	MatrixN &permutate(const Buffer<uint> &p)
	{
		if(p.size() != size()) return *this;
		for(uint i=0; i<size(); i++) {
			const VectorN v(val[i]);
			for(uint j=0; j<size(); j++) {
				val[i][j] = v[p[j]];
			}
		}
		return *this;
	}
	bool decomposeLU(MatrixN &l, MatrixN &u, Buffer<uint> &p) const {
		// resize l and u
		uint i, j, k;
		const uint n = size();

		l.val.resize(n);
		u = *this;
		p.resize(n);
		for(i=0; i<n; i++) {
			l[i].val.resize(n);
			p[i] = i;
		}

		// decompose
		for(i=0; i<n; i++)
		{
			// find largest absolute value of |u[i]|
			k = i;
			double kabs = std::abs(u[i][k]);
			for(j=i+1; j<n; j++) {
				const double jabs = std::abs(u[i][j]);
				if(jabs > kabs) {
					k = j;
					kabs = jabs;
				}
			}
			if(kabs < 1e-15) return false; // matrix is not invertible

			if(i != k) { // apply permutation
				for(j=0; j<n; j++) {
					const double ujk = u[j][k];
					u[j][k] = u[j][i];
					u[j][i] = ujk;
					if(p[j] == i) p[j] = k;
					else if(p[j] == k) p[j] = i;
				}
			}

			// modify matrices l and u
			l[i][i] = 1.0;
			for(j=i+1; j<n; j++) {

				double &lji = l[j][i];
				lji = u[j][i] / u[i][i];
				l[i][j] = 0.0;
				u[j][i] = 0.0;
				for(k=i+1; k<n; k++) u[j][k] -= lji * u[i][k];
			}
		}
		return true; // return true only if the matrix is invertible
	}

	Buffer<VectorN> val;

};



const Matrix2 IDENTITYMATRIX2(1,0,0,1);
const Matrix2 ZEROMATRIX2(0,0,0,0);

const Matrix3 IDENTITYMATRIX3(1,0,0,0,1,0,0,0,1);
const Matrix3 ZEROMATRIX3(0,0,0,0,0,0,0,0,0);

const Matrix4 IDENTITYMATRIX4(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
const Matrix4 ZEROMATRIX4(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

const SymMatrix2 IDENTITYSYMMATRIX2(1,0,1);
const SymMatrix2 ZEROSYMMATRIX2(0,0,0);

const SymTwoMatrix2 IDENTITYSYMTWOMATRIX2(1);
const SymTwoMatrix2 ZEROSYMTWOMATRIX2(0);

const SymMatrix3 IDENTITYSYMMATRIX3(1,0,1,0,0,1);
const SymMatrix3 ZEROSYMMATRIX3(0,0,0,0,0,0);

const SymTwoMatrix3 IDENTITYSYMTWOMATRIX3(1,0,1,0,0,1);
const SymTwoMatrix3 ZEROSYMTWOMATRIX3(0,0,0,0,0,0);

const SymThreeMatrix3 IDENTITYSYMTHREEMATRIX3(1);
const SymThreeMatrix3 ZEROSYMTHREEMATRIX3(0);

const SymMatrix4 IDENTITYSYMMATRIX4(1,0,1,0,0,1,0,0,0,1);
const SymMatrix4 ZEROSYMMATRIX4(0,0,0,0,0,0,0,0,0,0);

const SymTwoMatrix4 IDENTITYSYMTWOMATRIX4(1,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1);
const SymTwoMatrix4 ZEROSYMTWOMATRIX4(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

const SymThreeMatrix4 IDENTITYSYMTHREEMATRIX4(1,0,1,0,0,1,0,0,0,1);
const SymThreeMatrix4 ZEROSYMTHREEMATRIX4(0,0,0,0,0,0,0,0,0,0);

const SymFourMatrix4 IDENTITYSYMFOURMATRIX4(1);
const SymFourMatrix4 ZEROSYMFOURMATRIX4(0);

}

#endif //_MATRIX_HPP_INCLUDED_
