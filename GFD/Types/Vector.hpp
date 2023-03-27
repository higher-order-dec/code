/*
Vector.hpp implements 2, 3, and 4 dimensional vectors, two-vectors (= wedge product of two vectors), three-vectors, and four-vectors.
The space components are denoted by x, y, z, and t, respectively.
The vector classes are called Vector2, Vector3, and Vector4.
The two-vector classes are called TwoVector2, TwoVector3, and TwoVector4.
The three-vector classes are called ThreeVector3, and ThreeVector4.
The four-vector class is called FourVector4.
The classes implement basic vector operators (+,-,*,/), dot product, wedge product, Hodge dual, length, and normalization.
Vectors can be considered as diagonal matrices. Thus, the classes implement also diagonal matrix product, inverse, transpose, and determinant.
*/

#ifndef _VECTOR_HPP_INCLUDED_
#define _VECTOR_HPP_INCLUDED_

#include <cmath>
#include "Types.hpp"
#include "Buffer.hpp"

namespace gfd
{

template<class V> class Vector
{
public:
	// virtual functions
	virtual V copy() const { return V(*(V*)this); }
	virtual V transpose() const { return V(*(V*)this); }
	virtual V negation() const = 0;
	virtual V &increase(const V &v) = 0;
	virtual V &decrease(const V &v) = 0;
	virtual V &scale(const double d) = 0;
	virtual V product(const V &v) const = 0;
	virtual V inverse(const double d = 1.0) const = 0;
	virtual bool equals(const V &v) const = 0;
	virtual double dot(const V &v) const = 0;
	virtual double determinant() const = 0;

	// shared operators and functions
	V operator+() const { return copy(); }
	V &operator+=(const V &v) { return increase(v); }
	V operator+(const V &v) const { return copy().increase(v); }

	V operator-() const { return negation(); }
	V &operator-=(const V &v) { return decrease(v); }
	V operator-(const V &v) const { return copy().decrease(v); }

	V &operator*=(const double d) { return scale(d); }
	V operator*(const double d) const { return copy().scale(d); }
	V &operator*=(const V &v) { return (*(V*)this) = product(v); }
	V operator*(const V &v) const { return product(v); }
	friend V operator*(const double d, const V &v) { return v.copy().scale(d); }

	V &operator/=(const double d) { return scale(1.0 / d); }
	V operator/(const double d) const { return copy().scale(1.0 / d); }
	V &operator/=(const V &v) { return (*(V*)this) = product(v.inverse()); }
	V operator/(const V &v) const { return product(v.inverse()); }
	friend V operator/(const double d, const V &v) { return v.inverse(d); }

	bool operator==(const V &v) const { return equals(v); }
	bool operator!=(const V &v) const { return !equals(v); }

	double lensq() const { return dot(*(V*)this); }
	double len() const { return std::sqrt(lensq()); }
	V unit() const { return copy().scale(1.0 / len()); }
	V &normalize() { return scale(1.0 / len()); }
};

class TwoVector2;
class Matrix2;
class SymMatrix2;

class Vector2 : public Vector<Vector2>
{
public:
	Vector2() {}
	Vector2(const double vx, const double vy) {	x = vx;	y = vy;	}
	Vector2(const Vector2 &v) {	x = v.x; y = v.y; }

	Vector2 negation() const { return Vector2(-x, -y); }
	Vector2 &increase(const Vector2 &v) { x += v.x; y += v.y; return (*this); }
	Vector2 &decrease(const Vector2 &v) { x -= v.x; y -= v.y; return (*this); }
	Vector2 &scale(const double d) { x *= d; y *= d; return (*this); }
	Vector2 product(const Vector2 &v) const { return Vector2(x * v.x, y * v.y); }
	Vector2 inverse(const double d = 1.0) const { return Vector2(d / x, d / y); }
	bool equals(const Vector2 &v) const { return (x == v.x) && (y == v.y); }
	double dot(const Vector2 &v) const { return x * v.x + y * v.y; }
	double determinant() const { return x * y; }

	Vector2 dual() const { return Vector2(-y, x); }
	Vector2 dualof() const { return Vector2(y, -x); }

	Matrix2 outerProduct(const Vector2 &v) const;
	SymMatrix2 outerProduct() const;

	double x;
	double y;
};

class SymTwoMatrix2;

class TwoVector2 : public Vector<TwoVector2>
{
public:
	TwoVector2() {}
	TwoVector2(const double vxy) { xy = vxy; }
	TwoVector2(const TwoVector2 &v) { xy = v.xy; }
	TwoVector2(const Vector2 &v1, const Vector2 &v2) { xy = v1.x * v2.y - v1.y * v2.x; }

	TwoVector2 negation() const { return TwoVector2(-xy); }
	TwoVector2 &increase(const TwoVector2 &v) { xy += v.xy; return (*this); }
	TwoVector2 &decrease(const TwoVector2 &v) { xy -= v.xy; return (*this); }
	TwoVector2 &scale(const double d) { xy *= d; return (*this); }
	TwoVector2 product(const TwoVector2 &v) const { return TwoVector2(xy * v.xy); }
	TwoVector2 inverse(const double d = 1.0) const { return TwoVector2(d / xy); }
	bool equals(const TwoVector2 &v) const { return (xy == v.xy); }
	double dot(const TwoVector2 &v) const { return xy * v.xy; }
	double determinant() const { return xy; }

	double dual() const { return xy; }
	double dualof() const { return xy; }

	SymTwoMatrix2 outerProduct() const;

	double xy;
};

class TwoVector3;
class Matrix3;
class SymMatrix3;

class Vector3 : public Vector<Vector3>
{
public:
	Vector3() {}
	Vector3(const double vx, const double vy, const double vz) {	x = vx;	y = vy; z = vz;	}
	Vector3(const Vector2 &v, const double vz) { x = v.x; y = v.y; z = vz; }
	Vector3(const Vector3 &v) {	x = v.x; y = v.y; z = v.z; }

	Vector3 negation() const { return Vector3(-x, -y, -z); }
	Vector3 &increase(const Vector3 &v) { x += v.x; y += v.y; z += v.z; return (*this); }
	Vector3 &decrease(const Vector3 &v) { x -= v.x; y -= v.y; z -= v.z; return (*this); }
	Vector3 &scale(const double d) { x *= d; y *= d; z *= d; return (*this); }
	Vector3 product(const Vector3 &v) const { return Vector3(x * v.x, y * v.y, z * v.z); }
	Vector3 inverse(const double d = 1.0) const { return Vector3(d / x, d / y, d / z); }
	bool equals(const Vector3 &v) const { return (x == v.x) && (y == v.y) && (z == v.z); }
	double dot(const Vector3 &v) const { return x * v.x + y * v.y + z * v.z; }
	double determinant() const { return x * y * z; }

	TwoVector3 dual() const;
	TwoVector3 dualof() const;
	Vector2 toVector2() const {	return Vector2(x, y); }

	Matrix3 outerProduct(const Vector3 &v) const;
	SymMatrix3 outerProduct() const;

	double x;
	double y;
	double z;
};

class SymTwoMatrix3;

class TwoVector3 : public Vector<TwoVector3>
{
public:
	TwoVector3() {}
	TwoVector3(const double vxy, const double vxz, const double vyz) { xy = vxy; xz = vxz; yz = vyz; }
	TwoVector3(const TwoVector2 &v, const double vxz, const double vyz) { xy = v.xy; xz = vxz; yz = vyz; }
	TwoVector3(const TwoVector3 &v) { xy = v.xy; xz = v.xz; yz = v.yz; }
	TwoVector3(const double vxy, const Vector2 &vt) { xy = vxy; xz = vt.x; yz = vt.y; }
	TwoVector3(const Vector3 &v1, const Vector3 &v2)
	{
		xy = v1.x * v2.y - v1.y * v2.x;
		xz = v1.x * v2.z - v1.z * v2.x;
		yz = v1.y * v2.z - v1.z * v2.y;
	}

	TwoVector3 negation() const { return TwoVector3(-xy, -xz, -yz); }
	TwoVector3 &increase(const TwoVector3 &v) { xy += v.xy; xz += v.xz; yz += v.yz; return (*this); }
	TwoVector3 &decrease(const TwoVector3 &v) { xy -= v.xy; xz -= v.xz; yz -= v.yz; return (*this); }
	TwoVector3 &scale(const double d) { xy *= d; xz *= d; yz *= d; return (*this); }
	TwoVector3 product(const TwoVector3 &v) const { return TwoVector3(xy * v.xy, xz * v.xz, yz * v.yz); }
	TwoVector3 inverse(const double d = 1.0) const { return TwoVector3(d / xy, d / xz, d / yz); }
	bool equals(const TwoVector3 &v) const { return (xy == v.xy) && (xz == v.xz) && (yz == v.yz); }
	double dot(const TwoVector3 &v) const { return xy * v.xy + xz * v.xz + yz * v.yz; }
	double determinant() const { return xy * xz * yz; }

	Vector3 dual() const { return Vector3(yz, -xz, xy); }
	Vector3 dualof() const { return Vector3(yz, -xz, xy); }
	TwoVector2 toTwoVector2() const { return TwoVector2(xy); }

	SymTwoMatrix3 outerProduct() const;

	double xy;
	double xz;
	double yz;
};

class SymThreeMatrix3;

class ThreeVector3 : public Vector<ThreeVector3>
{
public:
	ThreeVector3() {}
	ThreeVector3(const double vxyz) { xyz = vxyz; }
	ThreeVector3(const ThreeVector3 &v) { xyz = v.xyz; }
	ThreeVector3(const TwoVector3 &v1, const Vector3 &v2) { initWedge(v1, v2); }
	ThreeVector3(const Vector3 &v1, const TwoVector3 &v2) { initWedge(v2, v1); }
	ThreeVector3(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3) { initWedge(TwoVector3(v1,v2),v3); }

	ThreeVector3 negation() const { return ThreeVector3(-xyz); }
	ThreeVector3 &increase(const ThreeVector3 &v) { xyz += v.xyz; return (*this); }
	ThreeVector3 &decrease(const ThreeVector3 &v) { xyz -= v.xyz; return (*this); }
	ThreeVector3 &scale(const double d) { xyz *= d; return (*this); }
	ThreeVector3 product(const ThreeVector3 &v) const { return ThreeVector3(v.xyz * xyz); }
	ThreeVector3 inverse(const double d = 1.0) const { return ThreeVector3(d / xyz); }
	bool equals(const ThreeVector3 &v) const { return (xyz == v.xyz); }
	double dot(const ThreeVector3 &v) const { return xyz * v.xyz; }
	double determinant() const { return xyz; }

	double dual() const { return xyz; }
	double dualof() const { return xyz; }

	SymThreeMatrix3 outerProduct() const;

	double xyz;
private:
	void initWedge(const TwoVector3 &v1, const Vector3 &v2) { xyz = v1.xy * v2.z - v1.xz * v2.y + v1.yz * v2.x; }
};

class TwoVector4;
class ThreeVector4;
class FourVector4;
class Matrix4;
class SymMatrix4;

class Vector4 : public Vector<Vector4>
{
public:
	Vector4() {}
	Vector4(const double vx, const double vy, const double vz, const double vt) { x = vx; y = vy; z = vz; t = vt; }
	Vector4(const Vector2 &v, const double vz, const double vt) { x = v.x; y = v.y; z = vz; t = vt; }
	Vector4(const Vector3 &v, const double vt) { x = v.x; y = v.y; z = v.z; t = vt; }
	Vector4(const Vector4 &v) { x = v.x; y = v.y; z = v.z; t = v.t; }

	Vector4 negation() const { return Vector4(-x, -y, -z, -t); }
	Vector4 &increase(const Vector4 &v) { x += v.x; y += v.y; z += v.z; t += v.t; return (*this); }
	Vector4 &decrease(const Vector4 &v) { x -= v.x; y -= v.y; z -= v.z; t -= v.t; return (*this); }
	Vector4 &scale(const double d) { x *= d; y *= d; z *= d; t *= d; return (*this); }
	Vector4 product(const Vector4 &v) const { return Vector4(x * v.x, y * v.y, z * v.z, t * v.t); }
	Vector4 inverse(const double d = 1.0) const { return Vector4(d / x, d / y, d / z, d / t); }
	bool equals(const Vector4 &v) const { return (x == v.x) && (y == v.y) && (z == v.z) && (t == v.t); }
	double dot(const Vector4 &v) const { return x * v.x + y * v.y + z * v.z + t * v.t; }
	double determinant() const { return x * y * z * t; }

	ThreeVector4 dual() const;
	ThreeVector4 dualof() const;
	Vector2 toVector2() const {	return Vector2(x, y); }
	Vector3 toVector3() const {	return Vector3(x, y, z); }

	Matrix4 outerProduct(const Vector4 &v) const;
	SymMatrix4 outerProduct() const;

	double x;
	double y;
	double z;
	double t;
};

class SymTwoMatrix4;

class TwoVector4 : public Vector<TwoVector4>
{
public:
	TwoVector4() {}
	TwoVector4(const double vxy, const double vxz, const double vyz, const double vxt, const double vyt, const double vzt) { xy = vxy; xz = vxz; yz = vyz; xt = vxt; yt = vyt; zt = vzt; }
	TwoVector4(const TwoVector2 &v, const double vxz, const double vyz, const double vxt, const double vyt, const double vzt) { xy = v.xy; xz = vxz; yz = vyz; xt = vxt; yt = vyt; zt = vzt; }
	TwoVector4(const TwoVector3 &v, const double vxt, const double vyt, const double vzt) { xy = v.xy; xz = v.xz; yz = v.yz; xt = vxt; yt = vyt; zt = vzt; }
	TwoVector4(const TwoVector4 &v) { xy = v.xy; xz = v.xz; yz = v.yz; xt = v.xt; yt = v.yt; zt = v.zt; }
	TwoVector4(const double vxy, const double vxz, const double vyz, const Vector3 &vt) { xy = vxy; xz = vxz; yz = vyz; xt = vt.x; yt = vt.y; zt = vt.z; }
	TwoVector4(const Vector4 &v1, const Vector4 &v2)
	{
		xy = v1.x * v2.y - v1.y * v2.x;
		xz = v1.x * v2.z - v1.z * v2.x;
		yz = v1.y * v2.z - v1.z * v2.y;
		xt = v1.x * v2.t - v1.t * v2.x;
		yt = v1.y * v2.t - v1.t * v2.y;
		zt = v1.z * v2.t - v1.t * v2.z;
	}

	TwoVector4 negation() const { return TwoVector4(-xy,-xz,-yz,-xt,-yt,-zt); }
	TwoVector4 &increase(const TwoVector4 &v) { xy += v.xy; xz += v.xz; yz += v.yz; xt += v.xt; yt += v.yt; zt += v.zt; return (*this); }
	TwoVector4 &decrease(const TwoVector4 &v) { xy -= v.xy; xz -= v.xz; yz -= v.yz; xt -= v.xt; yt -= v.yt; zt -= v.zt; return (*this); }
	TwoVector4 &scale(const double d) { xy *= d; xz *= d; yz *= d; xt *= d; yt *= d; zt *= d; return (*this); }
	TwoVector4 product(const TwoVector4 &v) const { return TwoVector4(xy * v.xy, xz * v.xz, yz * v.yz, xt * v.xt, yt * v.yt, zt * v.zt); }
	TwoVector4 inverse(const double d = 1.0) const { return TwoVector4(d / xy, d / xz, d / yz, d / xt, d / yt, d / zt); }
	bool equals(const TwoVector4 &v) const { return (xy == v.xy) && (xz == v.xz) && (yz == v.yz) && (xt == v.xt) && (yt == v.yt) && (zt == v.zt); }
	double dot(const TwoVector4 &v) const { return xy * v.xy + xz * v.xz + yz * v.yz + xt * v.xt + yt * v.yt + zt * v.zt; }
	double determinant() const { return xy * xz * yz * xt * yt * zt; }

	TwoVector4 dual() const { return TwoVector4(zt, -yt, xt, yz, -xz, xy); }
	TwoVector4 dualof() const { return TwoVector4(zt, -yt, xt, yz, -xz, xy); }
	TwoVector2 toTwoVector2() const { return TwoVector2(xy); }
	TwoVector3 toTwoVector3() const { return TwoVector3(xy, xz, yz); }

	SymTwoMatrix4 outerProduct() const;

	double xy;
	double xz;
	double yz;
	double xt;
	double yt;
	double zt;
};

class SymThreeMatrix4;

class ThreeVector4 : public Vector<ThreeVector4>
{
public:
	ThreeVector4() {}
	ThreeVector4(const double vxyz, const double vxyt, const double vxzt, const double vyzt) { xyz = vxyz; xyt = vxyt; xzt = vxzt; yzt = vyzt; }
	ThreeVector4(const ThreeVector3 &v, const double vxyt, const double vxzt, const double vyzt) { xyz = v.xyz; xyt = vxyt; xzt = vxzt; yzt = vyzt; }
	ThreeVector4(const ThreeVector4 &v) { xyz = v.xyz; xyt = v.xyt; xzt = v.xzt; yzt = v.yzt; }
	ThreeVector4(const double vxyz, const TwoVector3 &vt) { xyz = vxyz; xyt = vt.xy; xzt = vt.xz; yzt = vt.yz; }
	ThreeVector4(const double vxyz, const double vxyt, const Vector2 &vzt) { xyz = vxyz; xyt = vxyt; xzt = vzt.x; yzt = vzt.y; }
	ThreeVector4(const TwoVector4 &v1, const Vector4 &v2) { initWedge(v1, v2); }
	ThreeVector4(const Vector4 &v1, const TwoVector4 &v2) { initWedge(v2, v1); }
	ThreeVector4(const Vector4 &v1, const Vector4 &v2, const Vector4 &v3) { initWedge(TwoVector4(v1,v2),v3); }


	ThreeVector4 negation() const { return ThreeVector4(-xyz,-xyt,-xzt,-yzt); }
	ThreeVector4 &increase(const ThreeVector4 &v) { xyz += v.xyz; xyt += v.xyt; xzt += v.xzt; yzt += v.yzt; return (*this); }
	ThreeVector4 &decrease(const ThreeVector4 &v) { xyz -= v.xyz; xyt -= v.xyt; xzt -= v.xzt; yzt -= v.yzt; return (*this); }
	ThreeVector4 &scale(const double d) { xyz *= d; xyt *= d; xzt *= d; yzt *= d; return (*this); }
	ThreeVector4 product(const ThreeVector4 &v) const { return ThreeVector4(xyz * v.xyz, xyt * v.xyt, xzt * v.xzt, yzt * v.yzt); }
	ThreeVector4 inverse(const double d = 1.0) const { return ThreeVector4(d / xyz, d / xyt, d / xzt, d / yzt); }
	bool equals(const ThreeVector4 &v) const { return (xyz == v.xyz) && (xyt == v.xyt) && (xzt == v.xzt) && (yzt == v.yzt); }
	double dot(const ThreeVector4 &v) const { return xyz * v.xyz + xyt * v.xyt + xzt * v.xzt + yzt * v.yzt; }
	double determinant() const { return xyz * xyt * xzt * yzt; }

	Vector4 dual() const { return Vector4(-yzt, xzt, -xyt, xyz); }
	Vector4 dualof() const { return Vector4(yzt, -xzt, xyt, -xyz); }
	ThreeVector3 toThreeVector3() const { return ThreeVector3(xyz); }

	SymThreeMatrix4 outerProduct() const;

	double xyz;
	double xyt;
	double xzt;
	double yzt;
private:
	void initWedge(const TwoVector4 &v1, const Vector4 &v2)
	{
		xyz = v1.xy * v2.z - v1.xz * v2.y + v1.yz * v2.x;
		xyt = v1.xy * v2.t - v1.xt * v2.y + v1.yt * v2.x;
		xzt = v1.xz * v2.t - v1.xt * v2.z + v1.zt * v2.x;
		yzt = v1.yz * v2.t - v1.yt * v2.z + v1.zt * v2.y;
	}
};

class SymFourMatrix4;

class FourVector4 : public Vector<FourVector4>
{
public:
	FourVector4() {}
	FourVector4(const double vxyzt) { xyzt = vxyzt; }
	FourVector4(const FourVector4 &v) { xyzt = v.xyzt; }
	FourVector4(const ThreeVector4 &v1, const Vector4 &v2) { initWedge(v1, v2); }
	FourVector4(const TwoVector4 &v1, const TwoVector4 &v2) { xyzt = v1.xy * v2.zt - v1.xz * v2.yt + v1.yz * v2.xt + v1.xt * v2.yz - v1.yt * v2.xz + v1.zt * v2.xy; }
	FourVector4(const TwoVector4 &v1, const Vector4 &v2, const Vector4 &v3) { initWedge(ThreeVector4(v1,v2),v3); }
	FourVector4(const Vector4 &v1, const ThreeVector4 &v2) { xyzt = v1.x * v2.yzt - v1.y * v2.xzt + v1.z * v2.xyt - v1.t * v2.xyz; }
	FourVector4(const Vector4 &v1, const TwoVector4 &v2, const Vector4 &v3) { initWedge(ThreeVector4(v1,v2),v3); }
	FourVector4(const Vector4 &v1, const Vector4 &v2, const TwoVector4 &v3) { initWedge(ThreeVector4(v1,v3),v2); }
	FourVector4(const Vector4 &v1, const Vector4 &v2, const Vector4 &v3, const Vector4 &v4) { initWedge(ThreeVector4(TwoVector4(v1,v2),v3),v4); }

	FourVector4 negation() const { return FourVector4(-xyzt); }
	FourVector4 &increase(const FourVector4 &v) { xyzt += v.xyzt; return (*this); }
	FourVector4 &decrease(const FourVector4 &v) { xyzt -= v.xyzt; return (*this); }
	FourVector4 &scale(const double d) { xyzt *= d; return (*this); }
	FourVector4 product(const FourVector4 &v) const { return FourVector4(xyzt * v.xyzt); }
	FourVector4 inverse(const double d = 1.0) const { return FourVector4(d / xyzt); }
	bool equals(const FourVector4 &v) const { return (xyzt == v.xyzt); }
	double dot(const FourVector4 &v) const { return xyzt * v.xyzt; }
	double determinant() const { return xyzt; }

	double dual() const { return xyzt; }
	double dualof() const { return xyzt; }

	SymFourMatrix4 outerProduct() const;

	double xyzt;
private:
	void initWedge(const ThreeVector4 &v1, const Vector4 &v2) { xyzt = v1.xyz * v2.t - v1.xyt * v2.z + v1.xzt * v2.y - v1.yzt * v2.x; }
};

class MatrixN;

class VectorN : public Vector<VectorN>
{
public:
	VectorN() {}
	VectorN(const uint size, const double v) : val(size, v) {}
	VectorN(const double v) {
		val.resize(1);
		val[0] = v;
	}
	VectorN(const Vector2 &v) {
		val.resize(2);
		val[0] = v.x; val[1] = v.y;
	}
	VectorN(const Vector3 &v) {
		val.resize(3);
		val[0] = v.x; val[1] = v.y; val[2] = v.z;
	}
	VectorN(const Vector4 &v) {
		val.resize(4);
		val[0] = v.x; val[1] = v.y; val[2] = v.z; val[3] = v.t;
	}
	VectorN(const Buffer<double> &vval) { val = vval; }
	VectorN(const VectorN &v) {	val = v.val; }
	VectorN &operator=(const VectorN &v) {
		val = v.val;
		return (*this);
	}

	VectorN negation() const {
		VectorN v;
		v.val.resize(size());
		for(uint i=0; i<size(); i++) v[i] = -val[i];
		return v;
	}
	VectorN &increase(const VectorN &v) {
		if(v.size() > size()) toVectorN(v.size());
		for(uint i=0; i<v.size(); i++) val[i] += v[i];
		return (*this);
	}
	VectorN &decrease(const VectorN &v) {
		if(v.size() > size()) toVectorN(v.size());
		for(uint i=0; i<v.size(); i++) val[i] -= v[i];
		return (*this);
	}
	VectorN &scale(const double d) {
		for(uint i=0; i<size(); i++) val[i] *= d;
		return (*this);
	}
	VectorN product(const VectorN &v) const {
		VectorN r;
		r.val.resize(size() < v.size() ? size() : v.size());
		for(uint i=0; i<r.size(); i++) r[i] = val[i] * v[i];
		return r;
	}
	VectorN inverse(const double d = 1.0) const {
		VectorN v;
		v.val.resize(size());
		for(uint i=0; i<size(); i++) v[i] = d / val[i];
		return v;
	}
	bool equals(const VectorN &v) const {
		if(size() != v.size()) return false;
		for(uint i=0; i<size(); i++) {
			if(val[i] != v[i]) return false;
		}
		return true;
	}
	double dot(const VectorN &v) const {
		const uint is = (size() < v.size() ? size() : v.size());
		double sum = 0.0;
		for(uint i=0; i<is; i++) sum += val[i] * v[i];
		return sum;
	}
	double determinant() const {
		double sum = 1.0;
		for(uint i=0; i<size(); i++) sum *= val[i];
		return sum;
	}
	uint size() const { return val.size(); }
	double &operator [](const uint i) const { return val[i]; }
	void toVectorN(const uint size) {
		uint i = val.size();
		val.resize(size);
		while(i<size) val[i++] = 0.0;
	}

	MatrixN outerProduct(const VectorN &v) const;

	Buffer<double> val;
};

/*TwoVector2 operator*(Vector2 &l, Vector2 &r) { return TwoVector2(l, r); }
TwoVector3 operator*(Vector3 &l, Vector3 &r) { return TwoVector3(l, r); }
TwoVector4 operator*(Vector4 &l, Vector4 &r) { return TwoVector4(l, r); }

ThreeVector3 operator*(TwoVector3 &l, Vector3 &r) { return ThreeVector3(l, r); }
ThreeVector3 operator*(Vector3 &l, TwoVector3 &r) { return ThreeVector3(l, r); }
ThreeVector4 operator*(TwoVector4 &l, Vector4 &r) { return ThreeVector4(l, r); }
ThreeVector4 operator*(Vector4 &l, TwoVector4 &r) { return ThreeVector4(l, r); }

FourVector4 operator*(Vector4 &l, ThreeVector4 &r) { return FourVector4(l, r); }
FourVector4 operator*(TwoVector4 &l, TwoVector4 &r) { return FourVector4(l, r); }
FourVector4 operator*(ThreeVector4 &l, Vector4 &r) { return FourVector4(l, r); }
*/
}

#endif //_VECTOR_HPP_INCLUDED_
