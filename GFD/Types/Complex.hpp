/*
Complex.hpp implements complex numbers and vectors.
A complex number a+bi can be considered as 2x2-matrix:
a -b
b  a
Thus, the transpose() returns conjugate transpose and determinant() is a product of lensq()s.
The complex number class is called Complex.
The complex vector classes are called Complex2, Complex3, and Complex4.
*/

#ifndef _COMPLEX_HPP_INCLUDED_
#define _COMPLEX_HPP_INCLUDED_

#include "Types.hpp"
#include "Vector.hpp"

namespace gfd
{


class Complex : public Vector<Complex>
{
public:
	Complex() {}
	Complex(const double cx, const double cy) {	r = cx;	i = cy;	}
	Complex(const Complex &c) {	r = c.r; i = c.i; }

	Complex transpose() const { return Complex(r, -i); }
	Complex negation() const { return Complex(-r, -i); }
	Complex &increase(const Complex &c) { r += c.r; i += c.i; return (*this); }
	Complex &decrease(const Complex &c) { r -= c.r; i -= c.i; return (*this); }
	Complex &scale(const double d) { r *= d; i *= d; return (*this); }
	Complex product(const Complex &c) const { return Complex(r * c.r - i * c.i, i * c.r + r * c.i); }
	Complex inverse(const double d = 1.0) const { const double h = d / lensq(); return Complex(r * h, -i * h); }
	bool equals(const Complex &c) const { return (r == c.r) && (i == c.i); }
	double dot(const Complex &c) const { return r * c.r + i * c.i; }
	double determinant() const { return lensq(); }

	Complex times_i() const { return Complex(-i, r); }
	Complex per_i() const { return Complex(i, -r); }
	Complex sqroot() const
	{
		const double h = len();
		return Complex(h > -r ? std::sqrt(0.5 * (h + r)) : 0.0, h > r ? std::sqrt(0.5 * (h - r)) : 0.0);
	}
	Complex exponential() const
	{
		const double h = std::exp(r);
		return Complex(h * std::cos(i), h * std::sin(i));
	}

	double r;
	double i;
};

class Complex2 : public Vector<Complex2>
{
public:
	Complex2() {}
	Complex2(const double vxr, const double vxi, const double vyr, const double vyi) : x(vxr,vxi), y(vyr,vyi) {}
	Complex2(const Complex &vx, const Complex &vy) : x(vx), y(vy) {}
	Complex2(const Complex2 &v) : x(v.x), y(v.y) {}

	Complex2 negation() const { return Complex2(-x, -y); }
	Complex2 &increase(const Complex2 &v) { x += v.x; y += v.y; return (*this); }
	Complex2 &decrease(const Complex2 &v) { x -= v.x; y -= v.y; return (*this); }
	Complex2 &scale(const double d) { x *= d; y *= d; return (*this); }
	Complex2 product(const Complex2 &v) const { return Complex2(x * v.x, y * v.y); }
	Complex2 inverse(const double d = 1.0) const { return Complex2(d / x, d / y); }
	bool equals(const Complex2 &v) const { return (x == v.x) && (y == v.y); }
	double dot(const Complex2 &v) const { return x.dot(v.x) + y.dot(v.y); }
	double determinant() const { return x.lensq() * y.lensq(); }

	Complex x;
	Complex y;
};

class Complex3 : public Vector<Complex3>
{
public:
	Complex3() {}
	Complex3(const double vxr, const double vxi, const double vyr, const double vyi, const double vzr, const double vzi) : x(vxr,vxi), y(vyr,vyi), z(vzr,vzi) {}
	Complex3(const Complex &vx, const Complex &vy, const Complex &vz) : x(vx), y(vy), z(vz) {}
	Complex3(const Complex2 &v, const double vzr, const double vzi) : x(v.x), y(v.y), z(vzr,vzi) {}
	Complex3(const Complex3 &v) : x(v.x), y(v.y), z(v.z) {}

	Complex3 negation() const { return Complex3(-x, -y, -z); }
	Complex3 &increase(const Complex3 &v) { x += v.x; y += v.y; z += v.z; return (*this); }
	Complex3 &decrease(const Complex3 &v) { x -= v.x; y -= v.y; z -= v.z; return (*this); }
	Complex3 &scale(const double d) { x *= d; y *= d; z *= d; return (*this); }
	Complex3 product(const Complex3 &v) const { return Complex3(x * v.x, y * v.y, z * v.z); }
	Complex3 inverse(const double d = 1.0) const { return Complex3(d / x, d / y, d / z); }
	bool equals(const Complex3 &v) const { return (x == v.x) && (y == v.y) && (z == v.z); }
	double dot(const Complex3 &v) const { return x.dot(v.x) + y.dot(v.y) + z.dot(v.z); }
	double determinant() const { return x.lensq() * y.lensq() * z.lensq(); }

	Complex2 toComplex2() const { return Complex2(x, y); }

	Complex x;
	Complex y;
	Complex z;
};

class Complex4 : public Vector<Complex4>
{
public:
	Complex4() {}
	Complex4(const double vxr, const double vxi, const double vyr, const double vyi, const double vzr, const double vzi, const double vtr, const double vti)
	: x(vxr,vxi), y(vyr,vyi), z(vzr,vzi), t(vtr,vti) {}
	Complex4(const Complex &vx, const Complex &vy, const Complex &vz, const Complex &vt) : x(vx), y(vy), z(vz), t(vt) {}
	Complex4(const Complex2 &v, const double vzr, const double vzi, const double vtr, const double vti) : x(v.x), y(v.y), z(vzr,vzi), t(vtr,vti) {}
	Complex4(const Complex3 &v, const double vtr, const double vti) : x(v.x), y(v.y), z(v.z), t(vtr,vti) {}
	Complex4(const Complex4 &v) : x(v.x), y(v.y), z(v.z), t(v.t) {}

	Complex4 negation() const { return Complex4(-x, -y, -z, -t); }
	Complex4 &increase(const Complex4 &v) { x += v.x; y += v.y; z += v.z; t += v.t; return (*this); }
	Complex4 &decrease(const Complex4 &v) { x -= v.x; y -= v.y; z -= v.z; t -= v.t; return (*this); }
	Complex4 &scale(const double d) { x *= d; y *= d; z *= d; t *= d; return (*this); }
	Complex4 product(const Complex4 &v) const { return Complex4(x * v.x, y * v.y, z * v.z, t * v.t); }
	Complex4 inverse(const double d = 1.0) const { return Complex4(d / x, d / y, d / z, d / t); }
	bool equals(const Complex4 &v) const { return (x == v.x) && (y == v.y) && (z == v.z) && (t == v.t); }
	double dot(const Complex4 &v) const { return x.dot(v.x) + y.dot(v.y) + z.dot(v.z) + t.dot(v.t); }
	double determinant() const { return x.lensq() * y.lensq() * z.lensq() * t.lensq(); }

	Complex2 toComplex2() const { return Complex2(x, y); }
	Complex3 toComplex3() const { return Complex3(x, y, z); }

	Complex x;
	Complex y;
	Complex z;
	Complex t;
};

}

#endif //_COMPLEX_HPP_INCLUDED_
