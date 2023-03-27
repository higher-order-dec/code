/*
Types.hpp introduces most of the necessary data types, and constants to be used in the namespace gfd.
*/

#ifndef _TYPES_HPP_INCLUDED_
#define _TYPES_HPP_INCLUDED_

#include <cmath>

namespace gfd
{


template<typename T> class Couple {
public:
	Couple() { }
	Couple(const T &aa, const T &bb) { a = aa; b = bb; }
	Couple &operator+=(const Couple &c) { a += c.a; b += c.b; return (*this); }
	Couple operator+(const Couple &c) const { return Couple(a + c.a, b + c.b); }
	Couple &operator-=(const Couple &c) { a -= c.a; b -= c.b; return (*this); }
	Couple operator-(const Couple &c) const { return Couple(a - c.a, b - c.b); }
	bool operator==(const Couple &c) { return a == c.a && b == c.b; }
	bool operator<(const Couple &c) { return a < c.a || (a == c.a && b < c.b); }
	bool operator>(const Couple &c) { return a > c.a || (a == c.a && b > c.b); }
	T a;
	T b;
};


typedef char sign;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef long long llong;
typedef unsigned long long ullong;

const uint NONE = uint(-1);
const double PI = 3.1415926535897932385;
const double PIx2 = 6.2831853071795864770;

}

#endif //_TYPES_HPP_INCLUDED_
