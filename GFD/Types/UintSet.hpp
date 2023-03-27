/*
UintSet.hpp includes a class, which can be used to collect sections of unsigned integers.
By calling function .includes(uint i) one tests if i includes in the UintSet.
*/

#ifndef _UINTSET_HPP_INCLUDED_
#define _UINTSET_HPP_INCLUDED_

#include "Types.hpp"
#include "Buffer.hpp"

namespace gfd
{

class UintSet
{
public:
	UintSet() {}
	UintSet(const uint i) { insert(i, i, 0); }
	UintSet(const uint i0, const uint i1, const uint ip = 0) { insert(i0, i1, ip); }
	UintSet(const UintSet &us) { m_sec = us.m_sec; }
	UintSet &operator = (const UintSet &us)
	{
		m_sec = us.m_sec;
		return *this;
	}

	void insert(const uint i) { insert(i, i, 0); }
	void insert(const uint i0, const uint i1, const uint ip = 0)
	{
		m_sec.push_back(UintSection(i0, i1, ip));
	}
	void insert(const UintSet &us)
	{
		m_sec.combine(us.m_sec);
	}
	bool includes(const uint i) const
	{
		for(uint j=0; j<m_sec.size(); j++)
		{
			const UintSection &sec = m_sec[j];
			if(i < sec.i0) continue;
			const uint ij = (sec.ip == 0 ? i : ((i - sec.i0) % sec.ip) + sec.i0);
			if(ij <= sec.i1) return true;
		}
		return false;
	}

protected:
	struct UintSection
	{
		uint i0; // first id
		uint i1; // last id
		uint ip; // id period (0 = non-periodic section)
		UintSection(const uint ai0 = 0, const uint ai1 = 0, const uint aip = 0)
		{
			i0 = ai0;
			i1 = ai1;
			ip = aip;
		}
	};

	Buffer<UintSection> m_sec;
};

const UintSet UINTSETALL(0, uint(-1), 0);

}

#endif //_UINTSET_HPP_INCLUDED_
