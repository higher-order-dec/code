/*
Buffer.hpp includes a template class Buffer, which is a simple and light implementation of a general vector.
Similar to std::vector but a bit more memory efficient in some cases.
*/

#ifndef _BUFFER_HPP_INCLUDED_
#define _BUFFER_HPP_INCLUDED_

#include <memory>
#include <stdio.h>
#include <string.h>
#include "Types.hpp"

namespace gfd
{

template <typename T>
class Buffer
{
public:
	Buffer()
	{
		m_size = 0;
		m_data = NULL;
	}
	Buffer(const uint size)
	{
		m_size = size;
		if(m_size > 0) m_data = new T[m_size];
		else m_data = NULL;
	}
	Buffer(const uint size, const T &val)
	{
		m_size = size;
		if(m_size > 0) m_data = new T[m_size];
		else m_data = NULL;
		for(uint i=0; i<m_size; i++) m_data[i] = val;
	}
	Buffer(const Buffer &buf)
	{
		m_size = buf.m_size;
		if(m_size > 0)
		{
			m_data = new T[m_size];
			memcpy(m_data, buf.m_data, m_size * sizeof(T));
		}
		else m_data = NULL;
	}
	virtual ~Buffer() { clear(); }

	Buffer &copy(const Buffer &buf, const uint size)
	{
		if(m_size == size)
		{
			if(m_size > 0) memcpy(m_data, buf.m_data, m_size * sizeof(T));
			return *this;
		}
		if(m_data) delete[] m_data;
		m_size = size;
		if(m_size > 0)
		{
			m_data = new T[m_size];
			memcpy(m_data, buf.m_data, m_size * sizeof(T));
		}
		else m_data = NULL;
		return *this;
	}
	Buffer &copy(const Buffer &buf) { return copy(buf, buf.m_size); }
	void clear()
	{
		if(m_data)
		{
			delete[] m_data;
			m_data = NULL;
			m_size = 0;
		}
	}
	void push_back(const T &val)
	{
		T *data = new T[m_size + 1];
		if(m_data)
		{
			//for(uint i=0; i<m_size; i++) data[i] = m_data[i];
			memcpy(data, m_data, m_size * sizeof(T));
			delete[] m_data;
		}
		data[m_size++] = val;
		m_data = data;
	}
	void push_front(const T &val)
	{
		T *data = new T[m_size + 1];
		if(m_data)
		{
			//for(uint i=0; i<m_size; i++) data[i] = m_data[i];
			memcpy(&data[1], m_data, m_size * sizeof(T));
			delete[] m_data;
		}
		m_size++;
		data[0] = val;
		m_data = data;
	}
	void pop_back()
	{
		if(m_size <= 1)
		{
			clear();
			return;
		}
		T *data = new T[--m_size];
		if(m_data)
		{
			//for(uint i=0; i<m_size; i++) data[i] = m_data[i];
			memcpy(data, m_data, m_size * sizeof(T));
			delete[] m_data;
		}
		m_data = data;
	}
	void pop_front()
	{
		if(m_size <= 1)
		{
			clear();
			return;
		}
		T *data = new T[--m_size];
		if(m_data)
		{
			//for(uint i=0; i<m_size; i++) data[i] = m_data[i+1];
			memcpy(data, &(m_data[1]), m_size * sizeof(T));
			delete[] m_data;
		}
		m_data = data;
	}
	void insert(const T &val, const uint i)
	{
		if(i > m_size)
		{
			resize(i + 1);
			m_data[i] = val;
		}
		else if(i == m_size)
		{
			push_back(val);
		}
		else
		{
			T *data = new T[m_size + 1];
			//for(uint j=0; j<i; j++) data[j] = m_data[j];
			if(0 < i) memcpy(data, m_data, i * sizeof(T));
			data[i] = val;
			//for(uint j=i; j<m_size; j++) data[j+1] = m_data[j];
			memcpy(&data[i+1], &m_data[i], (m_size - i) * sizeof(T));
			delete[] m_data;
			m_data = data;
			m_size++;
		}
	}
	void erase(const uint i)
	{
		if(i + 1 >= m_size) pop_back();
		else
		{
			T *data = new T[m_size - 1];
			//for(uint j=0; j<i; j++) data[j] = m_data[j];
			if(0 < i) memcpy(data, m_data, i * sizeof(T));
			//for(uint j=i+1; j<m_size; j++) data[j-1] = m_data[j];
			memcpy(&data[i], &m_data[i+1], (m_size - 1 - i) * sizeof(T));
			delete[] m_data;
			m_data = data;
			m_size--;
		}
	}
	void resize(const uint size)
	{
		if(size == m_size) return;
		if(size == 0)
		{
		  clear();
		  return;
		}
		T *data = new T[size];
		if(m_data)
		{
			//const uint minsize = m_size < size ? m_size : size;
			//for(uint i=0; i<minsize; i++) data[i] = m_data[i];
			memcpy(data, m_data, (m_size < size ? m_size : size) * sizeof(T));
			delete[] m_data;
		}
		m_data = data;
		m_size = size;
	}
	void swap(Buffer &buf)
	{
		const uint bufsize = buf.m_size;
		T *bufdata = buf.m_data;
		buf.m_size = m_size;
		buf.m_data = m_data;
		m_size = bufsize;
		m_data = bufdata;
	}
	void combine(const Buffer &buf) // insert buf at the back of this
	{
		if(buf.m_size == 0) return;
		const uint size = m_size + buf.m_size;
		T *data = new T[size];
		if(m_data)
		{
			//for(uint i=0; i<m_size; i++) data[i] = m_data[i];
			memcpy(data, m_data, m_size * sizeof(T));
			delete[] m_data;
		}
		//for(uint i=0; i<buf.m_size; i++) data[i + m_size] = buf.m_data[i];
		memcpy(&data[m_size], buf.m_data, buf.m_size * sizeof(T));
		m_data = data;
		m_size = size;
	}
	void insert(const Buffer &buf, const uint i) // insert buf at desired location
	{
		if(buf.m_size == 0) return;
		if(i > m_size) {
			resize(i + buf.size());
			memcpy(&m_data[i], buf.m_data, buf.m_size * sizeof(T));
			return;
		}
		if(i == m_size) {
			combine(buf);
			return;
		}
		T *data = new T[m_size + buf.m_size];
		if(0 < i) memcpy(data, m_data, i * sizeof(T));
		memcpy(&data[i], buf.m_data, buf.m_size * sizeof(T));
		memcpy(&data[i+buf.m_size], &m_data[i], (m_size - i) * sizeof(T));
		delete[] m_data;
		m_data = data;
		m_size += buf.m_size;
	}
	void mergesum(const Buffer &buf) { // merges another buffer by summing up each common terms
		if(m_size < buf.m_size) {
			T *data = new T[buf.m_size];
			memcpy(data, buf.m_data, buf.m_size * sizeof(T));
			if(m_data) {
				for(uint i=0; i<m_size; i++) data[i] += m_data[i];
				delete[] m_data;
			}
			m_data = data;
			m_size = buf.m_size;
		}
		else {
			for(uint i=0; i<buf.m_size; i++) m_data[i] += buf.m_data[i];
		}
	}
	bool includes(const T &val, const uint size) const
	{
		for(uint i=0; i<size; i++)
		{
			if(val == m_data[i]) return true;
		}
		return false;
	}
	bool includes(const T &val) const { return includes(val, m_size); }
	Buffer getUnion(const Buffer &buf) const
	{
		uint i, j;
		uint size = 0;
		Buffer res(m_size + buf.m_size);
		for(i=0; i<m_size; i++)
		{
			for(j=0; j<size && res[j]!=m_data[i]; j++);
			if(j == size) res[size++] = m_data[i];
		}
		for(i=0; i<buf.m_size; i++)
		{
			for(j=0; j<size && res[j]!=buf.m_data[i]; j++);
			if(j == size) res[size++] = buf.m_data[i];
		}
		res.resize(size);
		return res;
	}
	Buffer getIntersection(const Buffer &buf) const
	{
		uint size = 0;
		Buffer res;
		for(uint i=0; i<m_size; i++)
		{
			for(uint j=0; j<buf.m_size; j++)
			{
				if(m_data[i] == buf.m_data[j])
				{
					res.gatherOnce(m_data[i], size);
					break;
				}
			}
		}
		res.resize(size);
		return res;
	}
	const T &getFirstIntersection(const Buffer &buf, const T &none) const
	{
		for(uint i=0; i<m_size; i++)
		{
			for(uint j=0; j<buf.m_size; j++)
			{
				if(m_data[i] == buf.m_data[j]) return m_data[i];
			}
		}
		return none;
	}
	Buffer getComplement(const Buffer &buf) const
	{
		uint size = 0;
		Buffer res(m_size);
		for(uint i=0; i<m_size; i++)
		{
			uint j = 0;
			while(j < buf.m_size && m_data[i] != buf.m_data[j]) j++;
			if(j == buf.m_size) res[size++] = m_data[i];
		}
		res.resize(size);
		return res;
	}
	const T &getFirstComplement(const Buffer &buf, const T &none) const
	{
		for(uint i=0; i<m_size; i++)
		{
			uint j = 0;
			while(j < buf.m_size && m_data[i] != buf.m_data[j]) j++;
			if(j == buf.m_size) return m_data[i];
		}
		return none;
	}
	bool isAnagram(const Buffer &buf) const // return true, if both buffers have exactly the same values
	{
		if(buf.m_size != m_size) return false;
		if(m_size == 0) return true;
		uint i, j;
		bool *chk = new bool[m_size];
		for(i=0; i<m_size; i++) chk[i] = false;
		for(i=0; i<m_size; i++)
		{
			const T &val = buf.m_data[i];
			j = 0;
			while(chk[j] || val != m_data[j])
			{
				if(++j >= m_size)
				{
					delete[] chk;
					return false;
				}
			}
			chk[j] = true;
		}
		delete[] chk;
		return true;
	}
	void gather(const T &val, uint &size) // remember to resize buffer with size after all gathered
	{
		if(size >= m_size) resize(2 * size + 1);
		m_data[size++] = val;
	}
	bool gatherOnce(const T &val, uint &size) // remember to resize buffer with size after all gathered
	{
		for(uint i=0; i<size; i++)
		{
			if(val == m_data[i]) return false;
		}
		gather(val, size);
		return true;
	}
	bool ungather(const T &val, uint &size)
	{
		for(uint i=0; i<size; i++)
		{
			if(val == m_data[i])
			{
				for(--size; i<size; i++) m_data[i] = m_data[i+1];
				return true;
			}
		}
		return false;
	}
	bool gatherOrUngather(const T &val, uint &size)
	{
		if(ungather(val, size)) return false;
		gather(val, size);
		return true;
	}
	uint findFirst(const T &val, const uint size) const // returns the index of the first instance of the value
	{
		for(uint i=0; i<size; i++)
		{
			if(m_data[i] == val) return i;
		}
		return size;
	}
	uint findFirst(const T &val) const { return findFirst(val, m_size); }
	bool eraseFirst(const T &val) // erase first instance of the value
	{
		for(uint i=0; i<m_size; i++)
		{
			if(val == m_data[i])
			{
				m_data[i] = m_data[m_size - 1];
				pop_back();
				return true;
			}
		}
		return false;
	}
	bool replaceFirst(const T &oldval, const T &newval) // replace first instance of the value with new value
	{
		for(uint i=0; i<m_size; i++)
		{
			if(oldval == m_data[i])
			{
				m_data[i] = newval;
				return true;
			}
		}
		return false;
	}
	void fill(const T &val)
	{
		for(uint i=0; i<m_size; i++) m_data[i] = val;
	}
	T &operator [](const uint i) const { return m_data[i]; }
	bool operator==(const Buffer &buf) const {
		if(m_size != buf.m_size) return false;
		for(uint i=0; i<m_size; i++) {
			if(m_data[i] != buf.m_data[i]) return false;
		}
		return true;
	}
	bool operator!=(const Buffer &buf) const { return !(*this == buf); }
	Buffer &operator=(const Buffer &buf) { return copy(buf, buf.m_size); }
	uint size() const { return m_size; }
	bool empty() const { return (m_size == 0); }
	T &back() const { return m_data[m_size - 1]; }
	T &front() const { return m_data[0]; }

protected:
  uint m_size;
  T *m_data;
};

}

#endif //_BUFFER_HPP_INCLUDED_
