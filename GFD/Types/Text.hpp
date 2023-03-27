/*
Text contains rows of strings
Text is able to print, load, and save
*/

#ifndef _TEXT_HPP_INCLUDED_
#define _TEXT_HPP_INCLUDED_

#include "Types.hpp"
#include "Buffer.hpp"
#include <vector>
#include <string>
#include <sstream>

namespace gfd
{

class Text : public std::stringstream
{
public:
  Text() {}
  virtual ~Text() {}

	void clear();

	// reads output from a text file
	bool load(const std::string &path);
	// saves output into a text file
	bool save(const std::string &path);

	// is there lines anymore
	bool hasRow();
	// get the next line
	const std::string getRow();

	// convert line to numbers
	template<typename T> static Buffer<T> toBuffer(const std::string &s) {
		Buffer<T> buf;
		T d;
		std::stringstream ss(s);
		while (ss >> d) buf.push_back(d);
		return buf;
	}

};

}

#endif //_TEXT_HPP_INCLUDED_
