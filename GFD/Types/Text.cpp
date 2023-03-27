#include "Text.hpp"
#include <fstream>
#include <cstdlib>

using namespace gfd;

static const uint MAXLINEWIDTH = 4096;

void Text::clear()
{
	str("");
}

bool Text::load(const std::string &path)
{
	clear();

	std::ifstream fs(path.c_str(), std::ios::in);
	if(fs.fail() != 0) return false;

	// get size of the file
	fs.seekg(0, std::ios::end);
	const uint size = uint(fs.tellg());
	fs.seekg(0, std::ios::beg);

	if(size > 0) {
		// read the content of *this
		char *buffer = new char[size];
		fs.read(buffer, size);

		// write the content to the file stream
		write(buffer, size);
		delete[] buffer;
	}

	// close the file and exit
	fs.close();
	return true;
}

bool Text::save(const std::string &path)
{
	std::ofstream fs(path.c_str(), std::ios_base::trunc);
	if(fs.fail() != 0) return false;

	// get size of *this
	seekg(0, std::ios::end);
	const uint size = uint(tellg());
	seekg(0, std::ios::beg);

	if(size > 0) {
		// read the content of *this
		char *buffer = new char[size];
		read(buffer, size);

		// write the content to the file stream
		fs.write(buffer, size);
		delete[] buffer;
	}

	// close the file and exit
	fs.close();
	return true;
}

bool Text::hasRow()
{
	if(eof())
	{
		clear();
		return false;
	}
	return true;
}

const std::string Text::getRow()
{
	char str[MAXLINEWIDTH];
	getline(str, MAXLINEWIDTH);
	return str;
}
/*
Buffer<double> Text::toDoubles(const std::string &s) {
	Buffer<double> buf;
	double d = 0.0;
	std::stringstream ss(s);
    while (ss >> d) buf.push_back(d);
	return buf;
}
*/