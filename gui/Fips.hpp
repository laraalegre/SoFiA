#ifndef FITS_PROCESSING_SOFTWARE_H
#define FITS_PROCESSING_SOFTWARE_H

#include <map>
#include <string>
#include <limits>

#define FITS_HEADER_UNIT_SIZE  2880
#define FITS_HEADER_ENTRY_SIZE   80
#define FITS_HEADER_KEY_SIZE      8

class Fips
{
public:
	Fips();
	Fips(const std::string &url);
	~Fips();
	
	// Member functions
	int readFile();
	int readFile(const std::string &url);
	int type();
	size_t dimension(const size_t axis = 0);
	double minimum();
	double maximum();
	double data(size_t *pos);
	std::string unit();

private:
	// Member functions
	void initializeMembers();
	
	// Helper functions
	void trimString(std::string &s);
	bool littleEndian();
	
	// Data members
	std::string  fileName;
	bool         verbose;
	int          bitpix;
	size_t       naxes;
	std::map <size_t, size_t> naxis;
	double       dataMin;
	double       dataMax;
	double       bzero;
	double       bscale;
	std::string  bunit;
	char        *dataArray;
	size_t       dataSize;
};

#endif
