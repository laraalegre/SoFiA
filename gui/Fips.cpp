// Compilation: g++ -std=c++11 -pedantic -Wall -Wextra -c Fips.cpp

#include "Fips.hpp"

#include <iostream>
#include <fstream>
#include <limits>
#include <new>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdint>



// ----------- //
// Constructor //
// ----------- //

Fips::Fips()
{
	if(verbose) std::cout << "Initialising without filename\n";
	initializeMembers();
	return;
}

Fips::Fips(const std::string &url)
{
	if(verbose) std::cout << "Initialising with filename\n";
	initializeMembers();
	fileName = url;
	
	if(not readFile())
	{
		if(verbose)
		{
			std::cout << "BITPIX   = " << bitpix  << "\n";
			std::cout << "NAXIS    = " << naxes   << "\n";
			
			for(size_t i = 1; i <= naxes; ++i)
			{
				std::map<size_t, size_t>::iterator it;
				it = naxis.find(i);
				if (it != naxis.end()) std::cout << "NAXIS" << i << "   = " << it->second << "\n";
			}
			
			std::cout << "DATAMIN  = " << dataMin << "\n";
			std::cout << "DATAMAX  = " << dataMax << "\n";
			std::cout << "BZERO    = " << bzero   << "\n";
			std::cout << "BSCALE   = " << bscale  << "\n\n";
		}
	}
	
	return;
}



// ---------- //
// Destructor //
// ---------- //

Fips::~Fips()
{
	if(dataArray != 0) delete[] dataArray;
	return;
}



// ------------------------- //
// Initialisation of members //
// ------------------------- //

void Fips::initializeMembers()
{
	verbose = true;
	
	bitpix  = 0;
	naxes   = 0;
	dataMin = std::numeric_limits<double>::quiet_NaN();
	dataMax = std::numeric_limits<double>::quiet_NaN();
	bzero   = 0.0;
	bscale  = 1.0;
	bunit.clear();
	
	dataArray = 0;
	dataSize  = 0;
	
	return;
}



// ------------------------------------- //
// Functions to return header parameters //
// ------------------------------------- //

int Fips::type()
{
	if(dataArray != 0) return bitpix;
	return 0;
}

size_t Fips::dimension(const size_t axis)
{
	if(dataArray != 0)
	{
		if(axis == 0) return naxes;
		if(axis <= naxes) return naxis[axis];
	}
	return 0;
}

double Fips::minimum()
{
	if(dataArray != 0) return dataMin;
	return std::numeric_limits<double>::quiet_NaN();
}

double Fips::maximum()
{
	if(dataArray != 0) return dataMax;
	return std::numeric_limits<double>::quiet_NaN();
}

std::string Fips::unit()
{
	if(dataArray != 0) return bunit;
	return std::string("");
}



// ----------- //
// Read header //
// ----------- //

int Fips::readFile()
{
	if(verbose) std::cout << "Attempting to read FITS file " << fileName << "\n";
	std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary);
	
	// Ensure that file is open
	if(file.is_open())
	{
		char *headerUnit = 0;
		bool endReached = false;
		
		// Ensure that EOF is not reached
		while(not file.eof() and not endReached)
		{
			headerUnit = new char[FITS_HEADER_UNIT_SIZE];
			
			// Read header unit(s)
			if(file.read(headerUnit, FITS_HEADER_UNIT_SIZE))
			{
				// Split header units into 80-byte lines
				for(size_t i = 0; i < FITS_HEADER_UNIT_SIZE / FITS_HEADER_ENTRY_SIZE; ++i)
				{
					// Extract i-th line
					std::string line(&headerUnit[i * FITS_HEADER_ENTRY_SIZE], FITS_HEADER_ENTRY_SIZE);
					
					// Extract key
					std::string key = line.substr(0, FITS_HEADER_KEY_SIZE);
					trimString(key);
					
					// Extract value
					std::string value;
					if(key == "COMMENT" or key == "HISTORY") value = line.substr(8);
					else value = line.substr(10);
					trimString(value);
					
					if(key == "BITPIX") bitpix = atoi(value.c_str());
					else if(key == "NAXIS") naxes = atoi(value.c_str());
					else if(key.substr(0, 5) == "NAXIS")
					{
						key = key.substr(5);
						size_t axis = atoi(key.c_str());
						if(axis > 0) naxis.insert(std::pair <size_t, size_t> (axis, atoi(value.c_str())));
					}
					else if(key == "DATAMIN") dataMin = atof(value.c_str());
					else if(key == "DATAMAX") dataMax = atof(value.c_str());
					else if(key == "BZERO")   bzero   = atof(value.c_str());
					else if(key == "BSCALE")  bscale  = atof(value.c_str());
					else if(key == "BUNIT")
					{
						size_t pos1 = value.find_first_of("'");
						size_t pos2 = value.find_last_of("'");
						
						if(pos1 < pos2 and pos1 != std::string::npos and pos2 != std::string::npos)
							bunit = value.substr(pos1 + 1, pos2 - pos1 - 1);
						else
							bunit.clear();
					}
					else if(key == "END")
					{
						endReached = true;
						break;
					}
				}
			}
		}
		
		
		// Clean up
		if(headerUnit) delete[] headerUnit;
		
		// Sanity checks
		if(bitpix != -64 and bitpix != -32 and bitpix != 8 and bitpix != 16 and bitpix != 32 and bitpix != 64)
		{
			std::cerr << "Error: Illegal value of BITPIX = " << bitpix << " found.\n";
			file.close();
			return 1;
		}
		
		if(naxes < 1)
		{
			std::cerr << "Error: No data unit in primary HDU of FITS file.\n";
			file.close();
			return 1;
		}
		
		if(naxes > 999)
		{
			std::cerr << "Error: FITS file claims to have more than 999 dimensions.\n";
			file.close();
			return 1;
		}
		
		for(size_t i = 1; i <= naxes; ++i)
		{
			std::map<size_t, size_t>::iterator it;
			it = naxis.find(i);
			if (it == naxis.end())
			{
				std::cerr << "Error: NAXIS" << i << " keyword missing from FITS header.\n";
				file.close();
				return 1;
			}
			else if(it->second < 1)
			{
				std::cerr << "Error: Invalid size of NAXIS" << i << " in FITS header.\n";
				file.close();
				return 1;
			}
		}
		
		if(not endReached)
		{
			std::cerr << "Error: No END keyword found in FITS header.\n";
			file.close();
			return 1;
		}
		
		// Read data, if present
		int overflowProtection = 0;
		
		dataSize = abs(bitpix / 8);
		overflowProtection += log(static_cast<double>(dataSize)) / log(2.0);
		
		for(size_t i = 1; i <= naxes; ++i)
		{
			dataSize *= naxis[i];
			overflowProtection += log(static_cast<double>(naxis[i])) / log(2.0);
		}
		
		// Integer overflow in array index detected?
		if(overflowProtection > log(std::numeric_limits<std::size_t>::max()) / log(2.0))
		{
			std::cerr << "Error: Integer overflow. Data array too large to be processed.\n";
			file.close();
			return 1;
		}
		
		// Allocate memory for data
		try
		{
			dataArray = new char[dataSize];
		}
		catch(std::bad_alloc &e)
		{
			std::cerr << "Error: Memory allocation error: " << e.what() << ".\n";
			dataArray = 0;
			file.close();
			return 1;
		}
		
		// Read in data
		if(not file.read(dataArray, dataSize))
		{
			std::cerr << "Error: Failed to read data unit of FITS file.\n";
			delete[] dataArray;
			dataArray = 0;
			file.close();
			return 1;
		}
		
		// Close file
		file.close();
	}
	else
	{
		std::cerr << "Error: Failed to open FITS file.\n";
		return 1;
	}
	
	return 0;
}


int Fips::readFile(const std::string &url)
{
	fileName = url;
	return readFile();
}



// --------------- //
// Get data values //
// --------------- //

double Fips::data(size_t *pos)
{
	if(dataArray == 0) return std::numeric_limits<double>::quiet_NaN();
	
	// Calculate index position in column-major order
	// as required by FITS standard
	size_t index = 0;
	
	for(size_t d = 1; d <= (naxes > 3 ? 3 : naxes); ++d)
	{
		size_t product = 1;
		
		if(d > 1)
			for(size_t i = 1; i < d; ++i)
				product *= naxis[i];
		
		index += product * pos[d - 1];
	}
	
	int bytes = abs(bitpix / 8);
	index *= bytes;
	
	if(index >= dataSize)
	{
		std::cerr << "Error: Data array index (" << index << ") out of range.\n";
		return std::numeric_limits<double>::quiet_NaN();
	}
	
	double result;
	
	// Copy value into result using correct data type
	// WARNING: These functions require fixed-width integers number from C++11
	//          as well as the built-in byte-swap functions of GCC!
	if(bitpix == -64 and sizeof(double) == bytes)
	{
		uint64_t tmp;
		memcpy(&tmp, &dataArray[index], bytes);
		if(littleEndian()) tmp = __builtin_bswap64(tmp);
		memcpy(&result, &tmp, bytes);
	}
	else if(bitpix == -32 and sizeof(float) == bytes)
	{
		uint32_t tmp;
		memcpy(&tmp, &dataArray[index], bytes);
		if(littleEndian()) tmp = __builtin_bswap32(tmp);
		float tmp2;
		memcpy(&tmp2, &tmp, bytes);
		result = static_cast<double>(tmp2);
	}
	else if(bitpix == 8 and sizeof(uint8_t) == bytes)
	{
		uint8_t tmp;
		memcpy(&tmp, &dataArray[index], bytes);
		result = static_cast<double>(tmp);
	}
	else if(bitpix == 16 and sizeof(uint16_t) == bytes)
	{
		uint16_t tmp;
		memcpy(&tmp, &dataArray[index], bytes);
		if(littleEndian()) tmp = __builtin_bswap16(tmp);
		result = static_cast<double>(tmp);
	}
	else if(bitpix == 32 and sizeof(uint32_t) == bytes)
	{
		uint32_t tmp;
		memcpy(&tmp, &dataArray[index], bytes);
		if(littleEndian()) tmp = __builtin_bswap32(tmp);
		result = static_cast<double>(tmp);
	}
	else if(bitpix == 64 and sizeof(uint64_t) == bytes)
	{
		uint64_t tmp;
		memcpy(&tmp, &dataArray[index], bytes);
		if(littleEndian()) tmp = __builtin_bswap64(tmp);
		result = static_cast<double>(tmp);
	}
	else
	{
		std::cerr << "Error: No native support for required data type on your system.\n";
		return std::numeric_limits<double>::quiet_NaN();
	}
	
	// Lastly, correct for BSCALE and BZERO, if necessary
	if(bzero != 0.0 or bscale != 1.0) result = bzero + bscale * result;
	
	return result;
}



// ---------------- //
// Trim std::string //
// ---------------- //

void Fips::trimString(std::string &s)
{
	// Use same whitespace characters as in PHP's trim() function
	std::string ws(" \t\n\r\0\v");
	size_t foundFirst;
	size_t foundLast;
	
	foundFirst = s.find_first_not_of(ws);
	foundLast  = s.find_last_not_of(ws);
	
	if(foundFirst != std::string::npos and foundLast != std::string::npos)
	{
		s.assign(s.substr(foundFirst, foundLast - foundFirst + 1));
	}
	else
	{
		// String contains only whitespace
		s.clear();
	}
	
	return;
}



// --------------------------- //
// Check endianness of machine //
// --------------------------- //

bool Fips::littleEndian()
{
	// Figure out endinanness of machine.
	// FITS data are big-endian according to standard,
	// so conversion to little-endian may be required.
	
	long n = 1;
	
	// Return 'true' on little-endian systems
	return *(char *)&n == 1;
}
