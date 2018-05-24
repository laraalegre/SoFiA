/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.2.0 (Fips.hpp) - Source Finding Application                  ///
/// Copyright (C) 2013-2018 Tobias Westmeier                             ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// Address:  Tobias Westmeier                                           ///
///           ICRAR M468                                                 ///
///           The University of Western Australia                        ///
///           35 Stirling Highway                                        ///
///           Crawley WA 6009                                            ///
///           Australia                                                  ///
///                                                                      ///
/// E-mail:   tobias.westmeier [at] uwa.edu.au                           ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// This program is free software: you can redistribute it and/or modify ///
/// it under the terms of the GNU General Public License as published by ///
/// the Free Software Foundation, either version 3 of the License, or    ///
/// (at your option) any later version.                                  ///
///                                                                      ///
/// This program is distributed in the hope that it will be useful,      ///
/// but WITHOUT ANY WARRANTY; without even the implied warranty of       ///
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ///
/// GNU General Public License for more details.                         ///
///                                                                      ///
/// You should have received a copy of the GNU General Public License    ///
/// along with this program. If not, see http://www.gnu.org/licenses/.   ///
/// ____________________________________________________________________ ///
///                                                                      ///

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
