/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.2.0 (main.cpp) - Source Finding Application                  ///
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

#include <QtGlobal>

// Import correct headers depending on Qt version:
#if QT_VERSION < 0x050000
	#include <QtGui/QApplication>
#else
	#include <QtWidgets/QApplication>
#endif

#include "SoFiA.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	
	SoFiA pipelineInterface(argc, argv);
	pipelineInterface.show();
	
	return app.exec();
}
