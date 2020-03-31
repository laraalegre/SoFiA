/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.3.3 (TableWidget.h) - Source Finding Application             ///
/// Copyright (C) 2014-2019 Tobias Westmeier                             ///
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

#ifndef TABLEWIDGET_H
#define TABLEWIDGET_H

#include <QtGlobal>

#include <QtCore/QString>
#include <QtCore/QList>
#include <QtCore/QDebug>

#include <QtGui/QClipboard>
#include <QtGui/QContextMenuEvent>
#include <QtGui/QKeyEvent>

// Import correct headers depending on Qt version:
#if QT_VERSION < 0x050000
	#include <QtGui/QApplication>
	#include <QtGui/QAction>
	#include <QtGui/QMenu>
	#include <QtGui/QWidget>
	#include <QtGui/QTableWidget>
#else
	#include <QtWidgets/QApplication>
	#include <QtWidgets/QAction>
	#include <QtWidgets/QMenu>
	#include <QtWidgets/QWidget>
	#include <QtWidgets/QTableWidget>
#endif

class TableWidget : public QTableWidget
{
	Q_OBJECT
	
public:
	TableWidget(QWidget *parent = 0);
	
private slots:
	void copy();
	void selectEverything();
	void selectNothing();
	void updateActions();
	
private:
	QIcon   iconCopy;
	QIcon   iconSelectAll;
	QIcon   iconClearSelection;
	
	QAction *actionCopy;
	QAction *actionSelectAll;
	QAction *actionClearSelection;
	
protected:
	void contextMenuEvent(QContextMenuEvent *event);
	void keyPressEvent(QKeyEvent *event);
};

#endif
