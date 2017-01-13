/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.0.0 (WidgetDataViewer.h) - Source Finding Application        ///
/// Copyright (C) 2016-2017 Tobias Westmeier                             ///
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

#ifndef WIDGETDATAVIEWER_H
#define WIDGETDATAVIEWER_H

#include <QtGlobal>
#include <QtCore/QVector>
#include <QtGui/QImage>
#include <QtGui/QPixmap>
#include <QtGui/QMouseEvent>
#include <QtGui/QClipboard>

#if QT_VERSION < 0x050000
	#include <QtGui/QApplication>
	#include <QtGui/QWidget>
	#include <QtGui/QLabel>
	#include <QtGui/QVBoxLayout>
	#include <QtGui/QHBoxLayout>
	#include <QtGui/QToolButton>
	#include <QtGui/QLineEdit>
	#include <QtGui/QSlider>
	#include <QtGui/QMenu>
#else
	#include <QtWidgets/QApplication>
	#include <QtWidgets/QWidget>
	#include <QtWidgets/QLabel>
	#include <QtWidgets/QVBoxLayout>
	#include <QtWidgets/QHBoxLayout>
	#include <QtWidgets/QToolButton>
	#include <QtWidgets/QLineEdit>
	#include <QtWidgets/QSlider>
	#include <QtWidgets/QMenu>
#endif

#include <string>

#include "Fips.hpp"

#define VIEWPORT_WIDTH  480
#define VIEWPORT_HEIGHT 480

#define GREYSCALE 0
#define RAINBOW   1
#define RANDOM    2

class WidgetDataViewer : public QWidget
{
	Q_OBJECT
	
public:
	WidgetDataViewer(const std::string &url, QWidget *parent = 0);
	~WidgetDataViewer();

private slots:
	void showPrevChannel();
	void showNextChannel();
	void showFirstChannel();
	void showLastChannel();
	void sliderChange(int value);
	void showContextMenu(const QPoint &where);
	void selectLutGreyscale();
	void selectLutRainbow();
	void selectLutRandom();
    void copy();

private:
	Fips   *fips;
	QImage *image;
	double  scale;
	double  offset;
	double dataMin, dataMax;
	double plotMin, plotMax;
	QVector<QRgb> lut;
	size_t currentChannel;
	
	QVBoxLayout *mainLayout;
	QHBoxLayout *layoutControls;
	
	QLabel *viewport;
	QLabel *status;
	QWidget *controls;
	
	QIcon iconGoPreviousView;
	QIcon iconGoNextView;
	QIcon iconGoFirstView;
	QIcon iconGoLastView;
	QIcon iconFillColor;
	QIcon iconDialogClose;
	QIcon iconEditCopy;
	
	QToolButton *buttonFirst;
	QToolButton *buttonLast;
	QToolButton *buttonPrev;
	QToolButton *buttonNext;
	QToolButton *buttonClose;
	QLineEdit   *fieldChannel;
	QSlider     *slider;
	
	void setUpInterface();
	void setUpLut(int type = GREYSCALE);
	int openFitsFile(const std::string &url);
	int plotChannelMap(size_t z);
	unsigned int flux2grey(double value);

protected:
	bool eventFilter(QObject *obj, QEvent *event);
	virtual void closeEvent(QCloseEvent *event);
};

#endif
