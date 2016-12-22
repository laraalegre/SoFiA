/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.0.0 (WidgetDataViewer.cpp) - Source Finding Application      ///
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

#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include "WidgetDataViewer.h"

// ----------- //
// CONSTRUCTOR //
// ----------- //

WidgetDataViewer::WidgetDataViewer(const std::string &url, QWidget *parent) : QWidget(parent, Qt::Window)
{
	this->setParent(parent);
	this->setWindowTitle("SoFiA Image Viewer");
	this->setAttribute(Qt::WA_DeleteOnClose);
	
	setUpInterface();
	
	scale  = 2.0;
	offset = 0.0;
	dataMin = 0.0;
	dataMax = 1.0;
	plotMin = 0.0;
	plotMax = 1.0;
	
	currentChannel = 0;
	fieldChannel->setText(QString::number(currentChannel));
	
	fips = new Fips;
	
	if(not openFitsFile(url))
	{
		plotMin = dataMin;
		plotMax = dataMax;
		plotChannelMap(currentChannel);
	}
	
	return;
}



// ---------- //
// DESTRUCTOR //
// ---------- //

WidgetDataViewer::~WidgetDataViewer()
{
	delete fips;
	return;
}



// ---------------------------- //
// FUNCTION to open a FITS file //
// ---------------------------- //

int WidgetDataViewer::openFitsFile(const std::string &url)
{
	if(fips->readFile(url)) return 1;
	
	// Redefine settings
	scale = std::max(static_cast<double>(viewport->width()) / static_cast<double>(fips->dimension(1)), static_cast<double>(viewport->height()) / static_cast<double>(fips->dimension(2)));
	
	if(std::isnan(fips->minimum()) or std::isnan(fips->maximum()) or fips->minimum() >= fips->maximum())
	{
		std::cerr << "Warning: DATAMIN or DATAMAX undefined; calculating values.\n";
		
		dataMin =  std::numeric_limits<double>::max();
		dataMax = -std::numeric_limits<double>::max();
		
		for(size_t z = 0; z < (fips->dimension() > 2 ? fips->dimension(3) : 1); ++z)
		{
			for(size_t y = 0; y < fips->dimension(2); ++y)
			{
				for(size_t x = 0; x < fips->dimension(1); ++x)
				{
					size_t pos[3] = {x, y, z};
					double value = fips->data(pos);
					
					if(not std::isnan(value))
					{
						if(value < dataMin) dataMin = value;
						if(value > dataMax) dataMax = value;
					}
				}
			}
		}
		
		if(dataMin >= dataMax)
		{
			std::cerr << "Warning: DATAMIN not greater than DATAMAX; using default values.\n";
			dataMin = 0.0;
			dataMax = 1.0;
		}
	}
	else
	{
		dataMin = fips->minimum();
		dataMax = fips->maximum();
	}
	
	//std::cout << "DATAMIN = " << dataMin << "\n";
	//std::cout << "DATAMAX = " << dataMax << "\n";
	
	// Adjust user interface
	slider->setMaximum(fips->dimension() < 3 ? 0 : fips->dimension(3) - 1);
	slider->setEnabled(fips->dimension() > 2);
	buttonFirst->setEnabled(fips->dimension() > 2);
	buttonPrev->setEnabled(fips->dimension() > 2);
	buttonNext->setEnabled(fips->dimension() > 2);
	buttonLast->setEnabled(fips->dimension() > 2);
	
	return 0;
}



// ------------------------------------------ //
// FUNCTION to draw channel map onto viewport //
// ------------------------------------------ //

int WidgetDataViewer::plotChannelMap(size_t z)
{
	if(not fips->dimension()) return 1;
	
	if(fips->dimension() < 3) z = 0;
	else if(z >= fips->dimension(3)) z = fips->dimension(3) - 1;
	
	for(size_t b = 0; b < VIEWPORT_HEIGHT; ++b)
	{
		for(size_t a = 0; a < VIEWPORT_WIDTH; ++a)
		{
			long x = floor(static_cast<double>(a) / scale - offset);
			long y = floor(static_cast<double>(viewport->height() - b - 1) / scale - offset);
			
			if(x >= 0 and static_cast<size_t>(x) < fips->dimension(1) and y >= 0 and static_cast<size_t>(y) < fips->dimension(2))
			{
				size_t position[3] = {static_cast<size_t>(x), static_cast<size_t>(y), z};
				double value = fips->data(position);
				unsigned int index = flux2grey(value);
				image->setPixel(a, b, index);
			}
			else
			{
				image->setPixel(a, b, 0);
			}
		}
	}
	
	viewport->setPixmap(QPixmap::fromImage(*image));
	
	return 0;
}



// ----------------------------- //
// SLOT to show previous channel //
// ----------------------------- //

void WidgetDataViewer::showPrevChannel()
{
	if(not fips->dimension()) return;
	if(currentChannel > 0) --currentChannel;
	fieldChannel->setText(QString::number(currentChannel));
	bool sliderBlocked = slider->blockSignals(true);
	slider->setValue(currentChannel);
	slider->blockSignals(sliderBlocked);
	plotChannelMap(currentChannel);
	return;
}



// ------------------------- //
// SLOT to show next channel //
// ------------------------- //

void WidgetDataViewer::showNextChannel()
{
	if(not fips->dimension()) return;
	if(currentChannel < fips->dimension(3) - 1) ++currentChannel;
	fieldChannel->setText(QString::number(currentChannel));
	bool sliderBlocked = slider->blockSignals(true);
	slider->setValue(currentChannel);
	slider->blockSignals(sliderBlocked);
	plotChannelMap(currentChannel);
	return;
}



// -------------------------- //
// SLOT to show first channel //
// -------------------------- //

void WidgetDataViewer::showFirstChannel()
{
	if(not fips->dimension()) return;
	currentChannel = 0;
	fieldChannel->setText(QString::number(currentChannel));
	bool sliderBlocked = slider->blockSignals(true);
	slider->setValue(currentChannel);
	slider->blockSignals(sliderBlocked);
	plotChannelMap(currentChannel);
	return;
}



// ------------------------- //
// SLOT to show last channel //
// ------------------------- //

void WidgetDataViewer::showLastChannel()
{
	if(not fips->dimension()) return;
	currentChannel = fips->dimension(3) - 1;
	fieldChannel->setText(QString::number(currentChannel));
	bool sliderBlocked = slider->blockSignals(true);
	slider->setValue(currentChannel);
	slider->blockSignals(sliderBlocked);
	plotChannelMap(currentChannel);
	return;
}



// -------------------------------- //
// SLOT to react to slider movement //
// -------------------------------- //

void WidgetDataViewer::sliderChange(int value)
{
	if(not fips->dimension()) return;
	currentChannel = value;
	fieldChannel->setText(QString::number(currentChannel));
	plotChannelMap(currentChannel);
	return;
}



// ------------------------------------------------- //
// FUNCTION to turn flux into 8-bit grey-scale value //
// ------------------------------------------------- //

unsigned int WidgetDataViewer::flux2grey(double value)
{
	if(std::isnan(value)) return 0;
	
	long grey = static_cast<int>(255.0 * (value - plotMin) / (plotMax - plotMin));
	
	if(grey < 0) grey = 0;
	else if(grey > 255) grey = 255;
	
	return static_cast<unsigned int>(grey);
}



// --------------------------------- //
// FUNCTION to set up user interface //
// --------------------------------- //

void WidgetDataViewer::setUpInterface()
{
	image = new QImage(VIEWPORT_WIDTH, VIEWPORT_HEIGHT, QImage::Format_Indexed8);
	setUpLut(RAINBOW);
	image->setColorTable(lut);
	image->fill(0);
	
	viewport = new QLabel(this);
	viewport->setPixmap(QPixmap::fromImage(*image));
	viewport->setMaximumWidth(VIEWPORT_WIDTH);
	viewport->setMinimumWidth(VIEWPORT_WIDTH);
	viewport->setMaximumHeight(VIEWPORT_HEIGHT);
	viewport->setMinimumHeight(VIEWPORT_HEIGHT);
	viewport->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	viewport->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(viewport, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(showContextMenu(const QPoint &)));
	
	status = new QLabel(this);
	status->setText("Status information");
	
	iconGoPreviousView.addFile(QString(":/icons/22/go-previous-view.png"), QSize(22, 22));
	iconGoPreviousView.addFile(QString(":/icons/16/go-previous-view.png"), QSize(16, 16));
	iconGoPreviousView  = QIcon::fromTheme("go-previous-view", iconGoPreviousView);
	
	iconGoNextView.addFile(QString(":/icons/22/go-next-view.png"), QSize(22, 22));
	iconGoNextView.addFile(QString(":/icons/16/go-next-view.png"), QSize(16, 16));
	iconGoNextView      = QIcon::fromTheme("go-next-view", iconGoNextView);
	
	iconGoFirstView.addFile(QString(":/icons/22/go-first-view.png"), QSize(22, 22));
	iconGoFirstView.addFile(QString(":/icons/16/go-first-view.png"), QSize(16, 16));
	iconGoFirstView     = QIcon::fromTheme("go-first-view", iconGoNextView);
	
	iconGoLastView.addFile(QString(":/icons/22/go-last-view.png"), QSize(22, 22));
	iconGoLastView.addFile(QString(":/icons/16/go-last-view.png"), QSize(16, 16));
	iconGoLastView      = QIcon::fromTheme("go-last-view", iconGoNextView);
	
	iconFillColor.addFile(QString(":/icons/22/fill-color.png"), QSize(22, 22));
	iconFillColor.addFile(QString(":/icons/16/fill-color.png"), QSize(16, 16));
	iconFillColor       = QIcon::fromTheme("fill-color", iconGoNextView);
	
	controls = new QWidget(this);
	buttonFirst  = new QToolButton(controls);
	buttonFirst->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonFirst->setEnabled(false);
	buttonFirst->setText("First");
	buttonFirst->setIcon(iconGoFirstView);
	buttonFirst->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	connect(buttonFirst, SIGNAL(clicked()), this, SLOT(showFirstChannel()));
	buttonLast   = new QToolButton(controls);
	buttonLast->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonLast->setEnabled(false);
	buttonLast->setText("Last");
	buttonLast->setIcon(iconGoLastView);
	buttonLast->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	connect(buttonLast, SIGNAL(clicked()), this, SLOT(showLastChannel()));
	buttonPrev   = new QToolButton(controls);
	buttonPrev->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonPrev->setEnabled(false);
	buttonPrev->setText("Prev.");
	buttonPrev->setIcon(iconGoPreviousView);
	buttonPrev->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	connect(buttonPrev, SIGNAL(clicked()), this, SLOT(showPrevChannel()));
	buttonNext   = new QToolButton(controls);
	buttonNext->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonNext->setEnabled(false);
	buttonNext->setText("Next");
	buttonNext->setIcon(iconGoNextView);
	buttonNext->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	connect(buttonNext, SIGNAL(clicked()), this, SLOT(showNextChannel()));
	fieldChannel = new QLineEdit(controls);
	fieldChannel->setMaximumWidth(50);
	fieldChannel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	fieldChannel->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	fieldChannel->setReadOnly(true);
	slider       = new QSlider(Qt::Horizontal, controls);
	slider->setEnabled(false);
	slider->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
	slider->setMinimum(0);
	slider->setMaximum(100);
	slider->setSingleStep(1);
	slider->setPageStep(10);
	connect(slider, SIGNAL(valueChanged(int)), this, SLOT(sliderChange(int)));
	
	layoutControls = new QHBoxLayout;
	layoutControls->addWidget(buttonFirst);
	layoutControls->addWidget(buttonPrev);
	layoutControls->addWidget(fieldChannel);
	layoutControls->addWidget(buttonNext);
	layoutControls->addWidget(buttonLast);
	layoutControls->addWidget(slider);
	layoutControls->setContentsMargins(0, 0, 0, 0);
	layoutControls->setSpacing(5);
	controls->setLayout(layoutControls);
	
	mainLayout = new QVBoxLayout;
	mainLayout->addWidget(viewport);
	mainLayout->addWidget(controls);
	mainLayout->addWidget(status);
	mainLayout->setContentsMargins(5, 5, 5, 5);
	mainLayout->setSpacing(5);
	this->setLayout(mainLayout);
	
	viewport->setMouseTracking(true);
	viewport->installEventFilter(this);
	
	this->show();
	
	return;
}



// --------------------------------- //
// FUNCTION to set up look-up tables //
// --------------------------------- //

void WidgetDataViewer::setUpLut(int type)
{
	lut.clear();
	srand(10);
	unsigned int valueR, valueG, valueB;
	
	for(unsigned int i = 0; i < 256; ++i)
	{
		switch(type)
		{
			case 1:
				// RAINBOW
				if (i <  50) lut.append(qRgb(0, static_cast<unsigned int>(255.0 * (static_cast<double>(i) - 0.0) / 50.0), 255));
				else if(i < 125) lut.append(qRgb(0, 255, static_cast<unsigned int>(255.0 * (125.0 - static_cast<double>(i)) / 75.0)));
				else if(i < 180) lut.append(qRgb(static_cast<unsigned int>(255 * (static_cast<double>(i) - 125) / 55), 255 ,0));
				else lut.append(qRgb(255, static_cast<unsigned int>(255 * (255 - static_cast<double>(i)) / 75), 0));
				break;
			
			case 2:
				// RANDOM
				valueR = static_cast<unsigned int>(255.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
				valueG = static_cast<unsigned int>(255.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
				valueB = static_cast<unsigned int>(255.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
				lut.append(qRgb(valueR, valueG, valueB));
				break;
			
			default:
				// GREYSCALE
				lut.append(qRgb(i, i, i));
		}
	}
	
	image->setColorTable(lut);
	return;
}

	

// ----------------------------- //
// SLOTs to change look-up table //
// ----------------------------- //

void WidgetDataViewer::selectLutGreyscale()
{
	setUpLut(GREYSCALE);
	plotChannelMap(currentChannel);
	return;
}

void WidgetDataViewer::selectLutRainbow()
{
	setUpLut(RAINBOW);
	plotChannelMap(currentChannel);
	return;
}

void WidgetDataViewer::selectLutRandom()
{
	setUpLut(RANDOM);
	plotChannelMap(currentChannel);
	return;
}



// ------------------------------------------ //
// FUNCTION to set up mouse move event filter //
// ------------------------------------------ //

bool WidgetDataViewer::eventFilter(QObject *obj, QEvent *event)
{
	if(event->type() == QEvent::MouseMove and fips->dimension())
	{
		QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
		
		int a = mouseEvent->pos().x();
		int b = mouseEvent->pos().y();
		
		long x = floor(static_cast<double>(a) / scale - offset);
		long y = floor(static_cast<double>(viewport->height() - b - 1) / scale - offset);
		
		QString text("Undefined");
		
		if(x >= 0 and static_cast<size_t>(x) < fips->dimension(1) and y >= 0 and static_cast<size_t>(y) < fips->dimension(2))
		{
			QString bunit = QString::fromStdString(fips->unit());
			size_t position[3] = {static_cast<size_t>(x), static_cast<size_t>(y), currentChannel};
			text = QString("Position: %1, %2   Value: %3   ").arg(x).arg(y).arg(fips->data(position));
			text.append(QString("Unit: "));
			text.append(bunit.isEmpty() ? QString("undefined") : bunit);
		}
		else text = QString("Position: undefined   Value: undefined");
		
		status->setText(text);
		return true;
	}
	
	return QObject::eventFilter(obj, event);
}



// --------------------------- //
// SLOT to create context menu //
// --------------------------- //

void WidgetDataViewer::showContextMenu(const QPoint &where)
{
	QMenu contextMenu("Context Menu", this);
	QMenu menuLut("Colour Scale", this);
	
	QAction actionLutGreyscale("Greyscale", this);
	QAction actionLutRainbow("Rainbow", this);
	QAction actionLutRandom("Random", this);
	/*QAction actionFirstChannel("First channel", this);
	QAction actionLastChannel("Last channel", this);
	QAction actionPrevChannel("Prev. channel", this);
	QAction actionNextChannel("Next channel", this);*/
	
	connect(&actionLutGreyscale, SIGNAL(triggered()), this, SLOT(selectLutGreyscale()));
	connect(&actionLutRainbow, SIGNAL(triggered()), this, SLOT(selectLutRainbow()));
	connect(&actionLutRandom, SIGNAL(triggered()), this, SLOT(selectLutRandom()));
	
	/*connect(&actionFirstChannel, SIGNAL(triggered()), this, SLOT(showFirstChannel()));
	connect(&actionLastChannel, SIGNAL(triggered()), this, SLOT(showLastChannel()));
	connect(&actionPrevChannel, SIGNAL(triggered()), this, SLOT(showPrevChannel()));
	connect(&actionNextChannel, SIGNAL(triggered()), this, SLOT(showNextChannel()));*/
	
	menuLut.addAction(&actionLutGreyscale);
	menuLut.addAction(&actionLutRainbow);
	menuLut.addAction(&actionLutRandom);
	menuLut.setIcon(iconFillColor);
	
	/*contextMenu.addAction(&actionPrevChannel);
	contextMenu.addAction(&actionNextChannel);
	contextMenu.addAction(&actionFirstChannel);
	contextMenu.addAction(&actionLastChannel);
	contextMenu.addSeparator();*/
	contextMenu.addMenu(&menuLut);
	
	contextMenu.exec(mapToGlobal(where));
	
	return;
}
