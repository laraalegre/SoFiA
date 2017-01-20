/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.1.0-beta (WidgetDataViewer.cpp) - Source Finding Application      ///
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
#include <ctime>
#include "WidgetDataViewer.h"

// ----------- //
// CONSTRUCTOR //
// ----------- //

WidgetDataViewer::WidgetDataViewer(const std::string &url, QWidget *parent) : QWidget(parent, Qt::Window)
{
	this->setParent(parent);
	this->setWindowTitle("SoFiA - Image Viewer");
	this->setAttribute(Qt::WA_DeleteOnClose);
	
	// Initialise random number generator
	srand(time(0));
	
	scale   = 2.0;
	offsetX = 0.0;
	offsetY = 0.0;
	dataMin = 0.0;
	dataMax = 1.0;
	plotMin = 0.0;
	plotMax = 1.0;
	
	revert = 0x00;
	invert = 0x00;
	currentLut = RAINBOW;
	transferFunction = LINEAR;
	currentChannel = 0;
	
	setUpInterface();
	fieldChannel->setText(QString::number(currentChannel));
	
	fips = new Fips;
	
	if(not openFitsFile(url))
	{
		plotMin = dataMin;
		plotMax = dataMax;
		plotChannelMap(currentChannel);
	}
	
	fieldLevelMin->setText(QString::number(plotMin));
	fieldLevelMax->setText(QString::number(plotMax));
	
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
	zoomToFit();
	
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
			long x = floor(static_cast<double>(a) / scale - offsetX);
			long y = floor(static_cast<double>(VIEWPORT_HEIGHT - b - 1) / scale - offsetY);
			
			if(x >= 0 and static_cast<size_t>(x) < fips->dimension(1) and y >= 0 and static_cast<size_t>(y) < fips->dimension(2))
			{
				size_t position[3] = {static_cast<size_t>(x), static_cast<size_t>(y), z};
				double value = fips->data(position);
				unsigned int index = flux2grey(value);
				
				image->setPixel(a, b, index xor revert);
			}
			else
			{
				image->setPixel(a, b, (((a + b) % 2) * 30) xor revert);
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
	if(std::isnan(value) or value < plotMin) return 0;
	else if(value > plotMax) return 255;
	
	switch(transferFunction)
	{
		case SQRT:
			return static_cast<unsigned int>(255.0 * sqrt((value - plotMin) / (plotMax - plotMin)));
		case LOG:
			return static_cast<unsigned int>(255.0 * log10(9.0 * (value - plotMin) / (plotMax - plotMin) + 1.0));
	}
	
	// Default for linear transfer function
	return static_cast<unsigned int>(255.0 * (value - plotMin) / (plotMax - plotMin));
}



// --------------------------------- //
// FUNCTION to set up user interface //
// --------------------------------- //

void WidgetDataViewer::setUpInterface()
{
	// Set up icons
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
	
	iconDialogClose.addFile(QString(":/icons/22/dialog-close.png"), QSize(22, 22));
	iconDialogClose.addFile(QString(":/icons/16/dialog-close.png"), QSize(16, 16));
	iconDialogClose = QIcon::fromTheme("dialog-close", iconDialogClose);
	
	iconEditCopy.addFile(QString(":/icons/22/edit-copy.png"), QSize(22, 22));
	iconEditCopy.addFile(QString(":/icons/16/edit-copy.png"), QSize(16, 16));
	iconEditCopy = QIcon::fromTheme("edit-copy", iconEditCopy);
	
	iconEditReset.addFile(QString(":/icons/22/edit-reset.png"), QSize(22, 22));
	iconEditReset.addFile(QString(":/icons/16/edit-reset.png"), QSize(16, 16));
	//iconEditReset = QIcon::fromTheme("edit-reset", iconEditReset);
	
	iconZoomIn.addFile(QString(":/icons/22/zoom-in.png"), QSize(22, 22));
	iconZoomIn.addFile(QString(":/icons/16/zoom-in.png"), QSize(16, 16));
	iconZoomIn = QIcon::fromTheme("zoom-in", iconEditCopy);
	
	iconZoomOut.addFile(QString(":/icons/22/zoom-out.png"), QSize(22, 22));
	iconZoomOut.addFile(QString(":/icons/16/zoom-out.png"), QSize(16, 16));
	iconZoomOut = QIcon::fromTheme("zoom-out", iconEditCopy);
	
	iconZoomFitBest.addFile(QString(":/icons/22/zoom-fit-best.png"), QSize(22, 22));
	iconZoomFitBest.addFile(QString(":/icons/16/zoom-fit-best.png"), QSize(16, 16));
	iconZoomFitBest = QIcon::fromTheme("zoom-fit-best", iconEditCopy);
	
	// Set up actions
	actionGroupColourScheme = new QActionGroup(this);
	actionGroupColourScheme->setExclusive(true);
	actionLutGreyscale = new QAction("Greyscale", actionGroupColourScheme);
	actionLutGreyscale->setCheckable(true);
	actionLutGreyscale->setShortcut(Qt::Key_1);
	connect(actionLutGreyscale, SIGNAL(triggered()), this, SLOT(selectLutGreyscale()));
	actionLutRainbow = new QAction("Rainbow", actionGroupColourScheme);
	actionLutRainbow->setCheckable(true);
	actionLutRainbow->setChecked(true);
	actionLutRainbow->setShortcut(Qt::Key_2);
	connect(actionLutRainbow, SIGNAL(triggered()), this, SLOT(selectLutRainbow()));
	actionLutRgb = new QAction("Velocity", actionGroupColourScheme);
	actionLutRgb->setCheckable(true);
	actionLutRgb->setShortcut(Qt::Key_3);
	connect(actionLutRgb, SIGNAL(triggered()), this, SLOT(selectLutRgb()));
	actionLutRandom = new QAction("Random", actionGroupColourScheme);
	actionLutRandom->setCheckable(true);
	actionLutRandom->setShortcut(Qt::Key_4);
	connect(actionLutRandom, SIGNAL(triggered()), this, SLOT(selectLutRandom()));
	
	actionRevert = new QAction("Revert", this);
	actionRevert->setToolTip("Revert colour scheme");
	actionRevert->setCheckable(true);
	actionRevert->setShortcut(Qt::Key_R);
	connect(actionRevert, SIGNAL(triggered()), this, SLOT(toggleRev()));
	actionInvert = new QAction("Invert", this);
	actionInvert->setToolTip("Invert colour scheme");
	actionInvert->setCheckable(true);
	actionInvert->setShortcut(Qt::Key_I);
	connect(actionInvert, SIGNAL(triggered()), this, SLOT(toggleInv()));
	
	actionFirst = new QAction("First", this);
	actionFirst->setToolTip("First channel");
	actionFirst->setIcon(iconGoFirstView);
	actionFirst->setShortcut(Qt::Key_Home);
	connect(actionFirst, SIGNAL(triggered()), this, SLOT(showFirstChannel()));
	actionPrev = new QAction("Previous", this);
	actionPrev->setToolTip("Previous channel");
	actionPrev->setIcon(iconGoPreviousView);
	actionPrev->setShortcut(Qt::Key_PageUp);
	connect(actionPrev, SIGNAL(triggered()), this, SLOT(showPrevChannel()));
	actionNext = new QAction("Next", this);
	actionNext->setToolTip("Next channel");
	actionNext->setIcon(iconGoNextView);
	actionNext->setShortcut(Qt::Key_PageDown);
	connect(actionNext, SIGNAL(triggered()), this, SLOT(showNextChannel()));
	actionLast = new QAction("Last", this);
	actionLast->setToolTip("Last channel");
	actionLast->setIcon(iconGoLastView);
	actionLast->setShortcut(Qt::Key_End);
	connect(actionLast, SIGNAL(triggered()), this, SLOT(showLastChannel()));
	
	actionZoomIn = new QAction("Zoom In", this);
	actionZoomIn->setToolTip("Zoom in");
	actionZoomIn->setIcon(iconZoomIn);
	actionZoomIn->setShortcut(Qt::Key_Plus);
	connect(actionZoomIn, SIGNAL(triggered()), this, SLOT(zoomIn()));
	actionZoomOut = new QAction("Zoom Out", this);
	actionZoomOut->setToolTip("Zoom out");
	actionZoomOut->setIcon(iconZoomOut);
	actionZoomOut->setShortcut(Qt::Key_Minus);
	connect(actionZoomOut, SIGNAL(triggered()), this, SLOT(zoomOut()));
	actionZoomToFit = new QAction("Zoom To Fit", this);
	actionZoomToFit->setToolTip("Zoom to fit");
	actionZoomToFit->setIcon(iconZoomFitBest);
	actionZoomToFit->setShortcut(Qt::Key_F);
	connect(actionZoomToFit, SIGNAL(triggered()), this, SLOT(zoomToFit()));
	
	
	actionCopy = new QAction("Copy", this);
	actionCopy->setIcon(iconEditCopy);
	actionCopy->setShortcut(QKeySequence::Copy);
	connect(actionCopy, SIGNAL(triggered()), this, SLOT(copy()));
	
	// The following is necessary to activate
	// the action's shortcut key bindings:
	this->addAction(actionCopy);
	this->addAction(actionLutGreyscale);
	this->addAction(actionLutRainbow);
	this->addAction(actionLutRgb);
	this->addAction(actionLutRandom);
	this->addAction(actionRevert);
	this->addAction(actionInvert);
	
	// Set up data display image
	image = new QImage(VIEWPORT_WIDTH, VIEWPORT_HEIGHT, QImage::Format_Indexed8);
	setUpLut(RAINBOW);
	image->setColorTable(lut);
	image->fill(0);
	
	// Set up image viewport:
	viewport = new QLabel(this);
	viewport->setFocusPolicy(Qt::StrongFocus);
	viewport->setPixmap(QPixmap::fromImage(*image));
	viewport->setMaximumWidth(VIEWPORT_WIDTH);
	viewport->setMinimumWidth(VIEWPORT_WIDTH);
	viewport->setMaximumHeight(VIEWPORT_HEIGHT);
	viewport->setMinimumHeight(VIEWPORT_HEIGHT);
	viewport->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	viewport->setContextMenuPolicy(Qt::CustomContextMenu);
	viewport->setFrameStyle(QFrame::Panel | QFrame::Sunken);
	viewport->setLineWidth(1);
	connect(viewport, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(showContextMenu(const QPoint &)));
	
	// Set up settings widget
	settings = new QWidget(this);
	labelLevelMin = new QLabel(settings);
	labelLevelMin->setText("Range:");
	labelLevelMin->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	labelLevelMax = new QLabel(settings);
	labelLevelMax->setText(QChar(0x2013));
	labelLevelMax->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	fieldLevelMin = new QLineEdit(settings);
	fieldLevelMin->setMaximumWidth(200);
	fieldLevelMin->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
	fieldLevelMin->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	fieldLevelMin->setToolTip("Lower intensity threshold");
	connect(fieldLevelMin, SIGNAL(editingFinished()), this, SLOT(setLevelMin()));
	fieldLevelMax = new QLineEdit(settings);
	fieldLevelMax->setMaximumWidth(200);
	fieldLevelMax->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
	fieldLevelMax->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	fieldLevelMax->setToolTip("Upper intensity threshold");
	connect(fieldLevelMax, SIGNAL(editingFinished()), this, SLOT(setLevelMax()));
	fieldTransFunc = new QComboBox(settings);
	fieldTransFunc->insertItem(LINEAR, "linear");
	fieldTransFunc->insertItem(SQRT, "sqrt");
	fieldTransFunc->insertItem(LOG, "log");
	fieldTransFunc->setToolTip("Transfer function");
	connect(fieldTransFunc, SIGNAL(currentIndexChanged(int)), this, SLOT(setTransferFunction(int)));
	
	buttonZoomIn = new QToolButton(settings);
	buttonZoomIn->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonZoomIn->setDefaultAction(actionZoomIn);
	buttonZoomIn->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	buttonZoomToFit = new QToolButton(settings);
	buttonZoomToFit->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonZoomToFit->setDefaultAction(actionZoomToFit);
	buttonZoomToFit->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	buttonZoomOut = new QToolButton(settings);
	buttonZoomOut->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonZoomOut->setDefaultAction(actionZoomOut);
	buttonZoomOut->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	
	buttonReset = new QToolButton(settings);
	buttonReset->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonReset->setText("Reset");
	buttonReset->setToolTip("Reset display");
	buttonReset->setIcon(iconEditReset);
	buttonReset->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	connect(buttonReset, SIGNAL(clicked()), this, SLOT(resetDisplaySettings()));
	
	layoutSettings = new QHBoxLayout;
	layoutSettings->addWidget(labelLevelMin);
	layoutSettings->addWidget(fieldLevelMin);
	layoutSettings->addWidget(labelLevelMax);
	layoutSettings->addWidget(fieldLevelMax);
	layoutSettings->addWidget(fieldTransFunc);
	layoutSettings->addStretch();
	layoutSettings->addWidget(buttonZoomIn);
	layoutSettings->addWidget(buttonZoomToFit);
	layoutSettings->addWidget(buttonZoomOut);
	layoutSettings->addStretch();
	layoutSettings->addWidget(buttonReset);
	layoutSettings->setContentsMargins(0, 0, 0, 0);
	layoutSettings->setSpacing(5);
	settings->setLayout(layoutSettings);
	
	// Set up status widget
	status = new QLabel(this);
	status->setText("Status information");
	
	// Set up control widget
	controls = new QWidget(this);
	buttonFirst = new QToolButton(controls);
	buttonFirst->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonFirst->setEnabled(false);
	buttonFirst->setDefaultAction(actionFirst);
	buttonFirst->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	buttonPrev = new QToolButton(controls);
	buttonPrev->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonPrev->setEnabled(false);
	buttonPrev->setDefaultAction(actionPrev);
	buttonPrev->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	fieldChannel = new QLineEdit(controls);
	fieldChannel->setMaximumWidth(50);
	fieldChannel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	fieldChannel->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	fieldChannel->setReadOnly(true);
	buttonNext = new QToolButton(controls);
	buttonNext->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonNext->setEnabled(false);
	buttonNext->setDefaultAction(actionNext);
	buttonNext->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	buttonLast = new QToolButton(controls);
	buttonLast->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonLast->setEnabled(false);
	buttonLast->setDefaultAction(actionLast);
	buttonLast->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	slider = new QSlider(Qt::Horizontal, controls);
	slider->setEnabled(false);
	slider->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
	slider->setMinimum(0);
	slider->setMaximum(100);
	slider->setSingleStep(1);
	slider->setPageStep(10);
	connect(slider, SIGNAL(valueChanged(int)), this, SLOT(sliderChange(int)));
	buttonClose = new QToolButton(controls);
	buttonClose->setToolButtonStyle(Qt::ToolButtonIconOnly);
	buttonClose->setEnabled(true);
	buttonClose->setText("Close");
	buttonClose->setToolTip("Close viewer");
	buttonClose->setIcon(iconDialogClose);
	buttonClose->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	connect(buttonClose, SIGNAL(clicked()), this, SLOT(close()));
	
	layoutControls = new QHBoxLayout;
	layoutControls->addWidget(buttonFirst);
	layoutControls->addWidget(buttonPrev);
	layoutControls->addWidget(fieldChannel);
	layoutControls->addWidget(buttonNext);
	layoutControls->addWidget(buttonLast);
	layoutControls->addWidget(slider);
	layoutControls->addWidget(buttonClose);
	layoutControls->setContentsMargins(0, 0, 0, 0);
	layoutControls->setSpacing(5);
	controls->setLayout(layoutControls);
	
	// Assemble main window layout
	mainLayout = new QVBoxLayout;
	mainLayout->addWidget(settings);
	mainLayout->addStretch();
	mainLayout->addWidget(viewport);
	mainLayout->addStretch();
	mainLayout->addWidget(controls);
	mainLayout->addWidget(status);
	mainLayout->setContentsMargins(5, 5, 5, 5);
	mainLayout->setSpacing(5);
	mainLayout->setAlignment(viewport, Qt::AlignHCenter);
	this->setLayout(mainLayout);
	
	// Set up event filters
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
	unsigned int valueR, valueG, valueB;
	
	for(unsigned int i = 0; i < 256; ++i)
	{
		switch(type)
		{
			case RAINBOW:
				if (i <  50)
				{
					valueR = 0;
					valueG = static_cast<unsigned int>(255.0 * (static_cast<double>(i) - 0.0) / 50.0);
					valueB = 255;
				}
				else if(i < 125)
				{
					valueR = 0;
					valueG = 255;
					valueB = static_cast<unsigned int>(255.0 * (125.0 - static_cast<double>(i)) / 75.0);
				}
				else if(i < 180)
				{
					valueR = static_cast<unsigned int>(255 * (static_cast<double>(i) - 125) / 55);
					valueG = 255;
					valueB = 0;
				}
				else
				{
					valueR = 255;
					valueG = static_cast<unsigned int>(255 * (255 - static_cast<double>(i)) / 75);
					valueB = 0;
				}
				break;
			
			case RGB:
				valueR = i;
				valueG = 255 - static_cast<unsigned int>(abs(2 * static_cast<int>(i) - 255));
				valueB = 255 - i;
				break;
			
			case RANDOM:
				valueR = static_cast<unsigned int>(255.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
				valueG = static_cast<unsigned int>(255.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
				valueB = static_cast<unsigned int>(255.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
				break;
			
			default:
				// GREYSCALE
				valueR = i;
				valueG = i;
				valueB = i;
		}
		
		if(i < 255) lut.append(qRgb(valueR xor invert, valueG xor invert, valueB xor invert));
		else lut.append(qRgb(valueR xor invert, valueG xor invert, valueB xor invert));
	}
	
	image->setColorTable(lut);
	return;
}

	

// ----------------------------- //
// SLOTs to change look-up table //
// ----------------------------- //

void WidgetDataViewer::selectLutGreyscale()
{
	currentLut = GREYSCALE;
	setUpLut(currentLut);
	plotChannelMap(currentChannel);
	return;
}

void WidgetDataViewer::selectLutRainbow()
{
	currentLut = RAINBOW;
	setUpLut(currentLut);
	plotChannelMap(currentChannel);
	return;
}

void WidgetDataViewer::selectLutRgb()
{
	currentLut = RGB;
	setUpLut(currentLut);
	plotChannelMap(currentChannel);
	return;
}

void WidgetDataViewer::selectLutRandom()
{
	currentLut = RANDOM;
	setUpLut(currentLut);
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
		
		long x = floor(static_cast<double>(a) / scale - offsetX);
		long y = floor(static_cast<double>(VIEWPORT_HEIGHT - b - 1) / scale - offsetY);
		
		QString text("Undefined");
		
		if(x >= 0 and static_cast<size_t>(x) < fips->dimension(1) and y >= 0 and static_cast<size_t>(y) < fips->dimension(2))
		{
			size_t position[3] = {static_cast<size_t>(x), static_cast<size_t>(y), currentChannel};
			double value = fips->data(position);
			
			if(std::isnan(value))
			{
				text = QString("Position: %1, %2   Value: blank").arg(x).arg(y);
			}
			else
			{
				text = QString("Position: %1, %2   Value: %3 ").arg(x).arg(y).arg(value);
				QString bunit = QString::fromStdString(fips->unit());
				if(not bunit.isEmpty()) text.append(bunit);
			}
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
	if(not fips->dimension()) return;
	
	QMenu contextMenu("Context Menu", this);
	
	QMenu menuLut("Colour Scheme", this);
	menuLut.setIcon(iconFillColor);
	
	menuLut.addAction(actionLutGreyscale);
	menuLut.addAction(actionLutRainbow);
	menuLut.addAction(actionLutRgb);
	menuLut.addAction(actionLutRandom);
	menuLut.addSeparator();
	menuLut.addAction(actionRevert);
	menuLut.addAction(actionInvert);
	
	contextMenu.addAction(actionCopy);
	contextMenu.addSeparator();
	contextMenu.addAction(actionZoomIn);
	contextMenu.addAction(actionZoomOut);
	contextMenu.addAction(actionZoomToFit);
	contextMenu.addSeparator();
	contextMenu.addMenu(&menuLut);
	//contextMenu.addSeparator();
	//contextMenu.addAction(actionPrev);
	//contextMenu.addAction(actionNext);
	//contextMenu.addSeparator();
	//contextMenu.addAction(actionFirst);
	//contextMenu.addAction(actionLast);
	
	contextMenu.exec(viewport->mapToGlobal(where));
	
	return;
}



// ------------------------- //
// SLOT to set minimum level //
// ------------------------- //

void WidgetDataViewer::setLevelMin()
{
	bool ok;
	double value = (fieldLevelMin->text()).toDouble(&ok);
	
	if(ok and value < plotMax)
	{
		plotMin = value;
		plotChannelMap(currentChannel);
	}
	else
	{
		fieldLevelMin->setText(QString::number(plotMin));
	}
	
	return;
}



// ------------------------- //
// SLOT to set maximum level //
// ------------------------- //

void WidgetDataViewer::setLevelMax()
{
	bool ok;
	double value = (fieldLevelMax->text()).toDouble(&ok);
	
	if(ok and value > plotMin)
	{
		plotMax = value;
		plotChannelMap(currentChannel);
	}
	else
	{
		fieldLevelMax->setText(QString::number(plotMax));
	}
	
	return;
}



// -------------------------------- //
// SLOT to change transfer function //
// -------------------------------- //

void WidgetDataViewer::setTransferFunction(int which)
{
	transferFunction = which;
	plotChannelMap(currentChannel);
	
	return;
}



// ------------------------------- //
// SLOT to toggle reverted colours //
// ------------------------------- //

void WidgetDataViewer::toggleRev()
{
	revert = actionRevert->isChecked() ? 0xFF : 0x00;
	plotChannelMap(currentChannel);
	
	return;
}



// ------------------------------- //
// SLOT to toggle inverted colours //
// ------------------------------- //

void WidgetDataViewer::toggleInv()
{
	invert = actionInvert->isChecked() ? 0xFF : 0x00;
	setUpLut(currentLut);
	plotChannelMap(currentChannel);
	
	return;
}



// --------------- //
// SLOT to zoom in //
// --------------- //

void WidgetDataViewer::zoomIn()
{
	if(scale < SCALE_MAX and fips->dimension())
	{
		double posX = static_cast<double>(VIEWPORT_WIDTH) / (2.0 * scale) - offsetX;
		double posY = static_cast<double>(VIEWPORT_HEIGHT) / (2.0 * scale) - offsetY;
		
		scale *= SCALE_FACTOR;
		
		offsetX = static_cast<double>(VIEWPORT_WIDTH) / (2.0 * scale) - posX;
		offsetY = static_cast<double>(VIEWPORT_HEIGHT) / (2.0 * scale) - posY;
		
		plotChannelMap(currentChannel);
	}
	
	return;
}



// ---------------- //
// SLOT to zoom out //
// ---------------- //

void WidgetDataViewer::zoomOut()
{
	if(scale > 1.0 / SCALE_MAX and fips->dimension())
	{
		double posX = static_cast<double>(VIEWPORT_WIDTH) / (2.0 * scale) - offsetX;
		double posY = static_cast<double>(VIEWPORT_HEIGHT) / (2.0 * scale) - offsetY;
		
		scale /= SCALE_FACTOR;
		
		offsetX = static_cast<double>(VIEWPORT_WIDTH) / (2.0 * scale) - posX;
		offsetY = static_cast<double>(VIEWPORT_HEIGHT) / (2.0 * scale) - posY;
		
		plotChannelMap(currentChannel);
	}
	
	return;
}



// ------------------- //
// SLOT to zoom to fit //
// ------------------- //

void WidgetDataViewer::zoomToFit()
{
	if(not fips->dimension()) return;
	
	scale = std::min(static_cast<double>(VIEWPORT_WIDTH) / static_cast<double>(fips->dimension(1)), static_cast<double>(VIEWPORT_HEIGHT) / static_cast<double>(fips->dimension(2)));
	offsetX = 0.0;
	offsetY = 0.0;
	
	plotChannelMap(currentChannel);
	
	return;
}



// ------------------------------ //
// SLOT to reset display settings //
// ------------------------------ //

void WidgetDataViewer::resetDisplaySettings()
{
	if(not fips->dimension()) return;
	
	plotMin = dataMin;
	plotMax = dataMax;
	if(actionRevert->isChecked()) actionRevert->trigger();
	if(actionInvert->isChecked()) actionInvert->trigger();
	actionLutRainbow->trigger();
	transferFunction = LINEAR;
	zoomToFit();
	
	fieldLevelMin->setText(QString::number(plotMin));
	fieldLevelMax->setText(QString::number(plotMax));
	fieldTransFunc->setCurrentIndex(transferFunction);
	
	return;
}



// ----------------------------------- //
// FUNCTION to copy image to clipboard //
// ----------------------------------- //

void WidgetDataViewer::copy()
{
    if(image and not image->isNull())
	{
		QClipboard *clipboard = QApplication::clipboard();
		clipboard->setImage(*image);
	}
	
	return;
}



// ----------------- //
// Mouse wheel event //
// ----------------- //

void WidgetDataViewer::wheelEvent(QWheelEvent *event)
{
	if(not fips->dimension()) return;
	
	int angle = event->delta();
	int vpX = event->x() - viewport->x();
	int vpY = VIEWPORT_HEIGHT - event->y() + viewport->y();
	
	if(vpX >= 0 and vpX < VIEWPORT_WIDTH and vpY >= 0 and vpY < VIEWPORT_HEIGHT)
	{
		double posX = static_cast<double>(vpX) / scale - offsetX;
		double posY = static_cast<double>(vpY) / scale - offsetY;
		
		if(angle < 0 and scale < SCALE_MAX) scale *= SCALE_FACTOR;
		else if(angle > 0 and scale > 1.0 / SCALE_MAX) scale /= SCALE_FACTOR;
		
		offsetX = static_cast<double>(vpX) / scale - posX;
		offsetY = static_cast<double>(vpY) / scale - posY;
		
		plotChannelMap(currentChannel);
		
		event->accept();
	}
	else
	{
		event->ignore();
	}
	
	return;
}



// ----------------- //
// Mouse press event //
// ----------------- //

void WidgetDataViewer::mousePressEvent(QMouseEvent *event)
{
	if(not fips->dimension()) return;
	
	if(event->button() == Qt::LeftButton)
	{
		int vpX = event->x() - viewport->x();
		int vpY = VIEWPORT_HEIGHT - event->y() + viewport->y();
		
		double posX = static_cast<double>(vpX) / scale - offsetX;
		double posY = static_cast<double>(vpY) / scale - offsetY;
		
		offsetX = static_cast<double>(VIEWPORT_WIDTH) / (2.0 * scale) - posX;
		offsetY = static_cast<double>(VIEWPORT_HEIGHT) / (2.0 * scale) - posY;
		
		plotChannelMap(currentChannel);
		
		event->accept();
	}
	else
	{
		event->ignore();
	}
	
	return;
}



// ----------- //
// Close event //
// ----------- //

void WidgetDataViewer::closeEvent(QCloseEvent *event)
{
	event->accept();
	return;
}
