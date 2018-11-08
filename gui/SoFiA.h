/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.2.1 (SoFiA.h) - Source Finding Application                   ///
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

#ifndef PIPELINEINTERFACE_H
#define PIPELINEINTERFACE_H

#define SOFIA_VERSION_NUMBER "1.2.1-beta"

#define SOFIA_DEFAULT_SETTINGS ":/SoFiA_default_input.txt"
#define SOFIA_PARSET_EXTRAGALACTIC ":/parsets/SoFiA_parset_extragalactic.par"
#define SOFIA_PARSET_EXTRAGALACTIC_ARTEFACTS ":/parsets/SoFiA_parset_extragalactic_artefacts.par"
#define SOFIA_PARSET_CONTINUUM ":/parsets/SoFiA_parset_continuum.par"
#define SOFIA_TEMP_FILE        "SoFiA.session"
#define SOFIA_RC_FILE          ".SoFiA"

#define SOFIA_RC_TOOLBAR       "showToolbar"
#define SOFIA_RC_PIPELINE      "showPipeline"
#define SOFIA_RC_PIPELINEPOS   "pipelinePosition"
#define SOFIA_RC_SESSION       "saveSession"
#define SOFIA_RC_WINWIDTH      "windowWidth"
#define SOFIA_RC_WINHEIGHT     "windowHeight"

#define MESSAGE_INFO    0
#define MESSAGE_WARNING 1
#define MESSAGE_ERROR   2

#define KERNEL_SCALE_FACTOR 100.0
#define RELMIN_SCALE_FACTOR 100.0

#include <iostream>
#include <ctime>

#include <QtGlobal>

#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QList>
#include <QtCore/QProcess>
#include <QtCore/QString>
#include <QtCore/QMimeData>
#include <QtCore/QVariant>
#include <QtCore/QRegExp>

#include <QtGui/QCloseEvent>
#include <QtGui/QDropEvent>

// Import correct headers depending on Qt version:
#if QT_VERSION < 0x050000
	#include <QtGui/QMainWindow>
	#include <QtGui/QWidget>
	#include <QtGui/QFrame>
	#include <QtGui/QDockWidget>
	#include <QtGui/QToolBox>
	#include <QtGui/QStyle>
	#include <QtGui/QMenuBar>
	#include <QtGui/QToolBar>
	#include <QtGui/QStatusBar>
	#include <QtGui/QAction>
	#include <QtGui/QLayout>
	#include <QtGui/QFormLayout>
	#include <QtGui/QGroupBox>
	#include <QtGui/QProgressBar>
	#include <QtGui/QLabel>
	#include <QtGui/QLineEdit>
	#include <QtGui/QTextEdit>
	#include <QtGui/QScrollBar>
	#include <QtGui/QPushButton>
	#include <QtGui/QCheckBox>
	#include <QtGui/QSpinBox>
	#include <QtGui/QComboBox>
	#include <QtGui/QSlider>
	#include <QtGui/QFileDialog>
	#include <QtGui/QMessageBox>
	#include <QtGui/QWhatsThis>
#else
	#include <QtWidgets/QMainWindow>
	#include <QtWidgets/QWidget>
	#include <QtWidgets/QFrame>
	#include <QtWidgets/QDockWidget>
	#include <QtWidgets/QToolBox>
	#include <QtWidgets/QStyle>
	#include <QtWidgets/QMenuBar>
	#include <QtWidgets/QToolBar>
	#include <QtWidgets/QStatusBar>
	#include <QtWidgets/QAction>
	#include <QtWidgets/QLayout>
	#include <QtWidgets/QFormLayout>
	#include <QtWidgets/QGroupBox>
	#include <QtWidgets/QProgressBar>
	#include <QtWidgets/QLabel>
	#include <QtWidgets/QLineEdit>
	#include <QtWidgets/QTextEdit>
	#include <QtWidgets/QScrollBar>
	#include <QtWidgets/QPushButton>
	#include <QtWidgets/QCheckBox>
	#include <QtWidgets/QSpinBox>
	#include <QtWidgets/QComboBox>
	#include <QtWidgets/QSlider>
	#include <QtWidgets/QFileDialog>
	#include <QtWidgets/QMessageBox>
	#include <QtWidgets/QWhatsThis>
#endif

#include "WidgetSpreadsheet.h"

class SoFiA : public QMainWindow
{
	Q_OBJECT
	
public:
	SoFiA(int argc, char *argv[]);
	~SoFiA();
	
private slots:
	void selectInputDataFile();
	void selectOutputDirectory();
	void loadSettings();
	void loadParsetExtragalactic();
	void loadParsetExtragalacticArtefacts();
	void loadParsetContinuum();
	void saveSettings();
	void saveSettingsAs();
	void saveLogAs();
	void clearLog();
	void parameterChanged();
	void resetToDefault();
	void showHandbook(const QString &page = QString("index.html"));
	void aboutSoFiA();
	void selectInputWeightsFile();
	void selectInputFlagFile();
	void selectInputMaskFile();
	void selectOpticalCatalogFile();
	void displayPrevTab();
	void displayNextTab();
	void updateFields();
	void runPipeline();
	void pipelineProcessReadStd();
	void pipelineProcessReadErr();
	void pipelineProcessStarted();
	void pipelineProcessFinished(int exitCode, QProcess::ExitStatus exitStatus);
	void pipelineProcessCancel();
	void pipelineProcessError(QProcess::ProcessError error);
	void showCatalogue();
	void showCube();
	void showFilteredCube();
	void showNoiseCube();
	void showMask();
	void showMom0();
	void showMom1();
	void toggleToolbar();
	void togglePipeline();
	void togglePipeline(bool state);
	void toggleSaveOnExit();
	void toggleFullScreen();
	void printPositivityWarning(bool checked);
	
private:
	QString currentFileName;
	QString settingsFileName;
	bool settingsChanged;
	QMap<QString, QString> parameters;
	bool shutdownInitiated;
	
	bool settingsPipeline;
	int  settingsPipelinePosition;
	bool settingsToolbar;
	bool settingsSession;
	int  settingsWindowWidth;
	int  settingsWindowHeight;
	
	QProcess   *pipelineProcess;
	
	WidgetSpreadsheet *spreadsheet;
	
	int  selectFile(QLineEdit *target, bool isDirectory = false);
	int  loadFile(const QString &fileName, QString &version);
	int  updateVariables();
	int  setFields();
	int  setDefaults(const QString &fileName = QString::fromUtf8(SOFIA_DEFAULT_SETTINGS));
	int  loadSession();
	int  saveSession();
	int  showMessage(int severity, QString &messageText, QString &statusText);
	void createInterface();
	void createWhatsThis();
	void updateActions();
	void updateWindowTitle();
	QString extractFileName(QString &extension);
	
	QIcon iconSoFiA;
	QIcon iconDocumentNew;
	QIcon iconDocumentOpen;
	QIcon iconTextCsv;
	QIcon iconDocumentSave;
	QIcon iconDocumentSaveAs;
	QIcon iconApplicationExit;
	QIcon iconDialogOkApply;
	QIcon iconDialogCancel;
	QIcon iconDialogClose;
	QIcon iconGoPreviousView;
	QIcon iconGoNextView;
	QIcon iconEditClearList;
	QIcon iconHelpContents;
	QIcon iconHelpAbout;
	QIcon iconTaskComplete;
	QIcon iconTaskReject;
	QIcon iconWhatsThis;
	QIcon iconFullScreen;
	QIcon iconFolderImage;
	QIcon iconDrawRectangle;
	QIcon iconDrawCube;
	QIcon iconImageXGeneric;
	
	QMenu      *menuFile;
	QMenu      *menuLoadParset;
	QMenu      *menuPipeline;
	QMenu      *menuView;
	QMenu      *menuShowImage;
	QMenu      *menuSettings;
	QMenu      *menuHelp;
	
	QToolBar   *toolBar;
	
	QAction    *actionOpen;
	QAction    *actionLoadParsetExtragalactic;
	QAction    *actionLoadParsetExtragalacticArtefacts;
	QAction    *actionLoadParsetContinuum;
	QAction    *actionSave;
	QAction    *actionSaveAs;
	QAction    *actionExit;
	QAction    *actionRun;
	QAction    *actionAbort;
	QAction    *actionDefault;
	QAction    *actionSaveLogAs;
	QAction    *actionClearLog;
	QAction    *actionShowCatalogue;
	QAction    *actionShowCube;
	QAction    *actionShowFilteredCube;
	QAction    *actionShowNoiseCube;
	QAction    *actionShowMask;
	QAction    *actionShowMom0;
	QAction    *actionShowMom1;
	QAction    *actionToolbar;
	QAction    *actionPipeline;
	QAction    *actionSaveOnExit;
	QAction    *actionFullScreen;
	QAction    *actionHelp;
	QAction    *actionWhatsThis;
	QAction    *actionAbout;
	
	QWidget     *widgetMain;
	QVBoxLayout *layoutMain;
	
	QTabWidget *tabs;
	QWidget    *tabInput;
	QWidget    *tabInFilter;
	QWidget    *tabSourceFinding;
	QWidget    *tabMerging;
	QWidget    *tabParametrisation;
	QWidget    *tabOutput;
	
	
	
	// Input tab
	
	QToolBox     *toolBoxIP;
	
	QVBoxLayout  *tabInputLayout;
	QHBoxLayout  *tabInputLayoutData;
	QWidget      *tabInputWidgetData;
	QHBoxLayout  *tabInputLayoutWeights;
	QWidget      *tabInputWidgetWeights;
	QHBoxLayout  *tabInputLayoutMask;
	QWidget      *tabInputWidgetMask;
	
	
	QGroupBox    *tabInputGroupBox1;
	QFormLayout  *tabInputForm1;
	QLineEdit    *tabInputFieldData;
	QPushButton  *tabInputButtonData;
	QCheckBox    *tabInputFieldInvertData;
	QLineEdit    *tabInputFieldMask;
	QLineEdit    *tabInputFieldSources;
	QPushButton  *tabInputButtonMask;
	QLineEdit    *tabInputFieldWeights;
	QPushButton  *tabInputButtonWeights;
	QLineEdit    *tabInputFieldWeightsFunction;
	
	QGroupBox    *tabInputGroupBox2;
	QFormLayout  *tabInputForm2;
	QLineEdit    *tabInputFieldCatalog;
	QPushButton  *tabInputButtonCatalog;
	QHBoxLayout  *tabInputLayoutCatalog;
	QWidget      *tabInputWidgetCatalog;
	QLineEdit    *tabInputFieldSpatialSize;
	QLineEdit    *tabInputFieldSpectralSize;
	QCheckBox    *tabInputFieldMultiCat;
	
	QGroupBox    *tabInputGroupBox3;
	QFormLayout  *tabInputForm3;
	QLineEdit    *tabInputFieldFlags;
	QHBoxLayout  *tabInputLayoutFlagCube;
	QWidget      *tabInputWidgetFlagCube;
	QLineEdit    *tabInputFieldFlagCube;
	QPushButton  *tabInputButtonFlagCube;
	
	
	QGroupBox    *tabInputGroupBox4;
	QFormLayout  *tabInputForm4;
	QLineEdit    *tabInputFieldSubcube;
	QComboBox    *tabInputFieldSubcubeMode;
	
	QHBoxLayout  *tabInputLayoutControls;
	QWidget      *tabInputWidgetControls;
	QPushButton  *tabInputButtonNext;
	
	
	
	// Input filter tab
	
	QToolBox     *toolBoxIF;
	
	QVBoxLayout  *tabInFilterLayout;
	
	QGroupBox    *tabInFilterGroupBox1;
	QFormLayout  *tabInFilterForm1;
	QLineEdit    *tabInFilterFieldSmoothingSpatialLon;
	QLineEdit    *tabInFilterFieldSmoothingSpatialLat;
	QLineEdit    *tabInFilterFieldSmoothingSpectral;
	QComboBox    *tabInFilterFieldKernel;
	QComboBox    *tabInFilterFieldBorder;
	
	QGroupBox    *tabInFilterGroupBox2;
	QFormLayout  *tabInFilterForm2;
	QWidget      *tabInFilterWidgetEdge;
	QHBoxLayout  *tabInFilterLayoutEdge;
	QLabel       *tabInFilterLabelEdgeX;
	QLabel       *tabInFilterLabelEdgeY;
	QLabel       *tabInFilterLabelEdgeZ;
	QSpinBox     *tabInFilterFieldEdgeX;
	QSpinBox     *tabInFilterFieldEdgeY;
	QSpinBox     *tabInFilterFieldEdgeZ;
	QComboBox    *tabInFilterFieldStatistic;
	QComboBox    *tabInFilterFieldFluxRange;
	QHBoxLayout  *tabInFilterLayoutScaleXYZ;
	QWidget      *tabInFilterWidgetScaleXYZ;
	QCheckBox    *tabInFilterFieldScaleX;
	QCheckBox    *tabInFilterFieldScaleY;
	QCheckBox    *tabInFilterFieldScaleZ;
	QComboBox    *tabInFilterFieldMethod;
	//QWidget      *tabInFilterWidgetGrid;
	//QHBoxLayout  *tabInFilterLayoutGrid;
	//QLabel       *tabInFilterLabelGridSpatial;
	//QLabel       *tabInFilterLabelGridSpectral;
	//QSpinBox     *tabInFilterFieldGridSpatial;
	//QSpinBox     *tabInFilterFieldGridSpectral;
	QWidget      *tabInFilterWidgetWindow;
	QHBoxLayout  *tabInFilterLayoutWindow;
	QLabel       *tabInFilterLabelWindowSpatial;
	QLabel       *tabInFilterLabelWindowSpectral;
	QSpinBox     *tabInFilterFieldWindowSpatial;
	QSpinBox     *tabInFilterFieldWindowSpectral;
	QComboBox    *tabInFilterFieldInterpolation;
	QFrame       *tabInFilterSeparator1;
	QFrame       *tabInFilterSeparator2;
	
	QGroupBox    *tabInFilterGroupBox3;
	QFormLayout  *tabInFilterForm3;
	QLineEdit    *tabInFilterField2d1dThreshold;
	QSpinBox     *tabInFilterField2d1dScaleXY;
	QSpinBox     *tabInFilterField2d1dScaleZ;
	QSpinBox     *tabInFilterField2d1dIterations;
	QCheckBox    *tabInFilterField2d1dPositivity;
	
	QHBoxLayout  *tabInFilterLayoutControls;
	QWidget      *tabInFilterWidgetControls;
	QPushButton  *tabInFilterButtonPrev;
	QPushButton  *tabInFilterButtonNext;
	
	
	
	// Source finding tab
	
	QToolBox     *toolBoxSF;
	
	QVBoxLayout  *tabSourceFindingLayout;
	
	QGroupBox    *tabSourceFindingGroupBox1;
	QWidget      *tabSourceFindingWidget1Left;
	QWidget      *tabSourceFindingWidget1Right;
	QFormLayout  *tabSourceFindingForm1Left;
	QFormLayout  *tabSourceFindingForm1Right;
	QHBoxLayout  *tabSourceFindingForm1Layout;
	QLineEdit    *tabSourceFindingFieldThreshold;
	QComboBox    *tabSourceFindingFieldEdgeMode;
	QComboBox    *tabSourceFindingFieldRmsMode;
	QComboBox    *tabSourceFindingFieldFluxRange;
	QComboBox    *tabSourceFindingFieldKunit;
	QTextEdit    *tabSourceFindingFieldKernels;
	
	QGroupBox    *tabSourceFindingGroupBox2;
	QFormLayout  *tabSourceFindingForm2;
	QLineEdit    *tabSourceFindingFieldThreshold2;
	QComboBox    *tabSourceFindingFieldRmsMode2;
	QComboBox    *tabSourceFindingFieldFluxRange2;
	QComboBox    *tabSourceFindingFieldClipMethod;
	
	QGroupBox    *tabSourceFindingGroupBox3;
	QFormLayout  *tabSourceFindingForm3;
	QLineEdit    *tabSourceFindingFieldPReq;
	QLineEdit    *tabSourceFindingFieldQReq;
	QSpinBox     *tabSourceFindingFieldMinScale;
	QSpinBox     *tabSourceFindingFieldMaxScale;
	QCheckBox    *tabSourceFindingMedianTest;
	QComboBox    *tabSourceFindingFieldVerbose;
	
	QHBoxLayout  *tabSourceFindingLayoutControls;
	QWidget      *tabSourceFindingWidgetControls;
	QPushButton  *tabSourceFindingButtonPrev;
	QPushButton  *tabSourceFindingButtonNext;
	
	
	
	// Merging tab
	
	QToolBox     *toolBoxME;
	
	QVBoxLayout  *tabMergingLayout;
	
	QGroupBox    *tabMergingGroupBox1;
	QWidget      *tabMergingWidget1Left;
	QWidget      *tabMergingWidget1Right;
	QFormLayout  *tabMergingForm1Left;
	QFormLayout  *tabMergingForm1Right;
	QHBoxLayout  *tabMergingForm1Layout;
	
	QSpinBox     *tabMergingFieldRadiusX;
	QSpinBox     *tabMergingFieldRadiusY;
	QSpinBox     *tabMergingFieldRadiusZ;
	QSpinBox     *tabMergingFieldMinSizeX;
	QSpinBox     *tabMergingFieldMinSizeY;
	QSpinBox     *tabMergingFieldMinSizeZ;
	
	QGroupBox    *tabMergingGroupBox2;
	QFormLayout  *tabMergingForm2;
	QCheckBox    *tabMergingButtonPositivity;
	QLabel       *tabMergingLabelWarning;
	
	QHBoxLayout  *tabMergingLayoutControls;
	QWidget      *tabMergingWidgetControls;
	QPushButton  *tabMergingButtonPrev;
	QPushButton  *tabMergingButtonNext;
	
	
	
	// Parametrisation tab
	
	QToolBox     *toolBoxPA;
	
	QVBoxLayout  *tabParametrisationLayout;
	QWidget      *tabParametrisationWidgetRelMin;
	QHBoxLayout  *tabParametrisationLayoutRelMin;
	QWidget      *tabParametrisationWidgetScaleKernel;
	QHBoxLayout  *tabParametrisationLayoutScaleKernel;
	
	QGroupBox    *tabParametrisationGroupBox1;
	QFormLayout  *tabParametrisationForm1;
	QCheckBox    *tabParametrisationButtonBusyFunction;
	QCheckBox    *tabParametrisationButtonMaskOpt;
	QCheckBox    *tabParametrisationButtonDilateMask;
	
	QGroupBox    *tabParametrisationGroupBox2;
	QFormLayout  *tabParametrisationForm2;
	QSlider      *tabParametrisationSliderRelMin;
	QLineEdit    *tabParametrisationFieldRelMin;
	QSlider      *tabParametrisationSliderScaleKernel;
	QLineEdit    *tabParametrisationFieldScaleKernel;
	QCheckBox    *tabParametrisationButtonRelPlot;
	
	QHBoxLayout  *tabParametrisationLayoutControls;
	QWidget      *tabParametrisationWidgetControls;
	QPushButton  *tabParametrisationButtonPrev;
	QPushButton  *tabParametrisationButtonNext;
	
	
	
	// Output tab
	
	QToolBox     *toolBoxOP;
	
	QVBoxLayout  *tabOutputLayout;
	
	QGroupBox    *tabOutputGroupBox1;
	QFormLayout  *tabOutputForm1;
	
	QHBoxLayout  *tabOutputLayoutProducts;
	QWidget      *tabOutputWidgetProducts;
	QHBoxLayout  *tabOutputLayoutProducts2;
	QWidget      *tabOutputWidgetProducts2;
	QCheckBox    *tabOutputButtonFilteredCube;
	QCheckBox    *tabOutputButtonNoiseCube;
	QCheckBox    *tabOutputButtonMask;
	QCheckBox    *tabOutputButtonMom0;
	QCheckBox    *tabOutputButtonMom1;
	QCheckBox    *tabOutputButtonCubelets;
	
	QHBoxLayout  *tabOutputLayoutFormat;
	QWidget      *tabOutputWidgetFormat;
	QLineEdit    *tabOutputFieldBaseName;
	QWidget      *tabOutputWidgetDirectory;
	QHBoxLayout  *tabOutputLayoutDirectory;
	QLineEdit    *tabOutputFieldDirectory;
	QPushButton  *tabOutputButtonDirectory;
	QLineEdit    *tabOutputFieldParameters;
	QLabel       *tabOutputLabelParameters;
	QCheckBox    *tabOutputButtonASCII;
	QCheckBox    *tabOutputButtonXML;
	QCheckBox    *tabOutputButtonSQL;
	QCheckBox    *tabOutputButtonCompress;
	QCheckBox    *tabOutputButtonOverwrite;
	
	QGroupBox    *tabOutputGroupBox2;
	QFormLayout  *tabOutputForm2;
	
	QHBoxLayout  *tabOutputLayoutControls;
	QWidget      *tabOutputWidgetControls;
	QPushButton  *tabOutputButtonPrev;
	QPushButton  *tabOutputButtonGo;
	
	
	
	// Pipeline output terminal
	
	QDockWidget  *dockWidgetOutput;
	QWidget      *widgetOutput;
	QTextEdit    *outputText;
	QProgressBar *outputProgress;
	QVBoxLayout  *outputLayout;
	
protected:
	void closeEvent(QCloseEvent *event);
	void dropEvent(QDropEvent *event);
	void dragMoveEvent(QDragMoveEvent *event);
	void dragEnterEvent(QDragEnterEvent *event);
};

#endif
