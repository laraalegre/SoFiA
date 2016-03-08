/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 0.5.0 (SoFiA.h) - Source Finding Application                   ///
/// Copyright (C) 2013-2016 Tobias Westmeier                             ///
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

#define SOFIA_DEFAULT_SETTINGS ":/SoFiA_default_input.txt"
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

#include <iostream>

#include <QtGlobal>

#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QList>
#include <QtCore/QProcess>
#include <QtCore/QString>
#include <QtCore/QMimeData>
#include <QtCore/QVariant>

#include <QtGui/QCloseEvent>
#include <QtGui/QDropEvent>

// Import correct headers depending on Qt version:
#if QT_VERSION < 0x050000
	#include <QtGui/QMainWindow>
	#include <QtGui/QWidget>
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
	#include <QtGui/QFileDialog>
	#include <QtGui/QMessageBox>
	#include <QtGui/QWhatsThis>
#else
	#include <QtWidgets/QMainWindow>
	#include <QtWidgets/QWidget>
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
	void saveSettings();
	void saveSettingsAs();
	void saveLogAs();
	void clearLog();
	void parameterChanged();
	void resetToDefault();
	void showHandbook(const QString &page = QString("index.html"));
	void aboutSoFiA();
	void selectInputWeightsFile();
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
	void toggleToolbar();
	void togglePipeline();
	void togglePipeline(bool state);
	void toggleSaveOnExit();
	void toggleFullScreen();
	
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
	int  loadFile(QString &fileName);
	int  updateVariables();
	int  setFields();
	int  setDefaults();
	int  loadSession();
	int  saveSession();
	int  showMessage(int severity, QString &messageText, QString &statusText);
	void createInterface();
	void createWhatsThis();
	void updateActions();
	void updateWindowTitle();
	
	QIcon iconSoFiA;
	QIcon iconDocumentNew;
	QIcon iconDocumentOpen;
	QIcon iconDocumentPreview;
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
	
	QMenu      *menuFile;
	QMenu      *menuPipeline;
	QMenu      *menuView;
	QMenu      *menuSettings;
	QMenu      *menuHelp;
	
	QToolBar   *toolBar;
	
	QAction    *actionOpen;
	QAction    *actionSave;
	QAction    *actionSaveAs;
	QAction    *actionExit;
	QAction    *actionRun;
	QAction    *actionAbort;
	QAction    *actionDefault;
	QAction    *actionSaveLogAs;
	QAction    *actionClearLog;
	QAction    *actionShowCatalogue;
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
	QWidget    *tabOutFilter;
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
	QLineEdit    *tabInputFieldMask;
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
	QSpinBox     *tabInFilterFieldEdgeX;
	QSpinBox     *tabInFilterFieldEdgeY;
	QSpinBox     *tabInFilterFieldEdgeZ;
	QComboBox    *tabInFilterFieldStatistic;
	QHBoxLayout  *tabInFilterLayoutScaleXYZ;
	QWidget      *tabInFilterWidgetScaleXYZ;
	QCheckBox    *tabInFilterFieldScaleX;
	QCheckBox    *tabInFilterFieldScaleY;
	QCheckBox    *tabInFilterFieldScaleZ;
	
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
	QComboBox    *tabSourceFindingFieldKunit;
	QTextEdit    *tabSourceFindingFieldKernels;
	
	QGroupBox    *tabSourceFindingGroupBox2;
	QFormLayout  *tabSourceFindingForm2;
	QLineEdit    *tabSourceFindingFieldThreshold2;
	QComboBox    *tabSourceFindingFieldRmsMode2;
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
	QFormLayout  *tabMergingForm1;
	
	QSpinBox     *tabMergingFieldRadiusX;
	QSpinBox     *tabMergingFieldRadiusY;
	QSpinBox     *tabMergingFieldRadiusZ;
	QSpinBox     *tabMergingFieldMinSizeX;
	QSpinBox     *tabMergingFieldMinSizeY;
	QSpinBox     *tabMergingFieldMinSizeZ;
	
	QHBoxLayout  *tabMergingLayoutControls;
	QWidget      *tabMergingWidgetControls;
	QPushButton  *tabMergingButtonPrev;
	QPushButton  *tabMergingButtonNext;
	
	
	
	// Parametrisation tab
	
	QToolBox     *toolBoxPA;
	
	QVBoxLayout  *tabParametrisationLayout;
	
	QGroupBox    *tabParametrisationGroupBox1;
	QFormLayout  *tabParametrisationForm1;
	QCheckBox    *tabParametrisationButtonBusyFunction;
	QCheckBox    *tabParametrisationButtonMaskOpt;
	QCheckBox    *tabParametrisationButtonDilateMask;
	
	QGroupBox    *tabParametrisationGroupBox2;
	QFormLayout  *tabParametrisationForm2;
	QLineEdit    *tabParametrisationFieldRelMin;
	//QLineEdit    *tabParametrisationFieldRelMax;
	//QLineEdit    *tabParametrisationFieldRelKernel;
	//QCheckBox    *tabParametrisationButtonAutoKernel;
	QCheckBox    *tabParametrisationButtonRelPlot;
	
	QHBoxLayout  *tabParametrisationLayoutControls;
	QWidget      *tabParametrisationWidgetControls;
	QPushButton  *tabParametrisationButtonPrev;
	QPushButton  *tabParametrisationButtonNext;
	
	
	
	// Output filter tab
	
	QToolBox     *toolBoxOF;
	
	QVBoxLayout  *tabOutFilterLayout;
	
	QGroupBox    *tabOutFilterGroupBox1;
	QFormLayout  *tabOutFilterForm1;
	QHBoxLayout  *tabOutFilterLayoutW50;
	QWidget      *tabOutFilterWidgetW50;
	QHBoxLayout  *tabOutFilterLayoutW20;
	QWidget      *tabOutFilterWidgetW20;
	QHBoxLayout  *tabOutFilterLayoutFpeak;
	QWidget      *tabOutFilterWidgetFpeak;
	QHBoxLayout  *tabOutFilterLayoutFint;
	QWidget      *tabOutFilterWidgetFint;
	QLineEdit    *tabOutFilterFieldW50Min;
	QLineEdit    *tabOutFilterFieldW50Max;
	QLineEdit    *tabOutFilterFieldW20Min;
	QLineEdit    *tabOutFilterFieldW20Max;
	QLineEdit    *tabOutFilterFieldFpeakMin;
	QLineEdit    *tabOutFilterFieldFpeakMax;
	QLineEdit    *tabOutFilterFieldFintMin;
	QLineEdit    *tabOutFilterFieldFintMax;
	QCheckBox    *tabOutFilterButtonW50;
	QCheckBox    *tabOutFilterButtonW20;
	QCheckBox    *tabOutFilterButtonFpeak;
	QCheckBox    *tabOutFilterButtonFint;
	
	QHBoxLayout  *tabOutFilterLayoutControls;
	QWidget      *tabOutFilterWidgetControls;
	QPushButton  *tabOutFilterButtonPrev;
	QPushButton  *tabOutFilterButtonNext;
	
	
	
	// Output tab
	
	QToolBox     *toolBoxOP;
	
	QVBoxLayout  *tabOutputLayout;
	
	QGroupBox    *tabOutputGroupBox1;
	QFormLayout  *tabOutputForm1;
	
	QHBoxLayout  *tabOutputLayoutProducts;
	QWidget      *tabOutputWidgetProducts;
	QCheckBox    *tabOutputButtonFilteredCube;
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
	QCheckBox    *tabOutputButtonASCII;
	QCheckBox    *tabOutputButtonXML;
	QCheckBox    *tabOutputButtonSQL;
	QCheckBox    *tabOutputButtonCompress;
	QCheckBox    *tabOutputButtonOverwrite;
	
	QGroupBox    *tabOutputGroupBox2;
	QFormLayout  *tabOutputForm2;
	QGridLayout  *tabOutputLayoutParameters;
	QWidget      *tabOutputWidgetParameters;
	QLabel       *tabOutputLabelParameters;
	QLabel       *tabOutputLabelParametersBasic;
	QLabel       *tabOutputLabelParametersBounds;
	QLabel       *tabOutputLabelParametersStat;
	QLabel       *tabOutputLabelParametersWCS;
	QLabel       *tabOutputLabelParametersPhysical;
	QLabel       *tabOutputLabelParametersBFFree;
	QLabel       *tabOutputLabelParametersBFPhys;
	
	QCheckBox    *tabOutputButtonParameter_id;
	QCheckBox    *tabOutputButtonParameter_name;
	QCheckBox    *tabOutputButtonParameter_x_geo;
	QCheckBox    *tabOutputButtonParameter_y_geo;
	QCheckBox    *tabOutputButtonParameter_z_geo;
	QCheckBox    *tabOutputButtonParameter_x;
	QCheckBox    *tabOutputButtonParameter_y;
	QCheckBox    *tabOutputButtonParameter_z;
	QCheckBox    *tabOutputButtonParameter_x_min;
	QCheckBox    *tabOutputButtonParameter_x_max;
	QCheckBox    *tabOutputButtonParameter_y_min;
	QCheckBox    *tabOutputButtonParameter_y_max;
	QCheckBox    *tabOutputButtonParameter_z_min;
	QCheckBox    *tabOutputButtonParameter_z_max;
	QCheckBox    *tabOutputButtonParameter_n_pix;
	//QCheckBox    *tabOutputButtonParameter_snr_min;
	//QCheckBox    *tabOutputButtonParameter_snr_max;
	//QCheckBox    *tabOutputButtonParameter_snr_sum;
	QCheckBox    *tabOutputButtonParameter_n_pos;
	QCheckBox    *tabOutputButtonParameter_n_neg;
	QCheckBox    *tabOutputButtonParameter_rel;
	QCheckBox    *tabOutputButtonParameter_bf_a;
	QCheckBox    *tabOutputButtonParameter_bf_b1;
	QCheckBox    *tabOutputButtonParameter_bf_b2;
	QCheckBox    *tabOutputButtonParameter_bf_c;
	QCheckBox    *tabOutputButtonParameter_bf_chi2;
	QCheckBox    *tabOutputButtonParameter_bf_flag;
	QCheckBox    *tabOutputButtonParameter_bf_f_int;
	QCheckBox    *tabOutputButtonParameter_bf_f_peak;
	QCheckBox    *tabOutputButtonParameter_bf_w;
	QCheckBox    *tabOutputButtonParameter_bf_w20;
	QCheckBox    *tabOutputButtonParameter_bf_w50;
	QCheckBox    *tabOutputButtonParameter_bf_xe;
	QCheckBox    *tabOutputButtonParameter_bf_xp;
	QCheckBox    *tabOutputButtonParameter_bf_z;
	QCheckBox    *tabOutputButtonParameter_ell_maj;
	QCheckBox    *tabOutputButtonParameter_ell_min;
	QCheckBox    *tabOutputButtonParameter_ell_pa;
	QCheckBox    *tabOutputButtonParameter_f_peak;
	QCheckBox    *tabOutputButtonParameter_f_int;
	QCheckBox    *tabOutputButtonParameter_f_wm50;
	QCheckBox    *tabOutputButtonParameter_rms;
	QCheckBox    *tabOutputButtonParameter_w20;
	QCheckBox    *tabOutputButtonParameter_w50;
	QCheckBox    *tabOutputButtonParameter_wm50;
	QCheckBox    *tabOutputButtonParameter_ra;
	QCheckBox    *tabOutputButtonParameter_dec;
	QCheckBox    *tabOutputButtonParameter_lon;
	QCheckBox    *tabOutputButtonParameter_lat;
	QCheckBox    *tabOutputButtonParameter_freq;
	QCheckBox    *tabOutputButtonParameter_velo;
	
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
