/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.1.0 (SoFiA.cpp) - Source Finding Application                 ///
/// Copyright (C) 2013-2017 Tobias Westmeier                             ///
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

#include "HelpBrowser.h"
#include "WidgetDataViewer.h"
#include "SoFiA.h"



// -----------
// Constructor
// -----------

SoFiA::SoFiA(int argc, char *argv[])
{
	// The following is required to distinguish between user- and shutdown-related changes:
	shutdownInitiated = false;
	
	// Set default settings:
	settingsPipeline = true;
	settingsPipelinePosition = Qt::TopDockWidgetArea;
	settingsToolbar  = true;
	settingsSession  = true;
	settingsWindowWidth  = 640;
	settingsWindowHeight = 480;
	
	// Read SoFiA settings file in home directory, if present:
	QString homeDirectory = QDir::homePath();
	
	if(not homeDirectory.isEmpty() and homeDirectory != "/")
	{
		settingsFileName = homeDirectory + "/" + SOFIA_RC_FILE;
		
		QFile inFile(settingsFileName);
		
		if(inFile.open(QIODevice::ReadOnly | QIODevice::Text))
		{
			std::cerr << "\nInformation: Reading user settings from file " << settingsFileName.toStdString() << "\n";
			
			QTextStream inStream(&inFile);
			QString keyname;
			QString value;
			
			while(not inStream.atEnd())
			{
				QString line = inStream.readLine().trimmed();
				
				if(not line.isEmpty() and not line.startsWith("#"))
				{
					keyname = line.section(QChar('='), 0, 0).trimmed();
					value   = line.section(QChar('='), 1).trimmed();
					
					if     (keyname == SOFIA_RC_PIPELINE)    settingsPipeline = (value == "true");
					else if(keyname == SOFIA_RC_TOOLBAR)     settingsToolbar  = (value == "true");
					else if(keyname == SOFIA_RC_SESSION)     settingsSession  = (value == "true");
					else if(keyname == SOFIA_RC_PIPELINEPOS)
					{
						settingsPipelinePosition = value.toInt();
						// Default to top if ill-defined:
						settingsPipelinePosition = settingsPipelinePosition == Qt::LeftDockWidgetArea or settingsPipelinePosition == Qt::RightDockWidgetArea or settingsPipelinePosition == Qt::TopDockWidgetArea or settingsPipelinePosition == Qt::BottomDockWidgetArea ? settingsPipelinePosition : Qt::TopDockWidgetArea;
					}
					else if(keyname == SOFIA_RC_WINWIDTH)
					{
						settingsWindowWidth = value.toInt();
						if(settingsWindowWidth < 480 or settingsWindowWidth > 2048) settingsWindowWidth = 640;
					}
					else if(keyname == SOFIA_RC_WINHEIGHT)
					{
						settingsWindowHeight = value.toInt();
						if(settingsWindowHeight < 360 or settingsWindowHeight > 2048) settingsWindowHeight = 480;
					}
				}
			}
			
			inFile.close();
		}
		else std::cerr << "\nWarning: No user settings found in home directory. Using default settings.\n";
	}
	else std::cerr << "\nWarning: Failed to locate home directory. Using default settings.\n";
	
	// Create pipeline process, but don't run it yet:
	pipelineProcess = new QProcess(this);
	connect(pipelineProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(pipelineProcessReadStd()));
	connect(pipelineProcess, SIGNAL(readyReadStandardError()), this, SLOT(pipelineProcessReadErr()));
	connect(pipelineProcess, SIGNAL(started()), this, SLOT(pipelineProcessStarted()));
	connect(pipelineProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(pipelineProcessFinished(int, QProcess::ExitStatus)));
	connect(pipelineProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(pipelineProcessError(QProcess::ProcessError)));
	
	// Create user interface:
	this->createInterface();
	this->setDefaults();
	this->setAcceptDrops(true);        // Enable drag and drop
	
	// Create catalogue viewer window, but don't show it yet:
	spreadsheet = new WidgetSpreadsheet(this);
	spreadsheet->hide();
	
	if(argc > 1)
	{
		// Load parameter file if specified:
		QString fileName = QString(argv[1]);
		
		if(not fileName.isEmpty() and fileName[0] != '/')
		{
			fileName.prepend('/');
			fileName.prepend(QDir::currentPath());
		}
		
		if(loadFile(fileName))
		{
			QString messageText = tr("<p>Failed to read input file %1.</p>").arg(fileName.section('/', -1));
			QString statusText  = tr("Failed to read input file %1.").arg(fileName.section('/', -1));
			showMessage(MESSAGE_ERROR, messageText, statusText);
		}
	}
	else
	{
		// Else try to restore previous session:
		loadSession();
	}
	
	settingsChanged = false;           // This will record whether any settings are changed by the user
	
	return;
}



// ----------
// Destructor
// ----------

SoFiA::~SoFiA()
{
	// Either save or remove session file:
	if(settingsSession)
	{
		saveSession();
	}
	else if(QFile::exists(SOFIA_TEMP_FILE) and not QFile::remove(SOFIA_TEMP_FILE))
	{
		std::cerr << "Warning: Failed to remove temporary session file on exit.\n";
	}
	
	// Save user settings:
	if(not settingsFileName.isEmpty())
	{
		QFile outFile(settingsFileName);
		
		if(not outFile.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			std::cerr << "\nError: Failed to save user settings to home directory.\n";
		}
		else
		{
			QTextStream outStream(&outFile);
			
			outStream << "# User-defined global settings for SoFiA\n\n";
			outStream << SOFIA_RC_TOOLBAR <<  "       =  " << (settingsToolbar  ? "true" : "false") << "\n";
			outStream << SOFIA_RC_PIPELINE <<  "      =  " << (settingsPipeline ? "true" : "false") << "\n";
			outStream << SOFIA_RC_PIPELINEPOS <<   "  =  " << this->dockWidgetArea(dockWidgetOutput) << "\n";
			outStream << SOFIA_RC_SESSION <<  "       =  " << (settingsSession  ? "true" : "false") << "\n";
			outStream << SOFIA_RC_WINWIDTH << "       =  " << (this->size()).width() << "\n";
			outStream << SOFIA_RC_WINHEIGHT << "      =  " << (this->size()).height() << endl;
			
			outFile.close();
		}
	}
	
	return;
}



// ---------------------------------
// Function to load previous session
// ---------------------------------

int SoFiA::loadSession()
{
	// Check if there is a temporary file from a previous session:
	if(QFile::exists(SOFIA_TEMP_FILE))
	{
		// Yes, there is one, so let’s load it to continue where we previously stopped:
		QString sessionFileName  = SOFIA_TEMP_FILE;
		QString originalFileName = currentFileName;        // Remember original current file name, if any
		
		if(loadFile(sessionFileName))
		{
			QString messageText = tr("<p>Failed to read previous session file.</p>");
			QString statusText = tr("Failed to read previous session file.");
			showMessage(MESSAGE_ERROR, messageText, statusText);
			return 1;
		}
		
		currentFileName = originalFileName;                // Revert current file name to original again
		updateWindowTitle();
		
		QString messageText = QString("Parameter settings from a previous SoFiA session were discovered in this directory and will be automatically restored.");
		QString statusText = tr("Previous session restored.");
		showMessage(MESSAGE_INFO, messageText, statusText);
	}
	
	return 0;
}



// --------------------------------
// Function to save current session
// --------------------------------

int SoFiA::saveSession()
{
	QString originalFileName = currentFileName;            // Remember original current file name, if any
	currentFileName = SOFIA_TEMP_FILE;                     // Set current file name to default temporary file
	saveSettings();                                        // Save settings under temporary file name
	
	currentFileName = originalFileName;                    // Revert current file name to original again
	updateWindowTitle();
	
	return 0;
}



// ----------------------------
// Function to display messages
// ----------------------------

int SoFiA::showMessage(int severity, QString &messageText, QString &statusText)
{
	if(severity < MESSAGE_INFO or severity > MESSAGE_ERROR) severity = MESSAGE_INFO;
	
	if(not statusText.isEmpty())
	{
		QString statusTitle[3] = {tr("Information"), tr("Warning"), tr("Error")};
		this->statusBar()->showMessage(QString("%1: %2").arg(statusTitle[severity]).arg(statusText));
	}
	
	if(not messageText.isEmpty())
	{
		QString titleText[3] = {tr("SoFiA - Information"), tr("SoFiA - Warning"), tr("SoFiA - Error")};
		QMessageBox messageBox;
		messageBox.setWindowTitle(titleText[severity]);
		messageBox.setText(messageText);
		if      (severity == MESSAGE_ERROR)   messageBox.setIcon(QMessageBox::Critical);
		else if (severity == MESSAGE_WARNING) messageBox.setIcon(QMessageBox::Warning);
		else                                  messageBox.setIcon(QMessageBox::Information);
		messageBox.exec();
	}
	
	return 0;
}



// ----------------------------------------------
// Function to select file or directory from disc
// ----------------------------------------------

int SoFiA::selectFile(QLineEdit *target, bool isDirectory)
{
	QString fileName;
	
	if(isDirectory == false) fileName = QFileDialog::getOpenFileName(this, tr("SoFiA - Select File"), QDir::currentPath());
	else fileName = QFileDialog::getExistingDirectory(this, tr("SoFiA - Select Directory"), QDir::currentPath());
	
	if(fileName.isEmpty() or target == 0) return 1;
	
	target->setText(fileName);
	
	updateFields();
	
	return 0;
}



// ------------------------
// Function to set defaults
// ------------------------

int SoFiA::setDefaults(const QString &fileName)
{
	if(loadFile(fileName))
	{
		QString messageText = tr("<p>Failed to load default parameters.</p><p>Please close the programme and check your installation. SoFiA will not function properly without the default parameters.</p>");
		QString statusText = tr("Failed to load default parameters.");
		showMessage(MESSAGE_ERROR, messageText, statusText);
		
		return 1;
	}
	
	currentFileName.clear();
	
	QString messageText = QString("");
	QString statusText = tr("Default parameters loaded.");
	showMessage(MESSAGE_INFO, messageText, statusText);
	
	settingsChanged = false;
	updateWindowTitle();
	
	return 0;
}



// ----------------------------
// Function to update variables
// ----------------------------

int SoFiA::updateVariables()
{
	// For each type of input field/button a separate set of commands is needed, because 
	// the access methods all differ (text(), value(), isChecked(), etc.)
	
	QList<QLineEdit*> widgetLineEdit = tabs->findChildren<QLineEdit*>();
	
	// WARNING: 'foreach()' is not a C++ statement, but a Qt macro!
	foreach(QLineEdit *w, widgetLineEdit)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			if(w->text() == "auto") parameters.insert(w->objectName(), "0.0"); // for scaleKernel
			else parameters.insert(w->objectName(), w->text());
		}
	}
	
	QList<QTextEdit*> widgetTextEdit = tabs->findChildren<QTextEdit*>();
	
	foreach(QTextEdit *w, widgetTextEdit)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			parameters.insert(w->objectName(), w->toPlainText());
		}
	}
	
	QList<QComboBox*> widgetComboBox = tabs->findChildren<QComboBox*>();
	
	foreach(QComboBox *w, widgetComboBox)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			//parameters.insert(w->objectName(), w->currentText());
			parameters.insert(w->objectName(), (w->itemData(w->currentIndex())).toString());
		}
	}
	
	QList<QSpinBox*> widgetSpinBox = tabs->findChildren<QSpinBox*>();
	
	foreach(QSpinBox *w, widgetSpinBox)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			QString value;
			value.setNum(w->value());
			
			parameters.insert(w->objectName(), value);
		}
	}
	
	QList<QAbstractButton*> widgetAbstractButton = tabs->findChildren<QAbstractButton*>();
	
	foreach(QAbstractButton *w, widgetAbstractButton)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			QString value = "false";
			if(w->isChecked()) value = "true";
			
			parameters.insert(w->objectName(), value);
		}
	}
	
	QList<QGroupBox*> widgetGroupBox = tabs->findChildren<QGroupBox*>();
	
	foreach(QGroupBox *w, widgetGroupBox)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			QString value = "false";
			if(w->isChecked()) value = "true";
			
			parameters.insert(w->objectName(), value);
		}
	}
	
	return 0;
}



//-----------------------
// Function to set fields
// ----------------------

int SoFiA::setFields()
{
	// For each type of input field/button a separate set of commands is needed, because 
	// the access methods all differ (text(), value(), isChecked(), etc.)
	
	QList<QLineEdit*> widgetLineEdit = tabs->findChildren<QLineEdit*>();
	
	foreach(QLineEdit *w, widgetLineEdit)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			w->setText(parameters.value(w->objectName()));
			if(w == tabParametrisationFieldScaleKernel) tabParametrisationSliderScaleKernel->setValue(static_cast<int>((parameters.value(w->objectName())).toDouble() * KERNEL_SCALE_FACTOR));
			else if(w == tabParametrisationFieldRelMin) tabParametrisationSliderRelMin->setValue(static_cast<int>((parameters.value(w->objectName())).toDouble() * RELMIN_SCALE_FACTOR));
		}
	}
	
	QList<QTextEdit*> widgetTextEdit = tabs->findChildren<QTextEdit*>();
	
	foreach(QTextEdit *w, widgetTextEdit)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			w->setPlainText(parameters.value(w->objectName()));
		}
	}
	
	QList<QComboBox*> widgetComboBox = tabs->findChildren<QComboBox*>();
	
	foreach(QComboBox *w, widgetComboBox)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			int index = w->findData(QVariant(parameters.value(w->objectName())));
			if(index >= 0) w->setCurrentIndex(index);
		}
	}
	
	QList<QSpinBox*> widgetSpinBox = tabs->findChildren<QSpinBox*>();
	
	foreach(QSpinBox *w, widgetSpinBox)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			QString value = parameters.value(w->objectName());
			
			w->setValue(value.toInt());
		}
	}
	
	QList<QAbstractButton*> widgetAbstractButton = tabs->findChildren<QAbstractButton*>();
	
	foreach(QAbstractButton *w, widgetAbstractButton)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			QString value = parameters.value(w->objectName());
			
			if(value == "True" or value == "true" or value == "TRUE" or value == "T" or value == "t" or value == "1" or value == "Yes" or value == "yes" or value == "YES" or value == "Y" or value == "y") w->setChecked(true);
			else w->setChecked(false);
		}
	}
	
	QList<QGroupBox*> widgetGroupBox = tabs->findChildren<QGroupBox*>();
	
	foreach(QGroupBox *w, widgetGroupBox)
	{
		if(parameters.contains(w->objectName()))      // Only existing parameters will get updated!
		{
			QString value = parameters.value(w->objectName());
			
			if(value == "True" or value == "true" or value == "TRUE" or value == "T" or value == "t" or value == "1" or value == "Yes" or value == "yes" or value == "YES" or value == "Y" or value == "y") w->setChecked(true);
			else w->setChecked(false);
		}
	}
	
	return 0;
}



// ---------------------
// Slot to load settings
// ---------------------

void SoFiA::loadSettings()
{
	QString newFileName = QFileDialog::getOpenFileName(this, tr("SoFiA - Load Parameters"), QDir::currentPath());
	
	if(not newFileName.isEmpty())
	{
		if(settingsChanged)
		{
			QMessageBox messageBox(this);
			messageBox.setWindowTitle(tr("SoFiA - Load Parameters"));
			messageBox.setText(tr("<p>This action will open a new file and override the current parameter settings. All unsaved changes will be lost.</p><p>Do you wish to open a new parameter file?</p>"));
			messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
			messageBox.setDefaultButton(QMessageBox::Ok);
			messageBox.setIcon(QMessageBox::Warning);
			int choice = messageBox.exec();
			
			if(choice != QMessageBox::Ok) return;
		}
		
		setDefaults();       // Load default settings first. This ensures that unset parameters will assume their default value.
		
		if(loadFile(newFileName))
		{
			QString messageText = tr("<p>Failed to read input file %1.</p>").arg(newFileName.section('/', -1));
			QString statusText = tr("Failed to read input file %1.").arg(newFileName.section('/', -1));
			showMessage(MESSAGE_ERROR, messageText, statusText);
		}
	}
	
	return;
}



// ---------------------
// Function to load file
// ---------------------

int SoFiA::loadFile(const QString &fileName)
{
	if(fileName.isEmpty()) return 0;
	
	QFile inFile(fileName);
	
	if(not inFile.open(QIODevice::ReadOnly | QIODevice::Text)) return 1;       // Error message should be generated by calling function.
	
	QTextStream inStream(&inFile);
	QString keyname;
	QString value;
	
	int counter = 0;
	
	while(not inStream.atEnd())
	{
		QString line = inStream.readLine().trimmed();
		
		if(not line.isEmpty() and not line.startsWith("#"))
		{
			keyname = line.section(QChar('='), 0, 0).trimmed();
			value   = line.section(QChar('='), 1).trimmed();
			
			if((fileName == SOFIA_DEFAULT_SETTINGS or parameters.contains(keyname)) and keyname.size() != 0)
			{
				parameters.insert(keyname, value);
				counter++;
			}
		}
	}
	
	inFile.close();
	
	if(counter == 0)
	{
		QString messageText = tr("<p>No valid parameters found in input file %1.</p>").arg(fileName.section('/', -1));
		QString statusText = tr("No valid parameters found in input file %1.").arg(fileName.section('/', -1));
		showMessage(MESSAGE_ERROR, messageText, statusText);
		return 1;
	}
	
	currentFileName = fileName;
	
	setFields();
	updateFields();
	
	QString messageText = QString("");
	QString statusText = tr("Parameters loaded from %1.").arg(currentFileName.section('/', -1));
	showMessage(MESSAGE_INFO, messageText, statusText);
	
	settingsChanged = false;
	updateWindowTitle();
	
	// Make appropriate activated sub-tabs visible after loading file:
	// Input:
	if     (tabInputGroupBox4->isChecked()) toolBoxIP->setCurrentIndex(toolBoxIP->indexOf(tabInputGroupBox4));
	else if(tabInputGroupBox3->isChecked()) toolBoxIP->setCurrentIndex(toolBoxIP->indexOf(tabInputGroupBox3));
	else if(tabInputGroupBox2->isChecked()) toolBoxIP->setCurrentIndex(toolBoxIP->indexOf(tabInputGroupBox2));
	else                                    toolBoxIP->setCurrentIndex(toolBoxIP->indexOf(tabInputGroupBox1));
	// Input filters:
	if     (tabInFilterGroupBox1->isChecked()) toolBoxIF->setCurrentIndex(toolBoxIF->indexOf(tabInFilterGroupBox1));
	else if(tabInFilterGroupBox2->isChecked()) toolBoxIF->setCurrentIndex(toolBoxIF->indexOf(tabInFilterGroupBox2));
	else if(tabInFilterGroupBox3->isChecked()) toolBoxIF->setCurrentIndex(toolBoxIF->indexOf(tabInFilterGroupBox3));
	else                                       toolBoxIF->setCurrentIndex(toolBoxIF->indexOf(tabInFilterGroupBox1));
	// Source finders:
	if     (tabSourceFindingGroupBox1->isChecked()) toolBoxSF->setCurrentIndex(toolBoxSF->indexOf(tabSourceFindingGroupBox1));
	else if(tabSourceFindingGroupBox2->isChecked()) toolBoxSF->setCurrentIndex(toolBoxSF->indexOf(tabSourceFindingGroupBox2));
	else if(tabSourceFindingGroupBox3->isChecked()) toolBoxSF->setCurrentIndex(toolBoxSF->indexOf(tabSourceFindingGroupBox3));
	else                                            toolBoxSF->setCurrentIndex(toolBoxSF->indexOf(tabSourceFindingGroupBox1));
	// Parameterisation:
	if     (tabParametrisationGroupBox1->isChecked()) toolBoxPA->setCurrentIndex(toolBoxPA->indexOf(tabParametrisationGroupBox1));
	else if(tabParametrisationGroupBox2->isChecked()) toolBoxPA->setCurrentIndex(toolBoxPA->indexOf(tabParametrisationGroupBox2));
	else                                              toolBoxPA->setCurrentIndex(toolBoxPA->indexOf(tabParametrisationGroupBox1));
	
	return 0;
}



// ---------------------
// Slot to save settings
// ---------------------

void SoFiA::saveSettings()
{
	if(currentFileName.isEmpty())
	{
		saveSettingsAs();
	}
	else
	{
		QFile outFile(currentFileName);
		
		if(not outFile.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			QString messageText = tr("<p>Failed to write to output file %1.</p>").arg(currentFileName.section('/', -1));
			QString statusText = tr("Failed to write to output file %1.").arg(currentFileName.section('/', -1));
			showMessage(MESSAGE_ERROR, messageText, statusText);
			
			currentFileName.clear();
			updateWindowTitle();
			
			return;
		}
		
		updateVariables();   // This will update all parameters before saving.
		
		QTextStream outStream(&outFile);
		
		// Write some header information first:
		time_t currentTime = time(0);
		outStream << "# SoFiA parameter file\n";
		outStream << "# Creator: SoFiA 1.2.0-beta\n";
		outStream << "# Date:    " << ctime(&currentTime) << '\n';
		
		// Then write the actual parameter settings:
		for(QMap<QString, QString>::iterator iter = parameters.begin(); iter != parameters.end(); iter++)
		{
			outStream << iter.key() << "\t=\t" << iter.value() << '\n';
		}
		
		outFile.close();
		
		QString messageText = QString("");
		QString statusText = tr("Parameters saved to %1.").arg(currentFileName.section('/', -1));
		showMessage(MESSAGE_INFO, messageText, statusText);
		
		if(currentFileName != SOFIA_TEMP_FILE)
		{
			settingsChanged = false;
			updateWindowTitle();
		}
	}
	
	return;
}



// ---------------------------
// Slot to save settings as...
// ---------------------------

void SoFiA::saveSettingsAs()
{
	QString fileName = QFileDialog::getSaveFileName(this, tr("SoFiA - Save Parameters"), QDir::currentPath());
	
	if(not fileName.isEmpty())
	{
		currentFileName = fileName;
		saveSettings();
	}
	
	return;
}



// ------------------------------------
// Slot to save pipeline messages as...
// ------------------------------------

void SoFiA::saveLogAs()
{
	QString fileName = QFileDialog::getSaveFileName(this, tr("SoFia - Save Pipeline Messages"), QDir::currentPath());
	
	if(fileName.isEmpty()) return;
	
	QFile file(fileName);
	
	if(file.open(QIODevice::WriteOnly))
	{
		QTextStream stream(&file);
		stream << outputText->toPlainText();
		
		file.close();
		
		QString messageText = QString("");
		QString statusText  = tr("Pipeline messages written to %1.").arg(fileName.section('/', -1));
		showMessage(MESSAGE_INFO, messageText, statusText);
	}
	else
	{
		QString messageText = tr("<p>Failed to write to file %1.</p><p>Pipeline messages not saved.</p>").arg(fileName.section('/', -1));
		QString statusText  = tr("Failed to write pipeline messages to %1.").arg(fileName.section('/', -1));
		showMessage(MESSAGE_ERROR, messageText, statusText);
	}
	
	return;
}



// -------------------------------
// Slot to clear pipeline messages
// -------------------------------

void SoFiA::clearLog()
{
	outputText->clear();
	outputProgress->setValue(0);
	
	QString messageText = tr("");
	QString statusText  = tr("Pipeline messages cleared.");
	showMessage(MESSAGE_INFO, messageText, statusText);
	
	updateActions();
	
	return;
}



// -------------------------
// Slot to update the fields
// -------------------------

void SoFiA::updateFields()
{
	// Activate or de-activate fields and buttons
	
	// Check ranges in optical source finding:
	double n = tabInputFieldSpatialSize->text().toDouble();
	if(n < 0.0) tabInputFieldSpatialSize->setText("0.0");
	n = tabInputFieldSpectralSize->text().toDouble();
	if(n < 0.0) tabInputFieldSpectralSize->setText("0.0");
	
	// Check ranges in CNHI finder:
	n = tabSourceFindingFieldPReq->text().toDouble();
	if(n < 0.0) tabSourceFindingFieldPReq->setText("0.0");
	n = tabSourceFindingFieldQReq->text().toDouble();
	if(n < 1.0) tabSourceFindingFieldQReq->setText("1.0");
	
	// Disable RMS mode and flux range fields in threshold finder if clip mode is 'absolute', and disable flux range for RMS mode 'negative':
	tabSourceFindingFieldRmsMode2->setEnabled(tabSourceFindingFieldClipMethod->currentIndex() == 0 and tabSourceFindingGroupBox2->isChecked());
	tabSourceFindingFieldFluxRange2->setEnabled(tabSourceFindingFieldClipMethod->currentIndex() == 0 and tabSourceFindingFieldRmsMode2->currentIndex() != 0 and tabSourceFindingGroupBox2->isChecked());
	
	// Disable flux range field in threshold finder if RMS mode is 'negative':
	tabSourceFindingFieldFluxRange->setEnabled(tabSourceFindingFieldRmsMode->currentIndex() != 0 and tabSourceFindingGroupBox1->isChecked());
	
	// Disable flux range field in noise scaling module if RMS statistic is 'negative':
	tabInFilterFieldFluxRange->setEnabled(tabInFilterFieldStatistic->currentIndex() != 0 and tabInFilterGroupBox2->isChecked());
	
	// Enable/disable writing of filtered cube when no filters are selected:
	tabOutputButtonFilteredCube->setEnabled(tabInFilterGroupBox1->isChecked() or tabInFilterGroupBox2->isChecked() or tabInFilterGroupBox3->isChecked() or not (tabInputFieldWeights->text()).isEmpty() or not (tabInputFieldWeightsFunction->text()).isEmpty());
	
	// Enable/disable writing of noise cube when noise scaling is not set:
	tabOutputButtonNoiseCube->setEnabled(tabInFilterGroupBox2->isChecked());
	
	// Check reliability threshold slider position:
	tabParametrisationFieldRelMin->setText(QString::number(static_cast<double>(tabParametrisationSliderRelMin->value()) / RELMIN_SCALE_FACTOR, 'f', 2));
	
	// Check kernel scale slider position:
	int tmpKernelScale = tabParametrisationSliderScaleKernel->value();
	if(tmpKernelScale == 0) tabParametrisationFieldScaleKernel->setText("auto");
	else tabParametrisationFieldScaleKernel->setText(QString::number(static_cast<double>(tmpKernelScale) / KERNEL_SCALE_FACTOR, 'f', 2));
	
	// Ensure that reliability range is within 0 and 1:
	n = tabParametrisationFieldRelMin->text().toDouble();
	if(n < 0.0) tabParametrisationFieldRelMin->setText("0.0");
	else if(n > 1.0) tabParametrisationFieldRelMin->setText("1.0");
	
	// Ensure that wavelet reconstruction threshold is >= 0:
	n = tabInFilterField2d1dThreshold->text().toDouble();
	if(n < 0.0) tabInFilterField2d1dThreshold->setText("0.0");
	
	// Ensure that S+C source finder threshold is >= 0:
	n = tabSourceFindingFieldThreshold->text().toDouble();
	if(n < 0.0) tabSourceFindingFieldThreshold->setText("0.0");
	
	// Ensure that threshold source finder threshold is >= 0:
	n = tabSourceFindingFieldThreshold2->text().toDouble();
	if(n < 0.0) tabSourceFindingFieldThreshold2->setText("0.0");
	
	// Ensure that smoothing lengths are >= 0.0:
	n = tabInFilterFieldSmoothingSpatialLon->text().toDouble();
	if(n < 0.0) tabInFilterFieldSmoothingSpatialLon->setText("0.0");
	n = tabInFilterFieldSmoothingSpatialLat->text().toDouble();
	if(n < 0.0) tabInFilterFieldSmoothingSpatialLat->setText("0.0");
	n = tabInFilterFieldSmoothingSpectral->text().toDouble();
	if(n < 0.0) tabInFilterFieldSmoothingSpectral->setText("0.0");
	
	// Noise scaling method
	//tabInFilterWidgetGrid->setEnabled(tabInFilterFieldMethod->currentIndex() == 1 and tabInFilterGroupBox2->isChecked());
	tabInFilterWidgetWindow->setEnabled(tabInFilterFieldMethod->currentIndex() == 1 and tabInFilterGroupBox2->isChecked());
	tabInFilterFieldInterpolation->setEnabled(tabInFilterFieldMethod->currentIndex() == 1 and tabInFilterGroupBox2->isChecked());
	tabInFilterWidgetEdge->setEnabled(tabInFilterFieldMethod->currentIndex() == 0 and tabInFilterGroupBox2->isChecked());
	tabInFilterWidgetScaleXYZ->setEnabled(tabInFilterFieldMethod->currentIndex() == 0 and tabInFilterGroupBox2->isChecked());
	
	//tabInFilterFieldGridSpatial->setValue(tabInFilterFieldGridSpatial->value() + tabInFilterFieldGridSpatial->value() % 2);
	//tabInFilterFieldGridSpectral->setValue(tabInFilterFieldGridSpectral->value() + tabInFilterFieldGridSpectral->value() % 2);
	tabInFilterFieldWindowSpatial->setValue(tabInFilterFieldWindowSpatial->value() + (1 - tabInFilterFieldWindowSpatial->value() % 2));
	tabInFilterFieldWindowSpectral->setValue(tabInFilterFieldWindowSpectral->value() + (1 - tabInFilterFieldWindowSpectral->value() % 2));
	
	// Disable source parameter list if not catalogue format selected:
	tabOutputFieldParameters->setEnabled(tabOutputButtonASCII->isChecked() or tabOutputButtonXML->isChecked() or tabOutputButtonSQL->isChecked());
	
	// Set icons on vertical tabs:
	if((tabInputFieldData->text()).isEmpty())    toolBoxIP->setItemIcon(0, iconTaskReject);
	else                                         toolBoxIP->setItemIcon(0, iconTaskComplete);
	if(tabInputGroupBox4->isChecked())           toolBoxIP->setItemIcon(1, iconTaskComplete);
	else                                         toolBoxIP->setItemIcon(1, iconTaskReject);
	if(tabInputGroupBox3->isChecked())           toolBoxIP->setItemIcon(2, iconTaskComplete);
	else                                         toolBoxIP->setItemIcon(2, iconTaskReject);
	if(tabInputGroupBox2->isChecked())           toolBoxIP->setItemIcon(3, iconTaskComplete);
	else                                         toolBoxIP->setItemIcon(3, iconTaskReject);
	
	if(tabInFilterGroupBox1->isChecked())        toolBoxIF->setItemIcon(0, iconTaskComplete);
	else                                         toolBoxIF->setItemIcon(0, iconTaskReject);
	if(tabInFilterGroupBox2->isChecked())        toolBoxIF->setItemIcon(1, iconTaskComplete);
	else                                         toolBoxIF->setItemIcon(1, iconTaskReject);
	if(tabInFilterGroupBox3->isChecked())        toolBoxIF->setItemIcon(2, iconTaskComplete);
	else                                         toolBoxIF->setItemIcon(2, iconTaskReject);
	
	if(tabSourceFindingGroupBox1->isChecked())   toolBoxSF->setItemIcon(0, iconTaskComplete);
	else                                         toolBoxSF->setItemIcon(0, iconTaskReject);
	if(tabSourceFindingGroupBox3->isChecked())   toolBoxSF->setItemIcon(1, iconTaskComplete);
	else                                         toolBoxSF->setItemIcon(1, iconTaskReject);
	if(tabSourceFindingGroupBox2->isChecked())   toolBoxSF->setItemIcon(2, iconTaskComplete);
	else                                         toolBoxSF->setItemIcon(2, iconTaskReject);
	
	if(tabMergingGroupBox1->isChecked())         toolBoxME->setItemIcon(0, iconTaskComplete);
	else                                         toolBoxME->setItemIcon(0, iconTaskReject);
	if(tabMergingButtonPositivity->isChecked())  toolBoxME->setItemIcon(1, iconTaskComplete);
	else                                         toolBoxME->setItemIcon(1, iconTaskReject);
	
	if(tabParametrisationGroupBox1->isChecked()) toolBoxPA->setItemIcon(0, iconTaskComplete);
	else                                         toolBoxPA->setItemIcon(0, iconTaskReject);
	if(tabParametrisationGroupBox2->isChecked()) toolBoxPA->setItemIcon(1, iconTaskComplete);
	else                                         toolBoxPA->setItemIcon(1, iconTaskReject);
	
	if(tabOutputButtonASCII->isChecked() or tabOutputButtonXML->isChecked() or tabOutputButtonSQL->isChecked() or (tabOutputButtonFilteredCube->isEnabled() and tabOutputButtonFilteredCube->isChecked()) or tabOutputButtonMask->isChecked() or tabOutputButtonMom0->isChecked() or tabOutputButtonMom1->isChecked() or tabOutputButtonCubelets->isChecked()) toolBoxOP->setItemIcon(0, iconTaskComplete);
	else toolBoxOP->setItemIcon(0, iconTaskReject);
	
	if(not (tabOutputFieldBaseName->text()).isEmpty() or not (tabOutputFieldDirectory->text()).isEmpty() or tabOutputButtonCompress->isChecked() or tabOutputButtonOverwrite->isChecked()) toolBoxOP->setItemIcon(1, iconTaskComplete);
	else toolBoxOP->setItemIcon(1, iconTaskReject);
	
	updateActions();
	
	return;
}



// --------------------------------------
// Function to update actions and buttons
// --------------------------------------

void SoFiA::updateActions()
{
	// Activate or de-activate actions and buttons
	
	tabOutputButtonGo->setEnabled(not (tabInputFieldData->text()).isEmpty() and pipelineProcess->state() == QProcess::NotRunning);
	actionRun        ->setEnabled(not (tabInputFieldData->text()).isEmpty() and pipelineProcess->state() == QProcess::NotRunning);
	
	actionAbort->setEnabled(pipelineProcess->state() == QProcess::Running);
	actionExit->setEnabled(pipelineProcess->state() == QProcess::NotRunning);
	
	actionSaveLogAs->setEnabled(outputText->toPlainText() != "" and pipelineProcess->state() == QProcess::NotRunning);
	actionClearLog->setEnabled(outputText->toPlainText() != "" and pipelineProcess->state() == QProcess::NotRunning);
	
	actionShowCatalogue->setEnabled(not (tabInputFieldData->text()).isEmpty());
	
	actionShowCube->setEnabled(not (tabInputFieldData->text()).isEmpty() and not tabOutputButtonCompress->isChecked());
	actionShowFilteredCube->setEnabled(not (tabInputFieldData->text()).isEmpty() and tabOutputButtonFilteredCube->isEnabled() and tabOutputButtonFilteredCube->isChecked() and not tabOutputButtonCompress->isChecked());
	actionShowNoiseCube->setEnabled(not (tabInputFieldData->text()).isEmpty() and tabOutputButtonNoiseCube->isEnabled() and tabOutputButtonNoiseCube->isChecked() and not tabOutputButtonCompress->isChecked());
	actionShowMask->setEnabled(not (tabInputFieldData->text()).isEmpty() and tabOutputButtonMask->isChecked() and not tabOutputButtonCompress->isChecked());
	actionShowMom0->setEnabled(not (tabInputFieldData->text()).isEmpty() and tabOutputButtonMom0->isChecked() and not tabOutputButtonCompress->isChecked());
	actionShowMom1->setEnabled(not (tabInputFieldData->text()).isEmpty() and tabOutputButtonMom1->isChecked() and not tabOutputButtonCompress->isChecked());
	menuShowImage->setEnabled(not (tabInputFieldData->text()).isEmpty() and not tabOutputButtonCompress->isChecked());
	
	return;
}



// -------------------------------------
// Slot to react to any parameter change
// -------------------------------------

void SoFiA::parameterChanged()
{
	settingsChanged = true;
	updateWindowTitle();
	
	return;
}



// ------------------------------
// Slot to select input data cube
// ------------------------------

void SoFiA::selectInputDataFile()
{
	selectFile(tabInputFieldData);
	
	return;
}



// ---------------------------------
// Slot to select input weights cube
// ---------------------------------

void SoFiA::selectInputWeightsFile()
{
	selectFile(tabInputFieldWeights);
	
	return;
}



// ------------------------------
// Slot to select input flag cube
// ------------------------------

void SoFiA::selectInputFlagFile()
{
	selectFile(tabInputFieldFlagCube);
	
	return;
}



// -------------------------------------
// Slot to select optical catalogue file
// -------------------------------------

void SoFiA::selectOpticalCatalogFile()
{
	selectFile(tabInputFieldCatalog);
	
	return;
}



// ------------------------------
// Slot to select input mask cube
// ------------------------------

void SoFiA::selectInputMaskFile()
{
	selectFile(tabInputFieldMask);
	
	return;
}



// -------------------------------
// Slot to select output directory
// -------------------------------

void SoFiA::selectOutputDirectory()
{
	selectFile(tabOutputFieldDirectory, true);
	
	return;
}



// ----------------------------
// Slot to display previous tab
// ----------------------------

void SoFiA::displayPrevTab()
{
	int i = tabs->currentIndex();
	
	if(i > 0) tabs->setCurrentIndex(i - 1);
	
	return;
}



// ------------------------
// Slot to display next tab
// ------------------------

void SoFiA::displayNextTab()
{
	int i = tabs->currentIndex();
	
	if(i < 6) tabs->setCurrentIndex(i + 1);
	
	return;
}



// ---------------------------------
// Slot to reset settings to default
// ---------------------------------

void SoFiA::resetToDefault()
{
	if(settingsChanged)
	{
		QMessageBox messageBox(this);
		messageBox.setWindowTitle(tr("SoFiA - New Parameter File"));
		messageBox.setText(tr("<p>This action will reset all parameters to their default values and close the current file. All unsaved changes will be lost.</p><p>Do you wish to create a new parameter file?</p>"));
		messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
		messageBox.setDefaultButton(QMessageBox::Ok);
		messageBox.setIcon(QMessageBox::Warning);
		
		int choice = messageBox.exec();
		if(choice != QMessageBox::Ok) return;
	}
	
	setDefaults();
	
	return;
}



// ------------------------------------------
// Slot to load parameter set 'Extragalactic'
// ------------------------------------------

void SoFiA::loadParsetExtragalactic()
{
	if(settingsChanged)
	{
		QMessageBox messageBox(this);
		messageBox.setWindowTitle(tr("SoFiA - New Parameter File"));
		messageBox.setText(tr("<p>This action will load the default parameter set %1 and close the current file. All unsaved changes will be lost.</p><p>Do you wish to create a new parameter file?</p>").arg(QString::fromUtf8("‘Extragalactic HI Cube’")));
		messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
		messageBox.setDefaultButton(QMessageBox::Ok);
		messageBox.setIcon(QMessageBox::Warning);
		
		int choice = messageBox.exec();
		if(choice != QMessageBox::Ok) return;
	}
	
	setDefaults(SOFIA_PARSET_EXTRAGALACTIC);
	return;
}



// --------------------------------------
// Slot to load parameter set 'Continuum'
// --------------------------------------

void SoFiA::loadParsetContinuum()
{
	if(settingsChanged)
	{
		QMessageBox messageBox(this);
		messageBox.setWindowTitle(tr("SoFiA - New Parameter File"));
		messageBox.setText(tr("<p>This action will load the default parameter set %1 and close the current file. All unsaved changes will be lost.</p><p>Do you wish to create a new parameter file?</p>").arg(QString::fromUtf8("‘Radio Continuum Image’")));
		messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
		messageBox.setDefaultButton(QMessageBox::Ok);
		messageBox.setIcon(QMessageBox::Warning);
		
		int choice = messageBox.exec();
		if(choice != QMessageBox::Ok) return;
	}
	
	setDefaults(SOFIA_PARSET_CONTINUUM);
	return;
}



// ---------------------------
// Slot to show SoFiA handbook
// ---------------------------

void SoFiA::showHandbook(const QString &page)
{
	if(page.isEmpty()) HelpBrowser::showPage("index.html");
	else HelpBrowser::showPage(page);
	
	return;
}



// --------------------------
// Slot to show About message
// --------------------------

void SoFiA::aboutSoFiA()
{
	QString messageText = tr("<h3>About SoFiA</h3><p>Version 1.2.0-beta (using Qt %1)</p><p>SoFiA, the <b>Source Finding Application</b>, is a 3D source finding pipeline designed to detect and parameterise galaxies in HI data cubes. The acronym SoFiA is based on the Greek word %2, which means wisdom.</p><p>SoFiA is free software: you can redistribute it and/or modify it under the terms of the <b>GNU General Public License</b> as published by the Free Software Foundation, either version 3 of the licence, or (at your option) any later version.</p><p>SoFiA is distributed in the hope that it will be useful, but <b>without any warranty</b>; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.</p><p>You should have received a copy of the GNU General Public License along with SoFiA. If not, see <a href=\"http://www.gnu.org/licenses/\">http://www.gnu.org/licenses/</a>.</p><p>SoFiA uses the Oxygen icon set which is licensed under version&nbsp;3 of the <a href=\"http://www.gnu.org/licenses/lgpl-3.0.txt\">GNU Lesser General Public License</a>. For more details please see the Oxygen section on <a href=\"https://techbase.kde.org/Projects/Oxygen/Licensing\">KDE TechBase</a> or the <a href=\"http://www.kde.org/\">KDE website</a>.</p><p>&copy; 2018 The SoFiA Authors</p>").arg(QString(qVersion())).arg(QString::fromUtf8("σοφία"));
	QString statusText = QString("");
	showMessage(MESSAGE_INFO, messageText, statusText);
	
	return;
}



// --------------------
// Slot to run pipeline
// --------------------

void SoFiA::runPipeline()
{
	if(pipelineProcess->state() != QProcess::NotRunning)
	{
		QString messageText = tr("<p>The pipeline is already running.</p>");
		QString statusText = tr("Pipeline already running.");
		showMessage(MESSAGE_INFO, messageText, statusText);
	}
	else
	{
		if(dockWidgetOutput->isHidden()) dockWidgetOutput->show();
		
		QString command("python");
		QStringList arguments;
		
		QString SOFIA_FULL_PATH = getenv("SOFIA_PIPELINE_PATH");
		
		// Replace name of pipeline if optically motivated source finding is requested:
		if(tabInputGroupBox2->isChecked()) SOFIA_FULL_PATH.replace("sofia_pipeline.py", "optical_find.py");
		
		// Save current parameter settings to session file:
		saveSession();
		
		// Run pipeline on session file:
		arguments << SOFIA_FULL_PATH.toUtf8().data() << SOFIA_TEMP_FILE;
		pipelineProcess->start(command, arguments);
	}
	
	return;
}



// -----------------------------------------
// Slot to read output from pipeline process
// -----------------------------------------

void SoFiA::pipelineProcessReadStd()
{
	QByteArray output = pipelineProcess->readAllStandardOutput();
	QString    outputStd(QString::fromUtf8(output));
	
	outputStd.remove(QChar('\r'));       // Get rid of carriage returns in the output
	
	if(not outputStd.isEmpty())
	{
		outputText->moveCursor(QTextCursor::End);
		outputText->setTextColor(Qt::black);
		outputText->insertPlainText(outputStd);
		outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
	}
	
	return;
}



// -----------------------------------------
// Slot to read errors from pipeline process
// -----------------------------------------

void SoFiA::pipelineProcessReadErr()
{
	QByteArray output = pipelineProcess->readAllStandardError();
	QString    outputErr(output);
	
	outputErr.remove(QChar('\r'));       // Get rid of carriage returns in the output
	
	if(not outputErr.isEmpty())
	{
		outputText->moveCursor(QTextCursor::End);
		outputText->setTextColor(Qt::red);
		outputText->insertPlainText(outputErr);
		outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
	}
	
	return;
}



// -------------------------------
// Slot to react to pipeline start
// -------------------------------

void SoFiA::pipelineProcessStarted()
{
	QString messageText("");
	QString statusText = tr("Pipeline started.");
	showMessage(MESSAGE_INFO, messageText, statusText);
	
	outputProgress->setMaximum(0);
	
	updateActions();
	
	return;
}



// --------------------------------
// Slot to react to pipeline finish
// --------------------------------

void SoFiA::pipelineProcessFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
	if(exitStatus == QProcess::NormalExit)
	{
		if(exitCode == 0)
		{
			// Pipeline finished successfully:
			QString messageText("");
			QString statusText = tr("Pipeline finished.");
			showMessage(MESSAGE_INFO, messageText, statusText);
			
			outputText->setTextColor(Qt::darkGreen);
			outputText->insertPlainText(QString("Pipeline finished with exit code %1.\n").arg(exitCode));
			outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
			
			outputProgress->setValue(100);
			
			if(not spreadsheet->isHidden()) showCatalogue();  // Reload catalogue if currently visible.
		}
		else
		{
			// Pipeline finished with error:
			QString messageText("");
			QString statusText = tr("Pipeline failed.");
			showMessage(MESSAGE_INFO, messageText, statusText);
			
			outputProgress->setValue(0);
			
			outputText->setTextColor(Qt::red);
			outputText->insertPlainText(QString("Pipeline failed with exit code %1.\n").arg(exitCode));
			outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
		}
	}
	else
	{
		// Pipeline was aborted or crashed:
		QString messageText("");
		QString statusText = tr("Pipeline aborted.");
		showMessage(MESSAGE_ERROR, messageText, statusText);
		
		outputProgress->setValue(0);
		
		outputText->setTextColor(Qt::red);
		outputText->insertPlainText(QString("Pipeline aborted with exit code %1.\n").arg(exitCode));
		outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
	}
	
	outputProgress->setMaximum(100);
	updateActions();
	
	return;
}



// ---------------------------
// Slot to cancel pipeline run
// ---------------------------

void SoFiA::pipelineProcessCancel()
{
	actionAbort->setEnabled(false);
	actionExit->setEnabled(true);
	
	if(pipelineProcess->state() == QProcess::Running)
	{
		pipelineProcess->terminate();                      // First try to terminate process.
	}
	
	if(not pipelineProcess->waitForFinished(30000))           // Give it 30 s to terminate.
	{
		if(pipelineProcess->state() == QProcess::Running)
		{
			pipelineProcess->kill();                       // If it doesn't, kill it.
		}
	}
	
	updateActions();
	
	return;
}



// ---------------------------------------
// Slot to react to pipeline process error
// ---------------------------------------

void SoFiA::pipelineProcessError(QProcess::ProcessError error)
{
	QString messageText("");
	QString statusText("");
	
	switch(error)
	{
		case QProcess::FailedToStart:
			outputText->setTextColor(Qt::red);
			outputText->insertPlainText(QString("Error: Failed to launch pipeline.\n"));
			outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
			
			messageText = tr("<p>Failed to launch pipeline.</p><p>Please ensure that the pipeline is installed on your computer and you have permission to execute it.</p>");
			statusText = tr("Failed to launch pipeline.");
			showMessage(MESSAGE_ERROR, messageText, statusText);
			break;
			
		case QProcess::Crashed:
			outputText->setTextColor(Qt::red);
			outputText->insertPlainText(QString("Error: Pipeline terminated prematurely.\n"));
			outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
			break;
			
		case QProcess::Timedout:
			outputText->setTextColor(Qt::red);
			outputText->insertPlainText(QString("Error: Pipeline timed out.\n"));
			outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
			break;
			
		case QProcess::WriteError:
			outputText->setTextColor(Qt::red);
			outputText->insertPlainText(QString("Error: Pipeline failed to receive input.\n"));
			outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
			break;
			
		case QProcess::ReadError:
			outputText->setTextColor(Qt::red);
			outputText->insertPlainText(QString("Error: Failed to receive output from pipeline.\n"));
			outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
			break;
			
		default:
			outputText->setTextColor(Qt::red);
			outputText->insertPlainText(QString("Error: An unspecified error occurred.\n"));
			outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
	}
	
	return;
}



// -----------------------------
// Function to extract file name
// -----------------------------

QString SoFiA::extractFileName(QString &extension)
{
	QString filename = (tabOutputFieldBaseName->text()).trimmed();
	QString dirname = (tabOutputFieldDirectory->text()).trimmed();
	
	QString fullFilePath = (tabInputFieldData->text()).trimmed();
	QString separator    = QString("/");                                       // WARNING: This might only work on Unix/Linux systems!
	int     separatorPos = fullFilePath.lastIndexOf(separator) + 1;
	
	// Determine file name:
	if(filename.isEmpty() or filename.contains("/") or filename.contains("\\") or filename == "." or filename == "..")
	{
		if(separatorPos > 0) filename = fullFilePath.right(fullFilePath.size() - separatorPos);
		else filename = fullFilePath;
	}
	if((filename.toLower()).endsWith(".fits") and filename.size() > 5) filename = filename.left(filename.size() - 5);
	filename.append(extension);
	
	// Determine directory name:
	if(dirname.isEmpty() and separatorPos > 0) dirname = fullFilePath.left(separatorPos);
	if(not (dirname.isEmpty() or dirname.endsWith("/"))) dirname.append("/");
	
	// Concatenate directory and file names:
	filename.prepend(dirname);
	
	return filename;
}



// -----------------------------
// Slot to show source catalogue
// -----------------------------

void SoFiA::showCatalogue()
{
	QString extension = QString("_cat.xml");
	QString filename = extractFileName(extension);
	
	// Load catalogue:
	if(spreadsheet->loadCatalog(filename))
	{
		QString messageText = tr("<p>Failed to load source catalogue:</p><p>\"%1 [.gz]\"</p>").arg(filename);
		QString statusText  = tr("Failed to load source catalogue.");
		showMessage(MESSAGE_ERROR, messageText, statusText);
	}
	else
	{
		spreadsheet->show();
		spreadsheet->raise();
	}
	
	return;
}



// ----------------------
// Slot to show data cube
// ----------------------

void SoFiA::showCube()
{
	QString filename = tabInputFieldData->text();
	
	WidgetDataViewer *widgetDataViewer = new WidgetDataViewer(filename.toUtf8().constData(), this);
	widgetDataViewer->show();
	
	return;
}



// --------------------------
// Slot to show filtered cube
// --------------------------

void SoFiA::showFilteredCube()
{
	QString extension = QString("_filtered.fits");
	QString filename = extractFileName(extension);
	
	WidgetDataViewer *widgetDataViewer = new WidgetDataViewer(filename.toUtf8().constData(), this);
	widgetDataViewer->show();
	
	return;
}



// -----------------------
// Slot to show noise cube
// -----------------------

void SoFiA::showNoiseCube()
{
	QString extension = QString("_noise.fits");
	QString filename = extractFileName(extension);
	
	WidgetDataViewer *widgetDataViewer = new WidgetDataViewer(filename.toUtf8().constData(), this);
	widgetDataViewer->show();
	
	return;
}



// ----------------------
// Slot to show mask cube
// ----------------------

void SoFiA::showMask()
{
	QString extension = QString("_mask.fits");
	QString filename = extractFileName(extension);
	
	WidgetDataViewer *widgetDataViewer = new WidgetDataViewer(filename.toUtf8().constData(), this);
	widgetDataViewer->show();
	
	return;
}



// ----------------------
// Slot to show mom 0 map
// ----------------------

void SoFiA::showMom0()
{
	QString extension = QString("_mom0.fits");
	QString filename = extractFileName(extension);
	
	WidgetDataViewer *widgetDataViewer = new WidgetDataViewer(filename.toUtf8().constData(), this);
	widgetDataViewer->show();
	
	return;
}



// ----------------------
// Slot to show mom 1 map
// ----------------------

void SoFiA::showMom1()
{
	QString extension = QString("_mom1.fits");
	QString filename = extractFileName(extension);
	
	WidgetDataViewer *widgetDataViewer = new WidgetDataViewer(filename.toUtf8().constData(), this);
	widgetDataViewer->show();
	
	return;
}



// -------------------------------
// Slot to toggle full-screen mode
// -------------------------------

void SoFiA::toggleFullScreen()
{
	if(this->isFullScreen())
	{
		this->showNormal();
		actionFullScreen->setChecked(false);
	}
	else
	{
		this->showFullScreen();
		actionFullScreen->setChecked(true);
	}
	
	return;
}



// ------------------------------
// Slot to toggle toolbar setting
// ------------------------------

void SoFiA::toggleToolbar()
{
	settingsToolbar = actionToolbar->isChecked();
	toolBar->setVisible(settingsToolbar);
	
	return;
}



// -------------------------------
// Slot to toggle pipeline setting
// -------------------------------

void SoFiA::togglePipeline(bool state)
{
	// React to external reasons for message window to (dis)appear
	// (but ignore programme shutdown as cause):
	if(not shutdownInitiated)
	{
		settingsPipeline = state;
		actionPipeline->setChecked(settingsPipeline);
	}
	
	return;
}

void SoFiA::togglePipeline()
{
	// React to regular change of settingsPipeline value:
	settingsPipeline = actionPipeline->isChecked();
	dockWidgetOutput->setVisible(settingsPipeline);
	
	return;
}



// -----------------------------------
// Slot to toggle save-on-exit setting
// -----------------------------------

void SoFiA::toggleSaveOnExit()
{
	settingsSession = actionSaveOnExit->isChecked();
	
	return;
}



// ------------------------------------------------------
// Slot to create pop-up warning message about positivity
// ------------------------------------------------------

void SoFiA::printPositivityWarning(bool checked)
{
	if(checked)
	{
		QString messageText = tr("<p>Enabling the positivity parameter is <strong>dangerous</strong> and will render some of SoFiA&rsquo;s most powerful algorithms useless, including mask optimisation and reliability calculation.</p><p>Only use this option if you are fully aware of the risks and consequences!</p>");
		QString statusText  = QString("");
		showMessage(MESSAGE_WARNING, messageText, statusText);
	}
	
	return;
}



// ------------------------------------
// Function to update window title text
// ------------------------------------

void SoFiA::updateWindowTitle()
{
	if(currentFileName.isEmpty())
	{
		if(settingsChanged) this->setWindowTitle(tr("SoFiA (modified)"));
		else this->setWindowTitle(tr("SoFiA"));
	}
	else
	{
		if(settingsChanged) this->setWindowTitle(tr("SoFiA - %1 (modified)").arg(currentFileName.section('/', -1)));
		else this->setWindowTitle(tr("SoFiA - %1").arg(currentFileName.section('/', -1)));
	}
	
	actionSave->setEnabled(settingsChanged);
	
	return;
}



// --------------------------------------------
// Function to create and set up user interface
// --------------------------------------------

void SoFiA::createInterface()
{
	// Load icons
	// ----------
	
	iconSoFiA.addFile(QString(":/icons/32/SoFiA.png"), QSize(32, 32));
	iconSoFiA.addFile(QString(":/icons/22/SoFiA.png"), QSize(22, 22));
	iconSoFiA.addFile(QString(":/icons/16/SoFiA.png"), QSize(16, 16));
	
	iconWhatsThis.addFile(QString(":/icons/22/whats-this.png"), QSize(22, 22));
	iconWhatsThis.addFile(QString(":/icons/16/whats-this.png"), QSize(16, 16));
	
	iconDocumentNew.addFile(QString(":/icons/22/document-new.png"), QSize(22, 22));
	iconDocumentNew.addFile(QString(":/icons/16/document-new.png"), QSize(16, 16));
	iconDocumentNew     = QIcon::fromTheme("document-new", iconDocumentNew);
	
	iconDocumentOpen.addFile(QString(":/icons/22/document-open.png"), QSize(22, 22));
	iconDocumentOpen.addFile(QString(":/icons/16/document-open.png"), QSize(16, 16));
	iconDocumentOpen    = QIcon::fromTheme("document-open", iconDocumentOpen);
	
	iconDocumentPreview.addFile(QString(":/icons/22/document-preview.png"), QSize(22, 22));
	iconDocumentPreview.addFile(QString(":/icons/16/document-preview.png"), QSize(16, 16));
	iconDocumentPreview = QIcon::fromTheme("document-preview", iconDocumentPreview);
	
	iconDocumentSave.addFile(QString(":/icons/22/document-save.png"), QSize(22, 22));
	iconDocumentSave.addFile(QString(":/icons/16/document-save.png"), QSize(16, 16));
	iconDocumentSave    = QIcon::fromTheme("document-save", iconDocumentSave);
	
	iconDocumentSaveAs.addFile(QString(":/icons/22/document-save-as.png"), QSize(22, 22));
	iconDocumentSaveAs.addFile(QString(":/icons/16/document-save-as.png"), QSize(16, 16));
	iconDocumentSaveAs  = QIcon::fromTheme("document-save-as", iconDocumentSaveAs);
	
	iconApplicationExit.addFile(QString(":/icons/22/application-exit.png"), QSize(22, 22));
	iconApplicationExit.addFile(QString(":/icons/16/application-exit.png"), QSize(16, 16));
	iconApplicationExit = QIcon::fromTheme("application-exit", iconApplicationExit);
	
	iconDialogOkApply.addFile(QString(":/icons/22/dialog-ok-apply.png"), QSize(22, 22));
	iconDialogOkApply.addFile(QString(":/icons/16/dialog-ok-apply.png"), QSize(16, 16));
	iconDialogOkApply   = QIcon::fromTheme("dialog-ok-apply", iconDialogOkApply);
	
	iconDialogCancel.addFile(QString(":/icons/22/dialog-cancel.png"), QSize(22, 22));
	iconDialogCancel.addFile(QString(":/icons/16/dialog-cancel.png"), QSize(16, 16));
	iconDialogCancel    = QIcon::fromTheme("dialog-cancel", iconDialogCancel);
	
	iconDialogClose.addFile(QString(":/icons/22/dialog-close.png"), QSize(22, 22));
	iconDialogClose.addFile(QString(":/icons/16/dialog-close.png"), QSize(16, 16));
	iconDialogClose     = QIcon::fromTheme("dialog-close", iconDialogClose);
	
	iconGoPreviousView.addFile(QString(":/icons/22/go-previous-view.png"), QSize(22, 22));
	iconGoPreviousView.addFile(QString(":/icons/16/go-previous-view.png"), QSize(16, 16));
	iconGoPreviousView  = QIcon::fromTheme("go-previous-view", iconGoPreviousView);
	
	iconGoNextView.addFile(QString(":/icons/22/go-next-view.png"), QSize(22, 22));
	iconGoNextView.addFile(QString(":/icons/16/go-next-view.png"), QSize(16, 16));
	iconGoNextView      = QIcon::fromTheme("go-next-view", iconGoNextView);
	
	iconEditClearList.addFile(QString(":/icons/22/edit-clear-list.png"), QSize(22, 22));
	iconEditClearList.addFile(QString(":/icons/16/edit-clear-list.png"), QSize(16, 16));
	iconEditClearList   = QIcon::fromTheme("edit-clear-list", iconEditClearList);
	
	iconFullScreen.addFile(QString(":/icons/22/view-fullscreen.png"), QSize(22, 22));
	iconFullScreen.addFile(QString(":/icons/16/view-fullscreen.png"), QSize(16, 16));
	iconFullScreen      = QIcon::fromTheme("view-fullscreen", iconFullScreen);
	
	iconHelpContents.addFile(QString(":/icons/22/help-contents.png"), QSize(22, 22));
	iconHelpContents.addFile(QString(":/icons/16/help-contents.png"), QSize(16, 16));
	iconHelpContents    = QIcon::fromTheme("help-contents", iconHelpContents);
	
	iconHelpAbout.addFile(QString(":/icons/22/help-about.png"), QSize(22, 22));
	iconHelpAbout.addFile(QString(":/icons/16/help-about.png"), QSize(16, 16));
	iconHelpAbout       = QIcon::fromTheme("help-about", iconHelpAbout);
	
	iconTaskComplete.addFile(QString(":/icons/16/task-complete.png"), QSize(16, 16));
	iconTaskComplete.addFile(QString(":/icons/22/task-complete.png"), QSize(22, 22));
	iconTaskComplete    = QIcon::fromTheme("task-complete", iconTaskComplete);
	
	iconTaskReject.addFile(QString(":/icons/16/task-reject.png"), QSize(16, 16));
	iconTaskReject.addFile(QString(":/icons/22/task-reject.png"), QSize(22, 22));
	iconTaskReject      = QIcon::fromTheme("task-reject", iconTaskReject);
	
	iconFolderImage.addFile(QString(":/icons/16/folder-image.png"), QSize(16, 16));
	iconFolderImage.addFile(QString(":/icons/22/folder-image.png"), QSize(22, 22));
	iconFolderImage     = QIcon::fromTheme("folder-image", iconFolderImage);
	
	// Create main widget that contains everything else
	// ------------------------------------------------
	
	widgetMain = new QWidget(this);
	
	// Set up tabs
	// -----------
	
	tabs = new QTabWidget(widgetMain);
	
	tabInput           = new QWidget(tabs);
	tabInFilter        = new QWidget(tabs);
	tabSourceFinding   = new QWidget(tabs);
	tabMerging         = new QWidget(tabs);
	tabParametrisation = new QWidget(tabs);
	tabOutput          = new QWidget(tabs);
	
	tabs->addTab(tabInput,           tr("Input"));
	tabs->addTab(tabInFilter,        tr("Filtering"));
	tabs->addTab(tabSourceFinding,   tr("Source Finding"));
	tabs->addTab(tabMerging,         tr("Merging"));
	tabs->addTab(tabParametrisation, tr("Parameterisation"));
	tabs->addTab(tabOutput,          tr("Output"));
	
	tabs->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Minimum);
	tabs->setUsesScrollButtons(false);
	
	// Note that additional spaces at the end of the display text for checkboxes 
	// and radio buttons are a workaround for rendering issues under Mac OS X 
	// whereby the last character of the text was partially cut off.
	
	// Set up input tab
	// ----------------
	
	toolBoxIP = new QToolBox(tabInput);
	
	tabInputLayout = new QVBoxLayout();
	
	// input files
	tabInputGroupBox1 = new QGroupBox(toolBoxIP);
	tabInputForm1 = new QFormLayout();
	
	tabInputWidgetData = new QWidget(tabInputGroupBox1);
	tabInputLayoutData = new QHBoxLayout();
	tabInputFieldData  = new QLineEdit(tabInputWidgetData);
	tabInputFieldData->setObjectName("import.inFile");
	connect(tabInputFieldData, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInputFieldData, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabInputButtonData = new QPushButton(tr("Select..."), tabInputWidgetData);
	connect(tabInputButtonData, SIGNAL(clicked()), this, SLOT(selectInputDataFile()));
	tabInputButtonData->setIcon(iconDocumentOpen);
	tabInputLayoutData->addWidget(tabInputFieldData);
	tabInputLayoutData->addWidget(tabInputButtonData);
	tabInputLayoutData->setContentsMargins(0, 0, 0, 0);
	tabInputWidgetData->setLayout(tabInputLayoutData);
	
	tabInputWidgetMask = new QWidget(tabInputGroupBox1);
	tabInputLayoutMask = new QHBoxLayout();
	tabInputFieldMask  = new QLineEdit(tabInputWidgetMask);
	tabInputFieldMask->setObjectName("import.maskFile");
	connect(tabInputFieldMask, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabInputButtonMask = new QPushButton(tr("Select..."), tabInputWidgetMask);
	connect(tabInputButtonMask, SIGNAL(clicked()), this, SLOT(selectInputMaskFile()));
	tabInputButtonMask->setIcon(iconDocumentOpen);
	tabInputLayoutMask->addWidget(tabInputFieldMask);
	tabInputLayoutMask->addWidget(tabInputButtonMask);
	tabInputLayoutMask->setContentsMargins(0, 0, 0, 0);
	tabInputWidgetMask->setLayout(tabInputLayoutMask);
	
	tabInputWidgetWeights = new QWidget(tabInputGroupBox1);
	tabInputLayoutWeights = new QHBoxLayout();
	tabInputFieldWeights  = new QLineEdit(tabInputWidgetWeights);
	tabInputFieldWeights->setObjectName("import.weightsFile");
	connect(tabInputFieldWeights, SIGNAL(textChanged(const QString &)), this, SLOT(updateFields()));
	connect(tabInputFieldWeights, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabInputButtonWeights = new QPushButton(tr("Select..."), tabInputWidgetWeights);
	connect(tabInputButtonWeights, SIGNAL(clicked()), this, SLOT(selectInputWeightsFile()));
	tabInputButtonWeights->setIcon(iconDocumentOpen);
	tabInputLayoutWeights->addWidget(tabInputFieldWeights);
	tabInputLayoutWeights->addWidget(tabInputButtonWeights);
	tabInputLayoutWeights->setContentsMargins(0, 0, 0, 0);
	tabInputWidgetWeights->setLayout(tabInputLayoutWeights);
	
	tabInputFieldWeightsFunction = new QLineEdit(tabInputGroupBox1);
	tabInputFieldWeightsFunction->setObjectName("import.weightsFunction");
	connect(tabInputFieldWeightsFunction, SIGNAL(textChanged(const QString &)), this, SLOT(updateFields()));
	connect(tabInputFieldWeightsFunction, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabInputForm1->addRow(tr("Data cube:"), tabInputWidgetData);
	tabInputForm1->addRow(tr("Mask cube:"), tabInputWidgetMask);
	tabInputForm1->addRow(tr("Weights cube:"), tabInputWidgetWeights);
	tabInputForm1->addRow(tr("Weights function:"), tabInputFieldWeightsFunction);
	tabInputForm1->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
	tabInputGroupBox1->setLayout(tabInputForm1);
	
	// subcube
	tabInputGroupBox4 = new QGroupBox(tr("Enable"), toolBoxIP);
	tabInputGroupBox4->setObjectName("steps.doSubcube");
	tabInputGroupBox4->setCheckable(true);
	tabInputGroupBox4->setChecked(false);
	connect(tabInputGroupBox4, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabInputGroupBox4, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabInputForm4 = new QFormLayout();
	
	tabInputFieldSubcube = new QLineEdit(tabInputGroupBox4);
	tabInputFieldSubcube->setObjectName("import.subcube");
	connect(tabInputFieldSubcube, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabInputFieldSubcubeMode = new QComboBox(tabInputGroupBox4);
	tabInputFieldSubcubeMode->setObjectName("import.subcubeMode");
	tabInputFieldSubcubeMode->addItem(tr("Pixels"), QVariant(QString("pixel")));
	tabInputFieldSubcubeMode->addItem(tr("World coordinates"), QVariant(QString("world")));
	connect(tabInputFieldSubcubeMode, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabInputForm4->addRow(tr("Range:"), tabInputFieldSubcube);
	tabInputForm4->addRow(tr("Mode:"), tabInputFieldSubcubeMode);
	tabInputForm4->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
	tabInputGroupBox4->setLayout(tabInputForm4);
	
	// flagging
	tabInputGroupBox3 = new QGroupBox(tr("Enable"), toolBoxIP);
	tabInputGroupBox3->setObjectName("steps.doFlag");
	tabInputGroupBox3->setCheckable(true);
	tabInputGroupBox3->setChecked(false);
	connect(tabInputGroupBox3, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabInputGroupBox3, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabInputForm3 = new QFormLayout();
	
	tabInputWidgetFlagCube = new QWidget(tabInputGroupBox3);
	tabInputLayoutFlagCube = new QHBoxLayout();
	tabInputFieldFlagCube  = new QLineEdit(tabInputWidgetFlagCube);
	tabInputFieldFlagCube->setObjectName("flag.file");
	connect(tabInputFieldFlagCube, SIGNAL(textChanged(const QString &)), this, SLOT(updateFields()));
	connect(tabInputFieldFlagCube, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabInputButtonFlagCube = new QPushButton(tr("Select..."), tabInputWidgetFlagCube);
	connect(tabInputButtonFlagCube, SIGNAL(clicked()), this, SLOT(selectInputFlagFile()));
	tabInputButtonFlagCube->setIcon(iconDocumentOpen);
	tabInputLayoutFlagCube->addWidget(tabInputFieldFlagCube);
	tabInputLayoutFlagCube->addWidget(tabInputButtonFlagCube);
	tabInputLayoutFlagCube->setContentsMargins(0, 0, 0, 0);
	tabInputWidgetFlagCube->setLayout(tabInputLayoutFlagCube);
	
	tabInputFieldFlags = new QLineEdit(tabInputGroupBox3);
	tabInputFieldFlags->setObjectName("flag.regions");
	connect(tabInputFieldFlags, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabInputForm3->addRow(tr("Flag cube:"), tabInputWidgetFlagCube);
	tabInputForm3->addRow(tr("Range:"), tabInputFieldFlags);
	tabInputForm3->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
	tabInputGroupBox3->setLayout(tabInputForm3);
	
	// optical source finding
	tabInputGroupBox2 = new QGroupBox(tr("Enable"), toolBoxIP);
	tabInputGroupBox2->setObjectName("steps.doOptical");
	tabInputGroupBox2->setCheckable(true);
	tabInputGroupBox2->setChecked(false);
	connect(tabInputGroupBox2, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabInputGroupBox2, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabInputForm2 = new QFormLayout();
	
	tabInputWidgetCatalog = new QWidget(tabInputGroupBox2);
	tabInputLayoutCatalog = new QHBoxLayout();
	tabInputFieldCatalog = new QLineEdit(tabInputWidgetCatalog);
	tabInputFieldCatalog->setObjectName("optical.sourceCatalogue");
	connect(tabInputFieldCatalog, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabInputButtonCatalog = new QPushButton(tr("Select..."), tabInputWidgetCatalog);
	connect(tabInputButtonCatalog, SIGNAL(clicked()), this, SLOT(selectOpticalCatalogFile()));
	tabInputButtonCatalog->setIcon(iconDocumentOpen);
	tabInputLayoutCatalog->addWidget(tabInputFieldCatalog);
	tabInputLayoutCatalog->addWidget(tabInputButtonCatalog);
	tabInputLayoutCatalog->setContentsMargins(0, 0, 0, 0);
	tabInputWidgetCatalog->setLayout(tabInputLayoutCatalog);
	
	tabInputFieldSpatialSize = new QLineEdit(tabInputGroupBox2);
	tabInputFieldSpatialSize->setObjectName("optical.spatSize");
	connect(tabInputFieldSpatialSize, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInputFieldSpatialSize, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabInputFieldSpectralSize = new QLineEdit(tabInputGroupBox2);
	tabInputFieldSpectralSize->setObjectName("optical.specSize");
	connect(tabInputFieldSpectralSize, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInputFieldSpectralSize, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabInputFieldMultiCat = new QCheckBox(tr("Create separate output catalogues"), tabInputGroupBox2);
	tabInputFieldMultiCat->setObjectName("optical.storeMultiCat");
	tabInputFieldMultiCat->setChecked(false);
	connect(tabInputFieldMultiCat, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabInputForm2->addRow(tr("Catalogue:"), tabInputWidgetCatalog);
	tabInputForm2->addRow(tr("Spatial size:"), tabInputFieldSpatialSize);
	tabInputForm2->addRow(tr("Spectral size:"), tabInputFieldSpectralSize);
	tabInputForm2->addRow(tr("Output:"), tabInputFieldMultiCat);
	tabInputForm2->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
	tabInputGroupBox2->setLayout(tabInputForm2);
	
	// controls
	tabInputButtonNext = new QPushButton(tr("Next"), tabInput);
	tabInputButtonNext->setIcon(iconGoNextView);
	connect(tabInputButtonNext, SIGNAL(clicked()), this, SLOT(displayNextTab()));
	tabInputLayoutControls = new QHBoxLayout();
	tabInputLayoutControls->setContentsMargins(0, 0, 0, 0);
	tabInputLayoutControls->setSpacing(0);
	tabInputLayoutControls->addStretch();
	tabInputLayoutControls->addWidget(tabInputButtonNext);
	tabInputWidgetControls = new QWidget(tabInput);
	tabInputWidgetControls->setLayout(tabInputLayoutControls);
	
	toolBoxIP->addItem(tabInputGroupBox1, iconTaskReject, tr("Input Data Products"));
	toolBoxIP->addItem(tabInputGroupBox4, iconTaskReject, tr("Subcube"));
	toolBoxIP->addItem(tabInputGroupBox3, iconTaskReject, tr("Flagging"));
	toolBoxIP->addItem(tabInputGroupBox2, iconTaskReject, tr("Catalogue-based Source Finding"));
	
	tabInputLayout->addWidget(toolBoxIP);
	tabInputLayout->addStretch();
	tabInputLayout->addWidget(tabInputWidgetControls);
	tabInput->setLayout(tabInputLayout);
	
	
	
	// Set up input filter tab
	// -----------------------
	
	toolBoxIF = new QToolBox(tabInFilter);
	
	tabInFilterLayout = new QVBoxLayout();
	
	tabInFilterGroupBox1 = new QGroupBox(tr("Enable"), toolBoxIF);
	tabInFilterGroupBox1->setObjectName("steps.doSmooth");
	tabInFilterGroupBox1->setCheckable(true);
	tabInFilterGroupBox1->setChecked(false);
	connect(tabInFilterGroupBox1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabInFilterGroupBox1, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	// smoothing
	tabInFilterForm1 = new QFormLayout();
	
	tabInFilterFieldKernel = new QComboBox(tabInFilterGroupBox1);
	tabInFilterFieldKernel->setObjectName("smooth.kernel");
	tabInFilterFieldKernel->addItem(tr("Gaussian"), QVariant(QString("gaussian")));
	tabInFilterFieldKernel->addItem(tr("Boxcar"), QVariant(QString("boxcar")));
	tabInFilterFieldKernel->addItem(tr("Median"), QVariant(QString("median")));
	connect(tabInFilterFieldKernel, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterFieldBorder = new QComboBox(tabInFilterGroupBox1);
	tabInFilterFieldBorder->setObjectName("smooth.edgeMode");
	tabInFilterFieldBorder->addItem(tr("Constant"), QVariant(QString("constant")));
	tabInFilterFieldBorder->addItem(tr("Reflect"), QVariant(QString("reflect")));
	tabInFilterFieldBorder->addItem(tr("Mirror"), QVariant(QString("mirror")));
	tabInFilterFieldBorder->addItem(tr("Nearest"), QVariant(QString("nearest")));
	tabInFilterFieldBorder->addItem(tr("Wrap"), QVariant(QString("wrap")));
	connect(tabInFilterFieldBorder, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterFieldSmoothingSpatialLon  = new QLineEdit(tabInFilterGroupBox1);
	tabInFilterFieldSmoothingSpatialLon->setObjectName("smooth.kernelX");
	tabInFilterFieldSmoothingSpatialLon->setMaxLength(10);
	tabInFilterFieldSmoothingSpatialLon->setMaximumWidth(100);
	connect(tabInFilterFieldSmoothingSpatialLon, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInFilterFieldSmoothingSpatialLon, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabInFilterFieldSmoothingSpatialLat  = new QLineEdit(tabInFilterGroupBox1);
	tabInFilterFieldSmoothingSpatialLat->setObjectName("smooth.kernelY");
	tabInFilterFieldSmoothingSpatialLat->setMaxLength(10);
	tabInFilterFieldSmoothingSpatialLat->setMaximumWidth(100);
	connect(tabInFilterFieldSmoothingSpatialLat, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInFilterFieldSmoothingSpatialLat, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabInFilterFieldSmoothingSpectral = new QLineEdit(tabInFilterGroupBox1);
	tabInFilterFieldSmoothingSpectral->setObjectName("smooth.kernelZ");
	tabInFilterFieldSmoothingSpectral->setMaxLength(10);
	tabInFilterFieldSmoothingSpectral->setMaximumWidth(100);
	connect(tabInFilterFieldSmoothingSpectral, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInFilterFieldSmoothingSpectral, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabInFilterForm1->addRow(tr("Kernel:"), tabInFilterFieldKernel);
	tabInFilterForm1->addRow(tr("Edge:"), tabInFilterFieldBorder);
	tabInFilterForm1->addRow(tr("Scale X:"), tabInFilterFieldSmoothingSpatialLon);
	tabInFilterForm1->addRow(tr("Scale Y:"), tabInFilterFieldSmoothingSpatialLat);
	tabInFilterForm1->addRow(tr("Scale Z:"), tabInFilterFieldSmoothingSpectral);
	
	
	
	tabInFilterGroupBox2 = new QGroupBox(tr("Enable"), toolBoxIF);
	tabInFilterGroupBox2->setObjectName("steps.doScaleNoise");
	tabInFilterGroupBox2->setCheckable(true);
	tabInFilterGroupBox2->setChecked(false);
	connect(tabInFilterGroupBox2, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabInFilterGroupBox2, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	// noise scaling
	tabInFilterForm2 = new QFormLayout();
	
	tabInFilterFieldMethod = new QComboBox(tabInFilterGroupBox2);
	tabInFilterFieldMethod->setObjectName("scaleNoise.method");
	tabInFilterFieldMethod->addItem(tr("Global"), QVariant(QString("global")));
	tabInFilterFieldMethod->addItem(tr("Local"), QVariant(QString("local")));
	connect(tabInFilterFieldMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFields()));
	connect(tabInFilterFieldMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterFieldStatistic = new QComboBox(tabInFilterGroupBox2);
	tabInFilterFieldStatistic->setObjectName("scaleNoise.statistic");
	tabInFilterFieldStatistic->addItem(tr("Gaussian fit to negative fluxes"), QVariant(QString("negative")));
	tabInFilterFieldStatistic->addItem(tr("Gaussian fit to flux histogram"), QVariant(QString("gauss")));
	tabInFilterFieldStatistic->addItem(tr("Median absolute deviation"), QVariant(QString("mad")));
	tabInFilterFieldStatistic->addItem(tr("Standard deviation"), QVariant(QString("std")));
	connect(tabInFilterFieldStatistic, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFields()));
	connect(tabInFilterFieldStatistic, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterFieldFluxRange = new QComboBox(tabInFilterGroupBox2);
	tabInFilterFieldFluxRange->setObjectName("scaleNoise.fluxRange");
	tabInFilterFieldFluxRange->addItem(tr("Negative"), QVariant(QString("negative")));
	tabInFilterFieldFluxRange->addItem(tr("Positive"), QVariant(QString("positive")));
	tabInFilterFieldFluxRange->addItem(tr("All"), QVariant(QString("all")));
	connect(tabInFilterFieldFluxRange, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterWidgetScaleXYZ = new QWidget(tabInFilterGroupBox2);
	tabInFilterLayoutScaleXYZ = new QHBoxLayout();
	tabInFilterFieldScaleX = new QCheckBox(tr("X"), tabInFilterWidgetScaleXYZ);
	tabInFilterFieldScaleX->setObjectName("scaleNoise.scaleX");
	tabInFilterFieldScaleX->setChecked(false);
	connect(tabInFilterFieldScaleX, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabInFilterFieldScaleY = new QCheckBox(tr("Y"), tabInFilterWidgetScaleXYZ);
	tabInFilterFieldScaleY->setObjectName("scaleNoise.scaleY");
	tabInFilterFieldScaleY->setChecked(false);
	connect(tabInFilterFieldScaleY, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabInFilterFieldScaleZ = new QCheckBox(tr("Z"), tabInFilterWidgetScaleXYZ);
	tabInFilterFieldScaleZ->setObjectName("scaleNoise.scaleZ");
	tabInFilterFieldScaleZ->setChecked(true);
	connect(tabInFilterFieldScaleZ, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabInFilterLayoutScaleXYZ->setContentsMargins(0, 0, 0, 0);
	tabInFilterLayoutScaleXYZ->addWidget(tabInFilterFieldScaleX);
	tabInFilterLayoutScaleXYZ->addWidget(tabInFilterFieldScaleY);
	tabInFilterLayoutScaleXYZ->addWidget(tabInFilterFieldScaleZ);
	tabInFilterWidgetScaleXYZ->setLayout(tabInFilterLayoutScaleXYZ);
	
	tabInFilterWidgetEdge = new QWidget(tabInFilterGroupBox2);
	tabInFilterLayoutEdge = new QHBoxLayout();
	tabInFilterFieldEdgeX  = new QSpinBox(tabInFilterWidgetEdge);
	tabInFilterFieldEdgeX->setObjectName("scaleNoise.edgeX");
	tabInFilterFieldEdgeX->setMaximumWidth(100);
	tabInFilterFieldEdgeX->setMinimum(0);
	tabInFilterFieldEdgeX->setMaximum(100);
	connect(tabInFilterFieldEdgeX, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabInFilterLabelEdgeX = new QLabel(tr("X"), tabInFilterWidgetEdge);
	tabInFilterFieldEdgeY = new QSpinBox(tabInFilterWidgetEdge);
	tabInFilterFieldEdgeY->setObjectName("scaleNoise.edgeY");
	tabInFilterFieldEdgeY->setMaximumWidth(100);
	tabInFilterFieldEdgeY->setMinimum(0);
	tabInFilterFieldEdgeY->setMaximum(100);
	connect(tabInFilterFieldEdgeY, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabInFilterLabelEdgeY = new QLabel(tr("Y"), tabInFilterWidgetEdge);
	tabInFilterFieldEdgeZ = new QSpinBox(tabInFilterWidgetEdge);
	tabInFilterFieldEdgeZ->setObjectName("scaleNoise.edgeZ");
	tabInFilterFieldEdgeZ->setMaximumWidth(100);
	tabInFilterFieldEdgeZ->setMinimum(0);
	tabInFilterFieldEdgeZ->setMaximum(100);
	connect(tabInFilterFieldEdgeZ, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabInFilterLabelEdgeZ = new QLabel(tr("Z"), tabInFilterWidgetEdge);
	tabInFilterLayoutEdge->addWidget(tabInFilterFieldEdgeX);
	tabInFilterLayoutEdge->addWidget(tabInFilterLabelEdgeX);
	tabInFilterLayoutEdge->addSpacing(10);
	tabInFilterLayoutEdge->addWidget(tabInFilterFieldEdgeY);
	tabInFilterLayoutEdge->addWidget(tabInFilterLabelEdgeY);
	tabInFilterLayoutEdge->addSpacing(10);
	tabInFilterLayoutEdge->addWidget(tabInFilterFieldEdgeZ);
	tabInFilterLayoutEdge->addWidget(tabInFilterLabelEdgeZ);
	tabInFilterLayoutEdge->setContentsMargins(0, 0, 0, 0);
	tabInFilterLayoutEdge->setSpacing(5);
	tabInFilterWidgetEdge->setLayout(tabInFilterLayoutEdge);
	
	/*tabInFilterWidgetGrid = new QWidget(tabInFilterGroupBox2);
	tabInFilterLayoutGrid = new QHBoxLayout();
	tabInFilterFieldGridSpatial = new QSpinBox(tabInFilterWidgetGrid);
	tabInFilterFieldGridSpatial->setObjectName("scaleNoise.gridSpatial");
	tabInFilterFieldGridSpatial->setMaximumWidth(100);
	tabInFilterFieldGridSpatial->setRange(2, 100);
	tabInFilterFieldGridSpatial->setSingleStep(2);
	connect(tabInFilterFieldGridSpatial, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInFilterFieldGridSpatial, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabInFilterLabelGridSpatial = new QLabel(tr("spatial"), tabInFilterWidgetGrid);
	tabInFilterFieldGridSpectral = new QSpinBox(tabInFilterWidgetGrid);
	tabInFilterFieldGridSpectral->setObjectName("scaleNoise.gridSpectral");
	tabInFilterFieldGridSpectral->setMaximumWidth(100);
	tabInFilterFieldGridSpectral->setRange(2, 100);
	tabInFilterFieldGridSpectral->setSingleStep(2);
	connect(tabInFilterFieldGridSpectral, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInFilterFieldGridSpectral, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabInFilterLabelGridSpectral = new QLabel(tr("spectral"), tabInFilterWidgetGrid);
	tabInFilterLayoutGrid->addWidget(tabInFilterFieldGridSpatial);
	tabInFilterLayoutGrid->addWidget(tabInFilterLabelGridSpatial);
	tabInFilterLayoutGrid->addSpacing(10);
	tabInFilterLayoutGrid->addWidget(tabInFilterFieldGridSpectral);
	tabInFilterLayoutGrid->addWidget(tabInFilterLabelGridSpectral);
	tabInFilterLayoutGrid->setContentsMargins(0, 0, 0, 0);
	tabInFilterLayoutGrid->setSpacing(5);
	tabInFilterWidgetGrid->setLayout(tabInFilterLayoutGrid);*/
	
	tabInFilterWidgetWindow = new QWidget(tabInFilterGroupBox2);
	tabInFilterLayoutWindow = new QHBoxLayout();
	tabInFilterFieldWindowSpatial = new QSpinBox(tabInFilterWidgetWindow);
	tabInFilterFieldWindowSpatial->setObjectName("scaleNoise.windowSpatial");
	tabInFilterFieldWindowSpatial->setMaximumWidth(100);
	tabInFilterFieldWindowSpatial->setRange(1, 999);
	tabInFilterFieldWindowSpatial->setSingleStep(2);
	connect(tabInFilterFieldWindowSpatial, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInFilterFieldWindowSpatial, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabInFilterLabelWindowSpatial = new QLabel(tr("spatial"), tabInFilterWidgetWindow);
	tabInFilterFieldWindowSpectral = new QSpinBox(tabInFilterWidgetWindow);
	tabInFilterFieldWindowSpectral->setObjectName("scaleNoise.windowSpectral");
	tabInFilterFieldWindowSpectral->setMaximumWidth(100);
	tabInFilterFieldWindowSpectral->setRange(1, 999);
	tabInFilterFieldWindowSpectral->setSingleStep(2);
	connect(tabInFilterFieldWindowSpectral, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInFilterFieldWindowSpectral, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabInFilterLabelWindowSpectral = new QLabel(tr("spectral"), tabInFilterWidgetWindow);
	tabInFilterLayoutWindow->addWidget(tabInFilterFieldWindowSpatial);
	tabInFilterLayoutWindow->addWidget(tabInFilterLabelWindowSpatial);
	tabInFilterLayoutWindow->addSpacing(10);
	tabInFilterLayoutWindow->addWidget(tabInFilterFieldWindowSpectral);
	tabInFilterLayoutWindow->addWidget(tabInFilterLabelWindowSpectral);
	tabInFilterLayoutWindow->setContentsMargins(0, 0, 0, 0);
	tabInFilterLayoutWindow->setSpacing(5);
	tabInFilterWidgetWindow->setLayout(tabInFilterLayoutWindow);
	
	tabInFilterFieldInterpolation = new QComboBox(tabInFilterGroupBox2);
	tabInFilterFieldInterpolation->setObjectName("scaleNoise.interpolation");
	tabInFilterFieldInterpolation->addItem(tr("None"), QVariant(QString("none")));
	tabInFilterFieldInterpolation->addItem(tr("Linear"), QVariant(QString("linear")));
	tabInFilterFieldInterpolation->addItem(tr("Cubic"), QVariant(QString("cubic")));
	connect(tabInFilterFieldInterpolation, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterSeparator1 = new QFrame(tabInFilterGroupBox2);
	tabInFilterSeparator1->setFrameShape(QFrame::HLine);
	tabInFilterSeparator1->setFrameShadow(QFrame::Sunken);
	tabInFilterSeparator2 = new QFrame(tabInFilterGroupBox2);
	tabInFilterSeparator2->setFrameShape(QFrame::HLine);
	tabInFilterSeparator2->setFrameShadow(QFrame::Sunken);
	
	tabInFilterForm2->addRow(tr("Method:"), tabInFilterFieldMethod);
	tabInFilterForm2->addRow(tr("Statistic:"), tabInFilterFieldStatistic);
	tabInFilterForm2->addRow(tr("Flux range:"), tabInFilterFieldFluxRange);
	tabInFilterForm2->addRow(tabInFilterSeparator1);
	tabInFilterForm2->addRow(tr("Dimensions:"), tabInFilterWidgetScaleXYZ);
	tabInFilterForm2->addRow(tr("Edge size:"), tabInFilterWidgetEdge);
	tabInFilterForm2->addRow(tabInFilterSeparator2);
	//tabInFilterForm2->addRow(tr("Grid size:"), tabInFilterWidgetGrid);
	tabInFilterForm2->addRow(tr("Window size:"), tabInFilterWidgetWindow);
	tabInFilterForm2->addRow(tr("Interpolation:"), tabInFilterFieldInterpolation);
	
	tabInFilterGroupBox3 = new QGroupBox(tr("Enable"), toolBoxIF);
	tabInFilterGroupBox3->setObjectName("steps.doWavelet");
	tabInFilterGroupBox3->setCheckable(true);
	tabInFilterGroupBox3->setChecked(false);
	connect(tabInFilterGroupBox3, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabInFilterGroupBox3, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	// 2D-1D wavelet filter
	tabInFilterForm3 = new QFormLayout();
	
	tabInFilterField2d1dThreshold = new QLineEdit(tabInFilterGroupBox3);
	tabInFilterField2d1dThreshold->setObjectName("wavelet.threshold");
	tabInFilterField2d1dThreshold->setMaximumWidth(100);
	tabInFilterField2d1dThreshold->setMaxLength(10);
	connect(tabInFilterField2d1dThreshold, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabInFilterField2d1dThreshold, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabInFilterField2d1dIterations = new QSpinBox(tabInFilterGroupBox3);
	tabInFilterField2d1dIterations->setObjectName("wavelet.iterations");
	tabInFilterField2d1dIterations->setMaximumWidth(100);
	tabInFilterField2d1dIterations->setMinimum(1);
	tabInFilterField2d1dIterations->setMaximum(50);
	connect(tabInFilterField2d1dIterations, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterField2d1dScaleXY = new QSpinBox(tabInFilterGroupBox3);
	tabInFilterField2d1dScaleXY->setObjectName("wavelet.scaleXY");
	tabInFilterField2d1dScaleXY->setMaximumWidth(100);
	tabInFilterField2d1dScaleXY->setMinimum(-1);
	tabInFilterField2d1dScaleXY->setMaximum(50);
	connect(tabInFilterField2d1dScaleXY, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterField2d1dScaleZ = new QSpinBox(tabInFilterGroupBox3);
	tabInFilterField2d1dScaleZ->setObjectName("wavelet.scaleZ");
	tabInFilterField2d1dScaleZ->setMaximumWidth(100);
	tabInFilterField2d1dScaleZ->setMinimum(-1);
	tabInFilterField2d1dScaleZ->setMaximum(50);
	connect(tabInFilterField2d1dScaleZ, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	
	tabInFilterField2d1dPositivity = new QCheckBox(tr("Enable"), tabInFilterGroupBox3);
	tabInFilterField2d1dPositivity->setObjectName("wavelet.positivity");
	tabInFilterField2d1dPositivity->setChecked(false);
	connect(tabInFilterField2d1dPositivity, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabInFilterForm3->addRow(tr("Threshold:"), tabInFilterField2d1dThreshold);
	tabInFilterForm3->addRow(tr("Iterations:"), tabInFilterField2d1dIterations);
	tabInFilterForm3->addRow(tr("Scale XY:"), tabInFilterField2d1dScaleXY);
	tabInFilterForm3->addRow(tr("Scale Z:"), tabInFilterField2d1dScaleZ);
	tabInFilterForm3->addRow(tr("Positivity:"), tabInFilterField2d1dPositivity);
	tabInFilterGroupBox3->setLayout(tabInFilterForm3);
	
	
	tabInFilterButtonPrev = new QPushButton(tr("Previous"), tabInFilter);
	tabInFilterButtonPrev->setIcon(iconGoPreviousView);
	connect(tabInFilterButtonPrev, SIGNAL(clicked()), this, SLOT(displayPrevTab()));
	tabInFilterButtonNext = new QPushButton(tr("Next"), tabInFilter);
	tabInFilterButtonNext->setIcon(iconGoNextView);
	connect(tabInFilterButtonNext, SIGNAL(clicked()), this, SLOT(displayNextTab()));
	tabInFilterLayoutControls = new QHBoxLayout();
	tabInFilterLayoutControls->setContentsMargins(0, 0, 0, 0);
	tabInFilterLayoutControls->setSpacing(0);
	tabInFilterLayoutControls->addWidget(tabInFilterButtonPrev);
	tabInFilterLayoutControls->addStretch();
	tabInFilterLayoutControls->addWidget(tabInFilterButtonNext);
	tabInFilterWidgetControls = new QWidget(tabInFilter);
	tabInFilterWidgetControls->setLayout(tabInFilterLayoutControls);
	
	toolBoxIF->addItem(tabInFilterGroupBox1, iconTaskReject, tr("Smoothing"));
	toolBoxIF->addItem(tabInFilterGroupBox2, iconTaskReject, tr("Noise Scaling"));
	toolBoxIF->addItem(tabInFilterGroupBox3, iconTaskReject, tr("2D-1D Wavelet Filter"));
	
	tabInFilterGroupBox1->setLayout(tabInFilterForm1);
	tabInFilterGroupBox2->setLayout(tabInFilterForm2);
	tabInFilterLayout->addWidget(toolBoxIF);
	tabInFilterLayout->addStretch();
	tabInFilterLayout->addWidget(tabInFilterWidgetControls);
	
	tabInFilter->setLayout(tabInFilterLayout);
	
	
	
	
	
	// Set up source finding tab
	// -------------------------
	
	toolBoxSF = new QToolBox(tabSourceFinding);
	
	tabSourceFindingLayout = new QVBoxLayout();
	
	//S+C
	tabSourceFindingGroupBox1 = new QGroupBox(tr("Enable"), toolBoxSF);
	tabSourceFindingGroupBox1->setObjectName("steps.doSCfind");
	tabSourceFindingGroupBox1->setCheckable(true);
	tabSourceFindingGroupBox1->setChecked(true);
	connect(tabSourceFindingGroupBox1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabSourceFindingGroupBox1, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabSourceFindingWidget1Left = new QWidget(tabSourceFindingGroupBox1);
	tabSourceFindingWidget1Right = new QWidget(tabSourceFindingGroupBox1);
	tabSourceFindingForm1Left = new QFormLayout();
	tabSourceFindingForm1Right = new QFormLayout();
	tabSourceFindingForm1Layout = new QHBoxLayout();
	
	tabSourceFindingFieldThreshold  = new QLineEdit(tabSourceFindingGroupBox1);
	tabSourceFindingFieldThreshold->setObjectName("SCfind.threshold");
	tabSourceFindingFieldThreshold->setMaximumWidth(100);
	tabSourceFindingFieldThreshold->setMaxLength(10);
	connect(tabSourceFindingFieldThreshold, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabSourceFindingFieldThreshold, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldRmsMode = new QComboBox(tabSourceFindingGroupBox1);
	tabSourceFindingFieldRmsMode->setObjectName("SCfind.rmsMode");
	tabSourceFindingFieldRmsMode->addItem(tr("Gaussian fit to negative fluxes"), QVariant(QString("negative")));
	tabSourceFindingFieldRmsMode->addItem(tr("Gaussian fit to flux histogram"), QVariant(QString("gauss")));
	tabSourceFindingFieldRmsMode->addItem(tr("Median absolute deviation"), QVariant(QString("mad")));
	tabSourceFindingFieldRmsMode->addItem(tr("Standard deviation"), QVariant(QString("std")));
	connect(tabSourceFindingFieldRmsMode, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFields()));
	connect(tabSourceFindingFieldRmsMode, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldFluxRange = new QComboBox(tabSourceFindingGroupBox1);
	tabSourceFindingFieldFluxRange->setObjectName("SCfind.fluxRange");
	tabSourceFindingFieldFluxRange->addItem(tr("Negative"), QVariant(QString("negative")));
	tabSourceFindingFieldFluxRange->addItem(tr("Positive"), QVariant(QString("positive")));
	tabSourceFindingFieldFluxRange->addItem(tr("All"), QVariant(QString("all")));
	connect(tabSourceFindingFieldFluxRange, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldEdgeMode = new QComboBox(tabSourceFindingGroupBox1);
	tabSourceFindingFieldEdgeMode->setObjectName("SCfind.edgeMode");
	tabSourceFindingFieldEdgeMode->addItem(tr("Constant"), QVariant(QString("constant")));
	tabSourceFindingFieldEdgeMode->addItem(tr("Reflect"), QVariant(QString("reflect")));
	tabSourceFindingFieldEdgeMode->addItem(tr("Mirror"), QVariant(QString("mirror")));
	tabSourceFindingFieldEdgeMode->addItem(tr("Nearest"), QVariant(QString("nearest")));
	tabSourceFindingFieldEdgeMode->addItem(tr("Wrap"), QVariant(QString("wrap")));
	connect(tabSourceFindingFieldEdgeMode, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldKunit = new QComboBox(tabSourceFindingGroupBox1);
	tabSourceFindingFieldKunit->setObjectName("SCfind.kernelUnit");
	tabSourceFindingFieldKunit->addItem(tr("Pixels"), QVariant(QString("pixel")));
	tabSourceFindingFieldKunit->addItem(tr("World coordinates"), QVariant(QString("world")));
	connect(tabSourceFindingFieldKunit, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldKernels = new QTextEdit(tabSourceFindingGroupBox1);
	tabSourceFindingFieldKernels->setObjectName("SCfind.kernels");
	tabSourceFindingFieldKernels->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	tabSourceFindingFieldKernels->setMaximumHeight(120);
	connect(tabSourceFindingFieldKernels, SIGNAL(textChanged()), this, SLOT(parameterChanged()));
	
	tabSourceFindingForm1Left->addRow(tr("Threshold:"), tabSourceFindingFieldThreshold);
	tabSourceFindingForm1Left->addRow(tr("RMS mode:"), tabSourceFindingFieldRmsMode);
	tabSourceFindingForm1Left->addRow(tr("Flux range:"), tabSourceFindingFieldFluxRange);
	tabSourceFindingForm1Left->addRow(tr("Edge mode:"), tabSourceFindingFieldEdgeMode);
	tabSourceFindingForm1Left->addRow(tr("Kernel units:"), tabSourceFindingFieldKunit);
	
	tabSourceFindingForm1Right->addRow(tr("Kernels:"), tabSourceFindingFieldKernels);
	tabSourceFindingForm1Right->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
	
	tabSourceFindingWidget1Left->setLayout(tabSourceFindingForm1Left);
	tabSourceFindingWidget1Right->setLayout(tabSourceFindingForm1Right);
	
	tabSourceFindingForm1Layout->addWidget(tabSourceFindingWidget1Left);
	tabSourceFindingForm1Layout->addWidget(tabSourceFindingWidget1Right);
	tabSourceFindingForm1Layout->setContentsMargins(0, 0, 0, 0);
	tabSourceFindingForm1Layout->setStretchFactor(tabSourceFindingWidget1Left, 0);
	tabSourceFindingForm1Layout->setStretchFactor(tabSourceFindingWidget1Right, 1);
	
	tabSourceFindingGroupBox1->setLayout(tabSourceFindingForm1Layout);
	
	
	
	//CNHI
	tabSourceFindingGroupBox3 = new QGroupBox(tr("Enable"), toolBoxSF);
	tabSourceFindingGroupBox3->setObjectName("steps.doCNHI");
	tabSourceFindingGroupBox3->setCheckable(true);
	tabSourceFindingGroupBox3->setChecked(false);
	connect(tabSourceFindingGroupBox3, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabSourceFindingGroupBox3, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabSourceFindingForm3 = new QFormLayout();
	
	tabSourceFindingFieldPReq  = new QLineEdit(tabSourceFindingGroupBox3);
	tabSourceFindingFieldPReq->setObjectName("CNHI.pReq");
	tabSourceFindingFieldPReq->setMaximumWidth(100);
	tabSourceFindingFieldPReq->setMaxLength(10);
	connect(tabSourceFindingFieldPReq, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabSourceFindingFieldPReq, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldQReq  = new QLineEdit(tabSourceFindingGroupBox3);
	tabSourceFindingFieldQReq->setObjectName("CNHI.qReq");
	tabSourceFindingFieldQReq->setMaximumWidth(100);
	tabSourceFindingFieldQReq->setMaxLength(10);
	connect(tabSourceFindingFieldQReq, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabSourceFindingFieldQReq, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldMinScale = new QSpinBox(tabSourceFindingGroupBox3);
	tabSourceFindingFieldMinScale->setObjectName("CNHI.minScale");
	tabSourceFindingFieldMinScale->setMaximumWidth(100);
	tabSourceFindingFieldMinScale->setMinimum(1);
	tabSourceFindingFieldMinScale->setMaximum(50);
	connect(tabSourceFindingFieldMinScale, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldMaxScale = new QSpinBox(tabSourceFindingGroupBox3);
	tabSourceFindingFieldMaxScale->setObjectName("CNHI.maxScale");
	tabSourceFindingFieldMaxScale->setMaximumWidth(100);
	tabSourceFindingFieldMaxScale->setMinimum(-1);
	tabSourceFindingFieldMaxScale->setMaximum(50);
	connect(tabSourceFindingFieldMaxScale, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingMedianTest = new QCheckBox(tr("Enable"), tabSourceFindingGroupBox3);
	tabSourceFindingMedianTest->setObjectName("CNHI.medianTest");
	connect(tabSourceFindingMedianTest, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldVerbose = new QComboBox(tabSourceFindingGroupBox3);
	tabSourceFindingFieldVerbose->setObjectName("CNHI.verbose");
	tabSourceFindingFieldVerbose->addItem(tr("None"), QVariant(QString("0")));
	tabSourceFindingFieldVerbose->addItem(tr("Minimal"), QVariant(QString("1")));
	tabSourceFindingFieldVerbose->addItem(tr("Maximal"), QVariant(QString("2")));
	connect(tabSourceFindingFieldVerbose, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingForm3->addRow(tr("Probability:"), tabSourceFindingFieldPReq);
	tabSourceFindingForm3->addRow(tr("Quality:"), tabSourceFindingFieldQReq);
	tabSourceFindingForm3->addRow(tr("Min. scale:"), tabSourceFindingFieldMinScale);
	tabSourceFindingForm3->addRow(tr("Max. scale:"), tabSourceFindingFieldMaxScale);
	tabSourceFindingForm3->addRow(tr("Median test:"), tabSourceFindingMedianTest);
	tabSourceFindingForm3->addRow(tr("Verbosity:"), tabSourceFindingFieldVerbose);
	
	tabSourceFindingGroupBox3->setLayout(tabSourceFindingForm3);
	
	
	
	//Threshold
	tabSourceFindingGroupBox2 = new QGroupBox(tr("Enable"), toolBoxSF);
	tabSourceFindingGroupBox2->setObjectName("steps.doThreshold");
	tabSourceFindingGroupBox2->setCheckable(true);
	tabSourceFindingGroupBox2->setChecked(false);
	connect(tabSourceFindingGroupBox2, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabSourceFindingGroupBox2, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabSourceFindingForm2 = new QFormLayout();
	
	tabSourceFindingFieldThreshold2 = new QLineEdit(tabSourceFindingGroupBox2);
	tabSourceFindingFieldThreshold2->setObjectName("threshold.threshold");
	tabSourceFindingFieldThreshold2->setMaximumWidth(100);
	tabSourceFindingFieldThreshold2->setMaxLength(10);
	connect(tabSourceFindingFieldThreshold2, SIGNAL(editingFinished()), this, SLOT(updateFields()));
	connect(tabSourceFindingFieldThreshold2, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldClipMethod = new QComboBox(tabSourceFindingGroupBox2);
	tabSourceFindingFieldClipMethod->setObjectName("threshold.clipMethod");
	tabSourceFindingFieldClipMethod->addItem(tr("Relative"), QVariant(QString("relative")));
	tabSourceFindingFieldClipMethod->addItem(tr("Absolute"), QVariant(QString("absolute")));
	connect(tabSourceFindingFieldClipMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFields()));
	connect(tabSourceFindingFieldClipMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldRmsMode2 = new QComboBox(tabSourceFindingGroupBox2);
	tabSourceFindingFieldRmsMode2->setObjectName("threshold.rmsMode");
	tabSourceFindingFieldRmsMode2->addItem(tr("Gaussian fit to negative fluxes"), QVariant(QString("negative")));
	tabSourceFindingFieldRmsMode2->addItem(tr("Gaussian fit to flux histogram"), QVariant(QString("gauss")));
	tabSourceFindingFieldRmsMode2->addItem(tr("Median absolute deviation"), QVariant(QString("mad")));
	tabSourceFindingFieldRmsMode2->addItem(tr("Standard deviation"), QVariant(QString("std")));
	connect(tabSourceFindingFieldRmsMode2, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFields()));
	connect(tabSourceFindingFieldRmsMode2, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingFieldFluxRange2 = new QComboBox(tabSourceFindingGroupBox2);
	tabSourceFindingFieldFluxRange2->setObjectName("threshold.fluxRange");
	tabSourceFindingFieldFluxRange2->addItem(tr("Negative"), QVariant(QString("negative")));
	tabSourceFindingFieldFluxRange2->addItem(tr("Positive"), QVariant(QString("positive")));
	tabSourceFindingFieldFluxRange2->addItem(tr("All"), QVariant(QString("all")));
	connect(tabSourceFindingFieldFluxRange2, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged()));
	
	tabSourceFindingForm2->addRow(tr("Threshold:"), tabSourceFindingFieldThreshold2);
	tabSourceFindingForm2->addRow(tr("Clip mode:"), tabSourceFindingFieldClipMethod);
	tabSourceFindingForm2->addRow(tr("RMS mode:"), tabSourceFindingFieldRmsMode2);
	tabSourceFindingForm2->addRow(tr("Flux range:"), tabSourceFindingFieldFluxRange2);
	tabSourceFindingGroupBox2->setLayout(tabSourceFindingForm2);
	
	
	
	tabSourceFindingButtonPrev = new QPushButton(tr("Previous"), tabSourceFinding);
	tabSourceFindingButtonPrev->setIcon(iconGoPreviousView);
	connect(tabSourceFindingButtonPrev, SIGNAL(clicked()), this, SLOT(displayPrevTab()));
	tabSourceFindingButtonNext = new QPushButton(tr("Next"), tabSourceFinding);
	tabSourceFindingButtonNext->setIcon(iconGoNextView);
	connect(tabSourceFindingButtonNext, SIGNAL(clicked()), this, SLOT(displayNextTab()));
	tabSourceFindingLayoutControls = new QHBoxLayout();
	tabSourceFindingLayoutControls->setContentsMargins(0, 0, 0, 0);
	tabSourceFindingLayoutControls->setSpacing(0);
	tabSourceFindingLayoutControls->addWidget(tabSourceFindingButtonPrev);
	tabSourceFindingLayoutControls->addStretch();
	tabSourceFindingLayoutControls->addWidget(tabSourceFindingButtonNext);
	tabSourceFindingWidgetControls = new QWidget(tabSourceFinding);
	tabSourceFindingWidgetControls->setLayout(tabSourceFindingLayoutControls);
	
	toolBoxSF->addItem(tabSourceFindingGroupBox1, iconTaskReject, tr("Smooth + Clip Finder"));
	toolBoxSF->addItem(tabSourceFindingGroupBox3, iconTaskReject, tr("CNHI Finder"));
	toolBoxSF->addItem(tabSourceFindingGroupBox2, iconTaskReject, tr("Threshold Finder"));
	
	tabSourceFindingLayout->addWidget(toolBoxSF);
	tabSourceFindingLayout->addStretch();
	tabSourceFindingLayout->addWidget(tabSourceFindingWidgetControls);
	tabSourceFinding->setLayout(tabSourceFindingLayout);
	
	
	
	// Set up merging tab
	// ------------------
	
	toolBoxME = new QToolBox(tabMerging);
	
	tabMergingLayout = new QVBoxLayout();
	
	tabMergingGroupBox1 = new QGroupBox(tr("Enable"), toolBoxME);
	tabMergingGroupBox1->setObjectName("steps.doMerge");
	tabMergingGroupBox1->setCheckable(true);
	tabMergingGroupBox1->setChecked(true);
	connect(tabMergingGroupBox1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabMergingGroupBox1, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabMergingForm1Layout = new QHBoxLayout();
	tabMergingForm1Left = new QFormLayout();
	tabMergingForm1Right = new QFormLayout();
	tabMergingWidget1Left = new QWidget(tabMergingGroupBox1);
	tabMergingWidget1Right = new QWidget(tabMergingGroupBox1);
	
	tabMergingFieldRadiusX = new QSpinBox(tabMergingWidget1Left);
	tabMergingFieldRadiusX->setObjectName("merge.radiusX");
	tabMergingFieldRadiusX->setMaximumWidth(100);
	tabMergingFieldRadiusX->setMinimum(0);
	tabMergingFieldRadiusX->setMaximum(50);
	connect(tabMergingFieldRadiusX, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabMergingFieldRadiusY = new QSpinBox(tabMergingWidget1Left);
	tabMergingFieldRadiusY->setObjectName("merge.radiusY");
	tabMergingFieldRadiusY->setMaximumWidth(100);
	tabMergingFieldRadiusY->setMinimum(0);
	tabMergingFieldRadiusY->setMaximum(50);
	connect(tabMergingFieldRadiusY, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabMergingFieldRadiusZ = new QSpinBox(tabMergingWidget1Left);
	tabMergingFieldRadiusZ->setObjectName("merge.radiusZ");
	tabMergingFieldRadiusZ->setMaximumWidth(100);
	tabMergingFieldRadiusZ->setMinimum(0);
	tabMergingFieldRadiusZ->setMaximum(50);
	connect(tabMergingFieldRadiusZ, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabMergingFieldMinSizeX = new QSpinBox(tabMergingWidget1Right);
	tabMergingFieldMinSizeX->setObjectName("merge.minSizeX");
	tabMergingFieldMinSizeX->setMaximumWidth(100);
	tabMergingFieldMinSizeX->setMinimum(1);
	tabMergingFieldMinSizeX->setMaximum(50);
	connect(tabMergingFieldMinSizeX, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabMergingFieldMinSizeY = new QSpinBox(tabMergingWidget1Right);
	tabMergingFieldMinSizeY->setObjectName("merge.minSizeY");
	tabMergingFieldMinSizeY->setMaximumWidth(100);
	tabMergingFieldMinSizeY->setMinimum(1);
	tabMergingFieldMinSizeY->setMaximum(50);
	connect(tabMergingFieldMinSizeY, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabMergingFieldMinSizeZ = new QSpinBox(tabMergingWidget1Right);
	tabMergingFieldMinSizeZ->setObjectName("merge.minSizeZ");
	tabMergingFieldMinSizeZ->setMaximumWidth(100);
	tabMergingFieldMinSizeZ->setMinimum(1);
	tabMergingFieldMinSizeZ->setMaximum(50);
	connect(tabMergingFieldMinSizeZ, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	
	tabMergingForm1Left->addRow(tr("Radius X:"), tabMergingFieldRadiusX);
	tabMergingForm1Left->addRow(tr("Radius Y:"), tabMergingFieldRadiusY);
	tabMergingForm1Left->addRow(tr("Radius Z:"), tabMergingFieldRadiusZ);
	tabMergingForm1Left->setContentsMargins(0, 0, 0, 0);
	tabMergingForm1Right->addRow(tr("Min. size X:"), tabMergingFieldMinSizeX);
	tabMergingForm1Right->addRow(tr("Min. size Y:"), tabMergingFieldMinSizeY);
	tabMergingForm1Right->addRow(tr("Min. size Z:"), tabMergingFieldMinSizeZ);
	tabMergingForm1Right->setContentsMargins(0, 0, 0, 0);
	tabMergingWidget1Left->setLayout(tabMergingForm1Left);
	tabMergingWidget1Right->setLayout(tabMergingForm1Right);
	
	tabMergingForm1Layout->addWidget(tabMergingWidget1Left);
	tabMergingForm1Layout->addWidget(tabMergingWidget1Right);
	tabMergingForm1Layout->addStretch();
	tabMergingForm1Layout->setSpacing(15);
	tabMergingGroupBox1->setLayout(tabMergingForm1Layout);
	
	
	tabMergingGroupBox2 = new QGroupBox(toolBoxME);
	tabMergingForm2 = new QFormLayout();
	
	tabMergingButtonPositivity = new QCheckBox(tr("Enable"), tabMergingGroupBox2);
	tabMergingButtonPositivity->setObjectName("merge.positivity");
	tabMergingButtonPositivity->setEnabled(true);
	tabMergingButtonPositivity->setChecked(false);
	connect(tabMergingButtonPositivity, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	connect(tabMergingButtonPositivity, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabMergingButtonPositivity, SIGNAL(clicked(bool)), this, SLOT(printPositivityWarning(bool)));
	
	tabMergingLabelWarning = new QLabel(tabMergingGroupBox2);
	tabMergingLabelWarning->setText(tr("<p><strong>Warning:</strong> Enabling positivity is dangerous and will render some of SoFiA&rsquo;s most powerful algorithms useless, including mask optimisation and reliability calculation. It can also create biases in certain source parameters measured by SoFiA.</p><p>This option should only be activated by expert users who are fully aware of the risks and consequences of discarding negative signals.</p>"));
	tabMergingLabelWarning->setWordWrap(true);
	
	tabMergingForm2->addRow(tr("Positivity:"), tabMergingButtonPositivity);
	tabMergingForm2->addRow(tabMergingLabelWarning);
	tabMergingGroupBox2->setLayout(tabMergingForm2);
	
	
	tabMergingButtonPrev = new QPushButton(tr("Previous"), tabMerging);
	tabMergingButtonPrev->setIcon(iconGoPreviousView);
	connect(tabMergingButtonPrev, SIGNAL(clicked()), this, SLOT(displayPrevTab()));
	tabMergingButtonNext = new QPushButton(tr("Next"), tabMerging);
	tabMergingButtonNext->setIcon(iconGoNextView);
	connect(tabMergingButtonNext, SIGNAL(clicked()), this, SLOT(displayNextTab()));
	tabMergingLayoutControls = new QHBoxLayout();
	tabMergingLayoutControls->setContentsMargins(0, 0, 0, 0);
	tabMergingLayoutControls->setSpacing(0);
	tabMergingLayoutControls->addWidget(tabMergingButtonPrev);
	tabMergingLayoutControls->addStretch();
	tabMergingLayoutControls->addWidget(tabMergingButtonNext);
	tabMergingWidgetControls = new QWidget(tabMerging);
	tabMergingWidgetControls->setLayout(tabMergingLayoutControls);
	
	toolBoxME->addItem(tabMergingGroupBox1, iconTaskReject, tr("Merging of Detections"));
	toolBoxME->addItem(tabMergingGroupBox2, iconTaskReject, tr("Additional Options"));
	
	tabMergingLayout->addWidget(toolBoxME);
	tabMergingLayout->addStretch();
	tabMergingLayout->addWidget(tabMergingWidgetControls);
	tabMerging->setLayout(tabMergingLayout);
	
	
	
	
	// Set up parametrisation tab
	// --------------------------
	
	toolBoxPA = new QToolBox(tabParametrisation);
	
	tabParametrisationLayout = new QVBoxLayout();
	
	tabParametrisationGroupBox1 = new QGroupBox(tr("Enable"), toolBoxPA);
	tabParametrisationGroupBox1->setObjectName("steps.doParameterise");
	tabParametrisationGroupBox1->setCheckable(true);
	tabParametrisationGroupBox1->setChecked(true);
	connect(tabParametrisationGroupBox1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabParametrisationGroupBox1, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabParametrisationForm1 = new QFormLayout();
	
	tabParametrisationButtonDilateMask = new QCheckBox(tr("Optimise mask (dilation)"), tabParametrisationGroupBox1);
	tabParametrisationButtonDilateMask->setObjectName("parameters.dilateMask");
	tabParametrisationButtonDilateMask->setEnabled(true);
	tabParametrisationButtonDilateMask->setChecked(false);
	connect(tabParametrisationButtonDilateMask, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabParametrisationButtonMaskOpt = new QCheckBox(tr("Optimise mask (ellipse)"), tabParametrisationGroupBox1);
	tabParametrisationButtonMaskOpt->setObjectName("parameters.optimiseMask");
	tabParametrisationButtonMaskOpt->setEnabled(true);
	tabParametrisationButtonMaskOpt->setChecked(false);
	connect(tabParametrisationButtonMaskOpt, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabParametrisationButtonBusyFunction = new QCheckBox(tr("Fit Busy Function"), tabParametrisationGroupBox1);
	tabParametrisationButtonBusyFunction->setObjectName("parameters.fitBusyFunction");
	tabParametrisationButtonBusyFunction->setEnabled(true);
	tabParametrisationButtonBusyFunction->setChecked(false);
	connect(tabParametrisationButtonBusyFunction, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabParametrisationForm1->addRow(tr(""), tabParametrisationButtonDilateMask);
	tabParametrisationForm1->addRow(tr(""), tabParametrisationButtonMaskOpt);
	tabParametrisationForm1->addRow(tr(""), tabParametrisationButtonBusyFunction);
	tabParametrisationGroupBox1->setLayout(tabParametrisationForm1);
	
	tabParametrisationGroupBox2 = new QGroupBox(tr("Enable"), toolBoxPA);
	tabParametrisationGroupBox2->setObjectName("steps.doReliability");
	tabParametrisationGroupBox2->setCheckable(true);
	tabParametrisationGroupBox2->setChecked(true);
	connect(tabParametrisationGroupBox2, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabParametrisationGroupBox2, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabParametrisationForm2 = new QFormLayout();
	
	tabParametrisationWidgetRelMin = new QWidget(tabParametrisationGroupBox1);
	tabParametrisationLayoutRelMin = new QHBoxLayout();
	tabParametrisationSliderRelMin = new QSlider(Qt::Horizontal, tabParametrisationWidgetRelMin);
	tabParametrisationSliderRelMin->setMinimum(0);
	tabParametrisationSliderRelMin->setMaximum(100);
	tabParametrisationSliderRelMin->setValue(0.1 * RELMIN_SCALE_FACTOR);
	tabParametrisationSliderRelMin->setSingleStep(1);
	tabParametrisationSliderRelMin->setPageStep(5);
	tabParametrisationSliderRelMin->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
	connect(tabParametrisationSliderRelMin, SIGNAL(valueChanged(int)), this, SLOT(updateFields()));
	connect(tabParametrisationSliderRelMin, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabParametrisationFieldRelMin = new QLineEdit(tabParametrisationGroupBox2);
	tabParametrisationFieldRelMin->setObjectName("reliability.threshold");
	tabParametrisationFieldRelMin->setMaxLength(10);
	tabParametrisationFieldRelMin->setReadOnly(true);
	tabParametrisationFieldRelMin->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	tabParametrisationLayoutRelMin->addWidget(tabParametrisationSliderRelMin);
	tabParametrisationLayoutRelMin->addWidget(tabParametrisationFieldRelMin);
	tabParametrisationLayoutRelMin->setContentsMargins(0, 0, 0, 0);
	tabParametrisationWidgetRelMin->setLayout(tabParametrisationLayoutRelMin);
	
	tabParametrisationWidgetScaleKernel = new QWidget(tabParametrisationGroupBox2);
	tabParametrisationLayoutScaleKernel = new QHBoxLayout();
	tabParametrisationSliderScaleKernel = new QSlider(Qt::Horizontal, tabParametrisationWidgetScaleKernel);
	tabParametrisationSliderScaleKernel->setMinimum(0);
	tabParametrisationSliderScaleKernel->setMaximum(100);
	tabParametrisationSliderScaleKernel->setSingleStep(1);
	tabParametrisationSliderScaleKernel->setPageStep(5);
	tabParametrisationSliderScaleKernel->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
	connect(tabParametrisationSliderScaleKernel, SIGNAL(valueChanged(int)), this, SLOT(updateFields()));
	connect(tabParametrisationSliderScaleKernel, SIGNAL(valueChanged(int)), this, SLOT(parameterChanged()));
	tabParametrisationFieldScaleKernel = new QLineEdit(tabParametrisationWidgetScaleKernel);
	tabParametrisationFieldScaleKernel->setObjectName("reliability.scaleKernel");
	tabParametrisationFieldScaleKernel->setMaxLength(10);
	tabParametrisationFieldScaleKernel->setReadOnly(true);
	tabParametrisationFieldScaleKernel->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	tabParametrisationLayoutScaleKernel->addWidget(tabParametrisationSliderScaleKernel);
	tabParametrisationLayoutScaleKernel->addWidget(tabParametrisationFieldScaleKernel);
	tabParametrisationLayoutScaleKernel->setContentsMargins(0, 0, 0, 0);
	tabParametrisationWidgetScaleKernel->setLayout(tabParametrisationLayoutScaleKernel);
	
	
	tabParametrisationButtonRelPlot = new QCheckBox(tr("Enable"), tabParametrisationGroupBox2);
	tabParametrisationButtonRelPlot->setObjectName("reliability.makePlot");
	tabParametrisationButtonRelPlot->setEnabled(true);
	tabParametrisationButtonRelPlot->setChecked(false);
	connect(tabParametrisationButtonRelPlot, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabParametrisationForm2->addRow(tr("Threshold:"), tabParametrisationWidgetRelMin);
	tabParametrisationForm2->addRow(tr("Kernel scale:"), tabParametrisationWidgetScaleKernel);
	tabParametrisationForm2->addRow(tr("Diagnostic plot:"), tabParametrisationButtonRelPlot);
	tabParametrisationForm2->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
	tabParametrisationGroupBox2->setLayout(tabParametrisationForm2);
	
	tabParametrisationButtonPrev = new QPushButton(tr("Previous"), tabParametrisation);
	tabParametrisationButtonPrev->setIcon(iconGoPreviousView);
	connect(tabParametrisationButtonPrev, SIGNAL(clicked()), this, SLOT(displayPrevTab()));
	tabParametrisationButtonNext = new QPushButton(tr("Next"), tabParametrisation);
	tabParametrisationButtonNext->setIcon(iconGoNextView);
	connect(tabParametrisationButtonNext, SIGNAL(clicked()), this, SLOT(displayNextTab()));
	tabParametrisationLayoutControls = new QHBoxLayout();
	tabParametrisationLayoutControls->setContentsMargins(0, 0, 0, 0);
	tabParametrisationLayoutControls->setSpacing(0);
	tabParametrisationLayoutControls->addWidget(tabParametrisationButtonPrev);
	tabParametrisationLayoutControls->addStretch();
	tabParametrisationLayoutControls->addWidget(tabParametrisationButtonNext);
	tabParametrisationWidgetControls = new QWidget(tabParametrisation);
	tabParametrisationWidgetControls->setLayout(tabParametrisationLayoutControls);
	
	toolBoxPA->addItem(tabParametrisationGroupBox1, iconTaskReject, tr("Source Parameterisation"));
	toolBoxPA->addItem(tabParametrisationGroupBox2, iconTaskReject, tr("Reliability Calculation"));
	
	tabParametrisationLayout->addWidget(toolBoxPA);
	tabParametrisationLayout->addStretch();
	tabParametrisationLayout->addWidget(tabParametrisationWidgetControls);
	tabParametrisation->setLayout(tabParametrisationLayout);
	
	
	
	// Set up output tab
	// -----------------
	
	toolBoxOP = new QToolBox(tabOutput);
	
	tabOutputLayout = new QVBoxLayout();
	
	tabOutputGroupBox1 = new QGroupBox(toolBoxOP);
		
	tabOutputButtonFilteredCube = new QCheckBox(tr("Filtered cube"), tabOutputGroupBox1);
	tabOutputButtonFilteredCube->setObjectName("steps.doWriteFilteredCube");
	tabOutputButtonFilteredCube->setChecked(false);
	connect(tabOutputButtonFilteredCube, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonFilteredCube, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabOutputButtonMask = new QCheckBox(tr("Mask"), tabOutputGroupBox1);
	tabOutputButtonMask->setObjectName("steps.doWriteMask");
	tabOutputButtonMask->setChecked(false);
	connect(tabOutputButtonMask, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonMask, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabOutputButtonMom0 = new QCheckBox(tr("Mom. 0"), tabOutputGroupBox1);
	tabOutputButtonMom0->setObjectName("steps.doMom0");
	tabOutputButtonMom0->setChecked(false);
	connect(tabOutputButtonMom0, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonMom0, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabOutputButtonMom1 = new QCheckBox(tr("Mom. 1"), tabOutputGroupBox1);
	tabOutputButtonMom1->setObjectName("steps.doMom1");
	tabOutputButtonMom1->setChecked(false);
	connect(tabOutputButtonMom1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonMom1, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabOutputButtonNoiseCube = new QCheckBox(tr("Noise cube"), tabOutputGroupBox1);
	tabOutputButtonNoiseCube->setObjectName("steps.doWriteNoiseCube");
	tabOutputButtonNoiseCube->setChecked(false);
	connect(tabOutputButtonNoiseCube, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonNoiseCube, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabOutputButtonCubelets = new QCheckBox(tr("Source products"), tabOutputGroupBox1);
	tabOutputButtonCubelets->setObjectName("steps.doCubelets");
	tabOutputButtonCubelets->setChecked(false);
	connect(tabOutputButtonCubelets, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonCubelets, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabOutputButtonASCII = new QCheckBox(tr("ASCII"), tabOutputGroupBox1);
	tabOutputButtonASCII->setObjectName("writeCat.writeASCII");
	tabOutputButtonASCII->setChecked(true);
	tabOutputButtonASCII->setEnabled(true);
	connect(tabOutputButtonASCII, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonASCII, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabOutputButtonXML = new QCheckBox(tr("VO table"), tabOutputGroupBox1);
	tabOutputButtonXML->setObjectName("writeCat.writeXML");
	tabOutputButtonXML->setChecked(false);
	tabOutputButtonXML->setEnabled(true);
	connect(tabOutputButtonXML, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonXML, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	tabOutputButtonSQL = new QCheckBox(tr("SQL"), tabOutputGroupBox1);
	tabOutputButtonSQL->setObjectName("writeCat.writeSQL");
	tabOutputButtonSQL->setChecked(false);
	tabOutputButtonSQL->setEnabled(true);
	connect(tabOutputButtonSQL, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonSQL, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabOutputFieldParameters = new QLineEdit(tabOutputGroupBox1);
	tabOutputFieldParameters->setObjectName("writeCat.parameters");
	tabOutputFieldParameters->setEnabled(true);
	connect(tabOutputFieldParameters, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabOutputLabelParameters = new QLabel(tr("<b>Note:</b> Depending on the actual pipeline settings, some selected parameters may not be written. Please see the SoFiA user manual or <a href=\"https://github.com/SoFiA-Admin/SoFiA/wiki/SoFiA-Source-Parameters\">wiki</a> for a complete list of source parameters."), tabOutputGroupBox1);
	tabOutputLabelParameters->setWordWrap(true);
	tabOutputLabelParameters->setContentsMargins(0, 0, 0, 0);
	tabOutputLabelParameters->setAlignment(Qt::AlignTop | Qt::AlignLeft);
	tabOutputLabelParameters->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Preferred);
	
	tabOutputWidgetProducts = new QWidget(tabOutputGroupBox1);
	tabOutputLayoutProducts = new QHBoxLayout();
	tabOutputLayoutProducts->setContentsMargins(0, 0, 0, 0);
	tabOutputLayoutProducts->setSpacing(10);
	tabOutputLayoutProducts->addWidget(tabOutputButtonFilteredCube);
	tabOutputLayoutProducts->addWidget(tabOutputButtonMask);
	tabOutputLayoutProducts->addWidget(tabOutputButtonMom0);
	tabOutputLayoutProducts->addWidget(tabOutputButtonMom1);
	tabOutputLayoutProducts->addStretch();
	tabOutputWidgetProducts->setLayout(tabOutputLayoutProducts);
	
	tabOutputWidgetProducts2 = new QWidget(tabOutputGroupBox1);
	tabOutputLayoutProducts2 = new QHBoxLayout();
	tabOutputLayoutProducts2->setContentsMargins(0, 0, 0, 0);
	tabOutputLayoutProducts2->setSpacing(10);
	tabOutputLayoutProducts2->addWidget(tabOutputButtonNoiseCube);
	tabOutputLayoutProducts2->addWidget(tabOutputButtonCubelets);
	tabOutputLayoutProducts2->addStretch();
	tabOutputWidgetProducts2->setLayout(tabOutputLayoutProducts2);
	
	tabOutputWidgetFormat = new QWidget(tabOutputGroupBox1);
	tabOutputLayoutFormat = new QHBoxLayout();
	tabOutputLayoutFormat->setContentsMargins(0, 0, 0, 0);
	tabOutputLayoutFormat->setSpacing(10);
	tabOutputLayoutFormat->addWidget(tabOutputButtonASCII);
	tabOutputLayoutFormat->addWidget(tabOutputButtonXML);
	tabOutputLayoutFormat->addWidget(tabOutputButtonSQL);
	tabOutputLayoutFormat->addStretch();
	tabOutputWidgetFormat->setLayout(tabOutputLayoutFormat);
	
	tabOutputForm1 = new QFormLayout();
	tabOutputForm1->addRow(tr("Data products:"), tabOutputWidgetProducts);
	tabOutputForm1->addRow(tr(""), tabOutputWidgetProducts2);
	tabOutputForm1->addRow(tr("Source catalogue:"), tabOutputWidgetFormat);
	tabOutputForm1->addRow(tr("Parameters:"), tabOutputFieldParameters);
	tabOutputForm1->addRow(tr(""), tabOutputLabelParameters);
	tabOutputForm1->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
	tabOutputGroupBox1->setLayout(tabOutputForm1);
	
	
	tabOutputGroupBox2 = new QGroupBox(toolBoxOP);
	
	tabOutputFieldBaseName = new QLineEdit(tabOutputGroupBox2);
	tabOutputFieldBaseName->setObjectName("writeCat.basename");
	tabOutputFieldBaseName->setEnabled(true);
	connect(tabOutputFieldBaseName, SIGNAL(textChanged(const QString &)), this, SLOT(updateFields()));
	connect(tabOutputFieldBaseName, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	
	tabOutputWidgetDirectory = new QWidget(tabOutputGroupBox2);
	tabOutputLayoutDirectory = new QHBoxLayout();
	tabOutputFieldDirectory  = new QLineEdit(tabOutputWidgetDirectory);
	tabOutputFieldDirectory->setObjectName("writeCat.outputDir");
	connect(tabOutputFieldDirectory, SIGNAL(textChanged(const QString &)), this, SLOT(updateFields()));
	connect(tabOutputFieldDirectory, SIGNAL(textChanged(const QString &)), this, SLOT(parameterChanged()));
	tabOutputButtonDirectory = new QPushButton(tr("Select..."), tabOutputWidgetDirectory);
	connect(tabOutputButtonDirectory, SIGNAL(clicked()), this, SLOT(selectOutputDirectory()));
	tabOutputButtonDirectory->setIcon(iconDocumentOpen);
	tabOutputLayoutDirectory->addWidget(tabOutputFieldDirectory);
	tabOutputLayoutDirectory->addWidget(tabOutputButtonDirectory);
	tabOutputLayoutDirectory->setContentsMargins(0, 0, 0, 0);
	tabOutputWidgetDirectory->setLayout(tabOutputLayoutDirectory);
	
	tabOutputButtonCompress = new QCheckBox(tr("Enable"), tabOutputGroupBox2);
	tabOutputButtonCompress->setObjectName("writeCat.compress");
	tabOutputButtonCompress->setChecked(false);
	connect(tabOutputButtonCompress, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonCompress, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabOutputButtonOverwrite = new QCheckBox(tr("Enable"), tabOutputGroupBox2);
	tabOutputButtonOverwrite->setObjectName("writeCat.overwrite");
	tabOutputButtonOverwrite->setChecked(false);
	connect(tabOutputButtonOverwrite, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
	connect(tabOutputButtonOverwrite, SIGNAL(toggled(bool)), this, SLOT(parameterChanged()));
	
	tabOutputForm2 = new QFormLayout();
	tabOutputForm2->addRow(tr("Base name:"), tabOutputFieldBaseName);
	tabOutputForm2->addRow(tr("Output directory:"), tabOutputWidgetDirectory);
	tabOutputForm2->addRow(tr("Compression:"), tabOutputButtonCompress);
	tabOutputForm2->addRow(tr("Overwrite files:"), tabOutputButtonOverwrite);
	tabOutputForm2->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
	tabOutputGroupBox2->setLayout(tabOutputForm2);
	
	tabOutputButtonPrev = new QPushButton(tr("Previous"), tabOutput);
	tabOutputButtonPrev->setIcon(iconGoPreviousView);
	connect(tabOutputButtonPrev, SIGNAL(clicked()), this, SLOT(displayPrevTab()));
	tabOutputButtonGo   = new QPushButton(tr("Run Pipeline"), tabOutput);
	tabOutputButtonGo->setIcon(iconDialogOkApply);
	connect(tabOutputButtonGo, SIGNAL(clicked()), this, SLOT(runPipeline()));
	tabOutputLayoutControls = new QHBoxLayout();
	tabOutputLayoutControls->setContentsMargins(0, 0, 0, 0);
	tabOutputLayoutControls->setSpacing(0);
	tabOutputLayoutControls->addWidget(tabOutputButtonPrev);
	tabOutputLayoutControls->addStretch();
	tabOutputLayoutControls->addWidget(tabOutputButtonGo);
	tabOutputWidgetControls = new QWidget(tabOutput);
	tabOutputWidgetControls->setLayout(tabOutputLayoutControls);
	
	toolBoxOP->addItem(tabOutputGroupBox1, iconTaskReject, tr("Output Data Products"));
	toolBoxOP->addItem(tabOutputGroupBox2, iconTaskReject, tr("Additional Options"));
	
	tabOutputLayout->addWidget(toolBoxOP);
	tabOutputLayout->addStretch();
	tabOutputLayout->addWidget(tabOutputWidgetControls);
	tabOutput->setLayout(tabOutputLayout);
	
	// Generate What's This? entries for all widgets
	// ---------------------------------------------
	
	createWhatsThis();
	
	// Set up output widget
	// --------------------
	
	widgetOutput = new QWidget(widgetMain);
	
	outputText = new QTextEdit(widgetOutput);
	outputText->setWhatsThis(tr("<h3>Pipeline Messages</h3><p>Any output messages and progress information produced by the SoFiA pipeline will be displayed here.</p>"));
	outputText->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
	outputText->setReadOnly(true);
	outputText->setLineWrapMode(QTextEdit::FixedColumnWidth);
	outputText->setLineWrapColumnOrWidth(80);
	outputText->setTabStopWidth(8 * outputText->fontMetrics().width("0"));     // setTabStopWidth() expects pixels!!!
	
	QFont font = QFont("Courier");
	font.setStyleHint(QFont::TypeWriter, QFont::PreferAntialias);
	font.setPointSize(10);
	font.setFixedPitch(true);
	font.setKerning(false);
	outputText->setFont(font);
	
	QFontMetrics fontMetrics(font);
	outputText->setMinimumSize(85 * fontMetrics.width(QChar('0')), 10 * fontMetrics.height());
	
	outputProgress = new QProgressBar(widgetOutput);
	outputProgress->setWhatsThis(tr("<h3>Progress Bar</h3><p>Indicates activity and success of the SoFiA pipeline. Will be moving while the pipeline is running and display 100% in the case of success or 0% in the case of failure once the pipeline run terminates.</p>"));
	outputProgress->setMinimum(0);
	outputProgress->setMaximum(100);
	outputProgress->setValue(0);
	
	outputLayout = new QVBoxLayout();
	outputLayout->addWidget(outputText);
	outputLayout->addWidget(outputProgress);
	outputLayout->setContentsMargins(0, 0, 0, 0);
	outputLayout->setSpacing(5);
	
	widgetOutput->setLayout(outputLayout);
	widgetOutput->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
	
	// Assemble main widget
	// --------------------
	
	layoutMain = new QVBoxLayout();
	layoutMain->addWidget(tabs);
	layoutMain->setContentsMargins(5, 5, 5, 0);
	layoutMain->setSpacing(5);
	
	widgetMain->setLayout(layoutMain);
	
	dockWidgetOutput = new QDockWidget(tr("Pipeline Messages"), this);
	dockWidgetOutput->setAllowedAreas(Qt::AllDockWidgetAreas);
	dockWidgetOutput->setWidget(widgetOutput);
	dockWidgetOutput->setFeatures(QDockWidget::DockWidgetClosable | QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
	dockWidgetOutput->setContentsMargins(5, 5, 5, 5);
	dockWidgetOutput->toggleViewAction()->setText(tr("Show Pipeline Messages"));
	dockWidgetOutput->setVisible(settingsPipeline);
	connect(dockWidgetOutput, SIGNAL(visibilityChanged(bool)), this, SLOT(togglePipeline(bool)));
	this->addDockWidget(Qt::DockWidgetArea(settingsPipelinePosition), dockWidgetOutput);
	
	// Set up actions
	// --------------
	
	actionOpen = new QAction(tr("Open..."), this);
	actionOpen->setShortcuts(QKeySequence::Open);
	actionOpen->setIcon(iconDocumentOpen);
	connect(actionOpen, SIGNAL(triggered()), this, SLOT(loadSettings()));
	
	actionLoadParsetExtragalactic = new QAction(tr("Extragalactic HI Cube"), this);
	connect(actionLoadParsetExtragalactic, SIGNAL(triggered()), this, SLOT(loadParsetExtragalactic()));
	
	actionLoadParsetContinuum = new QAction(tr("Radio Continuum Image"), this);
	connect(actionLoadParsetContinuum, SIGNAL(triggered()), this, SLOT(loadParsetContinuum()));
	
	actionSave = new QAction(tr("Save"), this);
	actionSave->setShortcuts(QKeySequence::Save);
	actionSave->setIcon(iconDocumentSave);
	connect(actionSave, SIGNAL(triggered()), this, SLOT(saveSettings()));
	
	actionSaveAs = new QAction(tr("Save As..."), this);
	actionSaveAs->setShortcuts(QKeySequence::SaveAs);
	actionSaveAs->setIcon(iconDocumentSaveAs);
	connect(actionSaveAs, SIGNAL(triggered()), this, SLOT(saveSettingsAs()));
	
	actionExit = new QAction(tr("Quit"), this);
	actionExit->setShortcuts(QKeySequence::Quit);
	actionExit->setIcon(iconApplicationExit);
	connect(actionExit, SIGNAL(triggered()), this, SLOT(close()));
	
	actionRun = new QAction(tr("Run Pipeline"), this);
	actionRun->setShortcut(Qt::Key_F2);
	actionRun->setIcon(iconDialogOkApply);
	connect(actionRun, SIGNAL(triggered()), this, SLOT(runPipeline()));
	
	actionAbort = new QAction(tr("Abort Pipeline"), this);
	actionAbort->setShortcut(Qt::Key_Escape);
	actionAbort->setIcon(iconDialogClose);
	connect(actionAbort, SIGNAL(triggered()), this, SLOT(pipelineProcessCancel()));
	
	actionDefault = new QAction(tr("New"), this);
	actionDefault->setShortcuts(QKeySequence::New);
	actionDefault->setIcon(iconDocumentNew);
	connect(actionDefault, SIGNAL(triggered()), this, SLOT(resetToDefault()));
	
	actionSaveLogAs = new QAction(tr("Save Messages As..."), this);
	actionSaveLogAs->setEnabled(false);
	actionSaveLogAs->setIcon(iconDocumentSaveAs);
	connect(actionSaveLogAs, SIGNAL(triggered()), this, SLOT(saveLogAs()));
	
	actionClearLog = new QAction(tr("Clear Messages"), this);
	actionClearLog->setEnabled(false);
	actionClearLog->setIcon(iconEditClearList);
	connect(actionClearLog, SIGNAL(triggered()), this, SLOT(clearLog()));
	
	actionShowCatalogue = new QAction(tr("View Catalogue"), this);
	actionShowCatalogue->setEnabled(false);
	actionShowCatalogue->setIcon(iconDocumentPreview);
	connect(actionShowCatalogue, SIGNAL(triggered()), this, SLOT(showCatalogue()));
	
	actionShowCube = new QAction(tr("Data Cube"), this);
	actionShowCube->setEnabled(false);
	connect(actionShowCube, SIGNAL(triggered()), this, SLOT(showCube()));
	
	actionShowFilteredCube = new QAction(tr("Filtered Cube"), this);
	actionShowFilteredCube->setEnabled(false);
	connect(actionShowFilteredCube, SIGNAL(triggered()), this, SLOT(showFilteredCube()));
	
	actionShowNoiseCube = new QAction(tr("Noise Cube"), this);
	actionShowNoiseCube->setEnabled(false);
	connect(actionShowNoiseCube, SIGNAL(triggered()), this, SLOT(showNoiseCube()));
	
	actionShowMask = new QAction(tr("Mask Cube"), this);
	actionShowMask->setEnabled(false);
	connect(actionShowMask, SIGNAL(triggered()), this, SLOT(showMask()));
	
	actionShowMom0 = new QAction(tr("Moment 0"), this);
	actionShowMom0->setEnabled(false);
	connect(actionShowMom0, SIGNAL(triggered()), this, SLOT(showMom0()));
	
	actionShowMom1 = new QAction(tr("Moment 1"), this);
	actionShowMom1->setEnabled(false);
	connect(actionShowMom1, SIGNAL(triggered()), this, SLOT(showMom1()));
	
	actionToolbar = new QAction(tr("Show Toolbar"), this);
	actionToolbar->setCheckable(true);
	actionToolbar->setChecked(settingsToolbar);            // Note that value of settingsToolbar must be known at this point!
	actionToolbar->setEnabled(true);
	connect(actionToolbar, SIGNAL(triggered()), this, SLOT(toggleToolbar()));
	
	actionPipeline = new QAction(tr("Show Pipeline Messages"), this);
	actionPipeline->setCheckable(true);
	actionPipeline->setChecked(settingsPipeline);          // Note that value of settingsPipeline must be known at this point!
	actionPipeline->setEnabled(true);
	connect(actionPipeline, SIGNAL(triggered()), this, SLOT(togglePipeline()));
	
	actionSaveOnExit = new QAction(tr("Save Session on Exit"), this);
	actionSaveOnExit->setCheckable(true);
	actionSaveOnExit->setChecked(settingsSession);         // Note that value of settingsSession must be known at this point!
	actionSaveOnExit->setEnabled(true);
	connect(actionSaveOnExit, SIGNAL(triggered()), this, SLOT(toggleSaveOnExit()));
	
	actionFullScreen = new QAction(tr("Full Screen"), this);
	actionFullScreen->setShortcut(Qt::Key_F11);
	actionFullScreen->setCheckable(true);
	actionFullScreen->setChecked(false);
	actionFullScreen->setEnabled(true);
	actionFullScreen->setIcon(iconFullScreen);
	connect(actionFullScreen, SIGNAL(triggered()), this, SLOT(toggleFullScreen()));
	
	actionHelp = new QAction(tr("SoFiA User Manual"), this);
	actionHelp->setShortcut(QKeySequence::HelpContents);
	actionHelp->setIcon(iconHelpContents);
	connect(actionHelp, SIGNAL(triggered()), this, SLOT(showHandbook()));
	
	actionWhatsThis = QWhatsThis::createAction(this);
	actionWhatsThis->setIcon(iconWhatsThis);
	
	actionAbout = new QAction(tr("About SoFiA"), this);
	actionAbout->setIcon(iconHelpAbout);
	connect(actionAbout, SIGNAL(triggered()), this, SLOT(aboutSoFiA()));
	
	// Set up toolbar
	// --------------
	
	toolBar = new QToolBar(tr("Toolbar"), this);
	toolBar->addAction(actionDefault);
	toolBar->addAction(actionOpen);
	toolBar->addAction(actionSave);
	toolBar->addSeparator();
	toolBar->addAction(actionRun);
	toolBar->addAction(actionAbort);
	toolBar->addAction(actionClearLog);
	toolBar->addSeparator();
	toolBar->addAction(actionWhatsThis);
	
	toolBar->setIconSize(QSize(22, 22));
	toolBar->setMovable(false);
	toolBar->setAllowedAreas(Qt::AllToolBarAreas);
	toolBar->toggleViewAction()->setText(tr("Show Toolbar"));
	toolBar->setVisible(settingsToolbar);
	
	// Set up menu
	// -----------
	
	menuLoadParset = new QMenu(tr("&Parameter Sets"), this);
	menuLoadParset->setIcon(iconDocumentNew);
	menuLoadParset->addAction(actionLoadParsetExtragalactic);
	menuLoadParset->addAction(actionLoadParsetContinuum);
	
	menuFile = new QMenu(tr("&File"), this);
	menuFile->addAction(actionDefault);
	menuFile->addMenu(menuLoadParset);
	menuFile->addSeparator();
	menuFile->addAction(actionOpen);
	menuFile->addSeparator();
	menuFile->addAction(actionSave);
	menuFile->addAction(actionSaveAs);
	menuFile->addSeparator();
	menuFile->addAction(actionExit);
	
	menuPipeline = new QMenu(tr("&Pipeline"), this);
	menuPipeline->addAction(actionRun);
	menuPipeline->addAction(actionAbort);
	menuPipeline->addSeparator();
	menuPipeline->addAction(actionSaveLogAs);
	menuPipeline->addAction(actionClearLog);
	
	menuShowImage = new QMenu(tr("&View Image"), this);
	menuShowImage->setIcon(iconFolderImage);
	menuShowImage->addAction(actionShowCube);
	menuShowImage->addAction(actionShowFilteredCube);
	menuShowImage->addAction(actionShowNoiseCube);
	menuShowImage->addAction(actionShowMask);
	menuShowImage->addAction(actionShowMom0);
	menuShowImage->addAction(actionShowMom1);
	
	menuView = new QMenu(tr("&Analysis"), this);
	menuView->addAction(actionShowCatalogue);
	menuView->addSeparator();
	menuView->addMenu(menuShowImage);
	
	menuSettings = new QMenu(tr("&Settings"), this);
	menuSettings->addAction(actionToolbar);
	menuSettings->addAction(actionPipeline);
	menuSettings->addAction(actionSaveOnExit);
	menuSettings->addSeparator();
	menuSettings->addAction(actionFullScreen);
	
	menuHelp = new QMenu(tr("&Help"), this);
	menuHelp->addAction(actionHelp);
	menuHelp->addAction(actionWhatsThis);
	menuHelp->addSeparator();
	menuHelp->addAction(actionAbout);
	
	this->menuBar()->addMenu(menuFile);
	this->menuBar()->addMenu(menuPipeline);
	this->menuBar()->addMenu(menuView);
	this->menuBar()->addMenu(menuSettings);
	this->menuBar()->addMenu(menuHelp);
	
	// Set up status bar
	// -----------------
	
	this->statusBar()->setSizeGripEnabled(false);
	this->statusBar()->showMessage("");
	
	// Set up main window
	// ------------------
	
	this->addToolBar(Qt::TopToolBarArea, toolBar);
	this->setCentralWidget(widgetMain);
	this->setContextMenuPolicy(Qt::NoContextMenu);    // This will prevent automatic context menus that allow hiding parts of the UI.
	this->resize(settingsWindowWidth, settingsWindowHeight);
	this->setWindowIcon(iconSoFiA);
	updateWindowTitle();
	
	return;
}



// ----------------------------------------
// Function to create What's This? entries:
// ----------------------------------------

void SoFiA::createWhatsThis()
{
	tabSourceFindingFieldMaxScale->setWhatsThis(tr("<h3>CNHI.maxScale</h3><p>The maximum size of test regions. The default value of <strong>-1</strong> is a flag value that sets the maximum size to half of the size of the spectral axis.</p>"));
	tabSourceFindingMedianTest->setWhatsThis(tr("<h3>CNHI.medianTest</h3><p>This parameter determines whether test regions need to have a median greater than that of the remaining data in order to be considered a source.</p>"));
	tabSourceFindingFieldMinScale->setWhatsThis(tr("<h3>CNHI.minScale</h3><p>The minimum size of test regions.</p>"));
	tabSourceFindingFieldPReq->setWhatsThis(tr("<h3>CNHI.pReq</h3><p>Minimum probability requirement for detections to be considered genuine. Sensible values are typically in the range of about 10<sup>&minus;3</sup> to 10<sup>&minus;5</sup>.</p>"));
	tabSourceFindingFieldQReq->setWhatsThis(tr("<h3>CNHI.qReq</h3><p>This is the Q value of the Kuiper test, which is a heuristic parameter for assessing the quality/accuracy of the probability calculated from the Kuiper test. The minimum scale that the CNHI source finder processes is increased until it is sufficiently large to ensure that the required Q value is achieved. This requirement supersedes the user-specified minimum scale (see parameter <strong>CNHI.minScale</strong>). The default value is <strong>3.8</strong>, which is the generally accepted minimally useful value.</p>"));
	tabSourceFindingFieldVerbose->setWhatsThis(tr("<h3>CNHI.verbose</h3><p>An integer value that indicates the level of verbosity of the CNHI finder. Values of <strong>0</strong>, <strong>1</strong> and <strong>2</strong> correspond to <strong>none</strong>, <strong>minimal</strong> and <strong>maximal</strong>, respectively.</p>"));
	tabInputFieldFlags->setWhatsThis(tr("<h3>flag.regions</h3><p>Pixel/channel range(s) to be flagged prior to source finding. Format:</p><p><code>[[x1, x2, y1, y2, z1, z2], ...]</code></p><p>A place holder, <code>''</code> (two single quotes), can be used for the upper range limit (<code>x2</code>, <code>y2</code>, and <code>z2</code>) to flag all the way to the end, e.g.</p><p><code>[[0, '', 0, '', 0, 19]]</code></p><p>will flag the first 20 channels of the entire cube. The default is an empty list, <code>[]</code>, which means to not flag anything.</p>"));
	tabInputFieldFlagCube->setWhatsThis(tr("<h3>flag.file</h3><p>Full path and file name of an optional flagging cube that must be of the same dimensions as the input cube. All pixels that are flagged (i.e. set to BLANK) in the flagging cube will also be flagged in the input cube prior to source finding. The default is to not apply flags.</p>"));
	tabInputFieldData->setWhatsThis(tr("<h3>import.inFile</h3><p>Full path and file name of the input data cube. This option is mandatory, and there is no default. Note that only <b>FITS</b> files are currently supported by SoFiA.</p>"));
	tabInputFieldMask->setWhatsThis(tr("<h3>import.maskFile</h3><p>Full path and file name of an optional file containing a mask of pixels identified as part of a source, e.g. from a previous run of SoFiA. This can be used to re-parametrise sources without repeating the source finding step or to add more sources from a second source finding run. The default is to not read a mask cube.</p>"));
	tabInputFieldSubcube->setWhatsThis(tr("<h3>import.subcube</h3><p>This parameter defines a subcube to be read in and processed by SoFiA. Depending on the value of import.subcubeMode, the range is either specified in pixels as</p><p><code>[x1, x2, y1, y2, z1, z2]</code></p><p>or in world coordinates as</p><p><code>[x, y, z, rx, ry, rz]</code></p><p>In the latter case, <code>x</code>, <code>y</code> and <code>z</code> define the centre of the subcube, and <code>rx</code>, <code>ry</code> and <code>rz</code> specify the half-widths in the three dimensions. If world coordinates are used, all parameters must be in the native format as defined in the header of the data cube; e.g. if <code>CUNIT3</code> is <code>'Hz'</code> then both <code>z</code> and <code>rz</code> must be given in hertz. The default is an empty list, <code>[]</code>, which means to read the entire cube.</p>"));
	tabInputFieldSubcubeMode->setWhatsThis(tr("<h3>import.subcubeMode</h3><p>This parameter defines whether import.subcube is specified in pixels (<b>pixel</b>) or in world coordinates (<b>world</b>).</p>"));
	tabInputFieldWeights->setWhatsThis(tr("<h3>import.weightsFile</h3><p>Full path and file name of an optional file containing the weights of pixels in the input cube. The weights will be applied before running the source finder. The default is to not apply weights.</p>"));
	tabInputFieldWeightsFunction->setWhatsThis(tr("<h3>import.weightsFunction</h3><p>Analytic function used to describe the data weights as a function of the three cube dimensions, <b>x</b>, <b>y</b>, and <b>z</b>. The default is to not apply weights. The following mathematical functions from Numpy are supported:<p><p style=\"font-family:monospace;\">sin(), cos(), tan(), arcsin(), arccos(), arctan(), arctan2(), sinh(), cosh(), tanh(), arcsinh(), arccosh(), arctanh(), exp(), log(), log10(), log2(), sqrt(), square(), power(), absolute(), fabs(), sign()</p><p>Note that the weights function is not applied whenever a weights cube is specified (see <b>import.weightsFile</b>).</p>"));
	tabMergingFieldMinSizeX->setWhatsThis(tr("<h3>merge.minSizeX</h3><p>Minimum required extent of genuine sources in first dimension. Sources smaller than this will be discarded.</p>"));
	tabMergingFieldMinSizeY->setWhatsThis(tr("<h3>merge.minSizeY</h3><p>Minimum required extent of genuine sources in second dimension. Sources smaller than this will be discarded.</p>"));
	tabMergingFieldMinSizeZ->setWhatsThis(tr("<h3>merge.minSizeZ</h3><p>Minimum required extent of genuine sources in third dimension. Sources smaller than this will be discarded.</p>"));
	tabMergingFieldRadiusX->setWhatsThis(tr("<h3>merge.radiusX</h3><p>Merging radius in first dimension in pixels.</p>"));
	tabMergingFieldRadiusY->setWhatsThis(tr("<h3>merge.radiusY</h3><p>Merging radius in second dimension in pixels.</p>"));
	tabMergingFieldRadiusZ->setWhatsThis(tr("<h3>merge.radiusZ</h3><p>Merging radius in third dimension in pixels.</p>"));
	tabMergingButtonPositivity->setWhatsThis(tr("<h3>merge.positivity</h3><p>Discard all negative signals and only merge positive signals into detections.</p>"));
	tabInputFieldCatalog->setWhatsThis(tr("<h3>optical.sourceCatalogue</h3><p>This defines the full path and file name of the input catalogue required for catalogue-based source finding (see parameter <b>steps.doOptical</b>). There is no default.</p><p>Catalogues must be comma-separated and contain at least four columns containing a unique ID number, right ascension, declination and frequency/velocity of the positions to be searched. All parameters must be specified in the native WCS units of the cube. In addition, a header line must be provided, with the four parameter columns above listed as <code>id</code>, <code>ra</code>, <code>dec</code> and <code>z</code>.</p>"));
	tabInputFieldSpatialSize->setWhatsThis(tr("<h3>optical.spatSize</h3><p>This defines the <b>spatial</b> size of the sub-cube to be searched around each catalogue position. The size must be specified in the <b>native units</b> of the data cube, e.g. in degrees.</p>"));
	tabInputFieldSpectralSize->setWhatsThis(tr("<h3>optical.specSize</h3><p>This defines the <b>spectral</b> size of the sub-cube to be searched around each catalogue position. The size must be specified in the <b>native units</b> of the data cube, e.g. in km/s or Hz.</p>"));
	tabInputFieldMultiCat->setWhatsThis(tr("<h3>optical.storeMultiCat</h3><p>If set to <b>true</b>, a separate output catalogue will be created for each input position, containing only the sources found in that subcube. In addition, a single, merged catalogue will also be created. By default this parameter is set to <b>false</b>, in which case only a single output catalogue file is generated.</p>"));
	tabParametrisationButtonDilateMask->setWhatsThis(tr("<h3>parameters.dilateMask</h3><p>Run the mask optimisation algorithm based on spatially <b>dilating</b> the initial mask to achieve more accurate flux measurements.</p>"));
	tabParametrisationButtonBusyFunction->setWhatsThis(tr("<h3>parameters.fitBusyFunction</h3><p>Fit the Busy Function (<a href=\"http://adsabs.harvard.edu/abs/2014MNRAS.438.1176W\">Westmeier et al. 2014</a>) to the integrated spectrum for more accurate parameterisation.</p>"));
	tabParametrisationButtonMaskOpt->setWhatsThis(tr("<h3>parameters.optimiseMask</h3><p>Run the mask optimisation algorithm based on fitting and <b>growing</b> ellipses to achieve more accurate flux measurements.</p>"));
	tabParametrisationFieldScaleKernel->setWhatsThis(tr("<h3>reliability.scaleKernel</h3><p>If <b>autoKernel</b> is set to <b>true</b>, then this parameter will determine whether the kernel size will be set automatically by SoFiA (<b>auto</b>) or scaled by the given factor (<b>&gt;&nbsp;0</b>). If <b>autoKernel</b> is <b>false</b>, this option will be ignored."));
	tabParametrisationSliderScaleKernel->setWhatsThis(tr("<h3>reliability.scaleKernel</h3><p>If <b>autoKernel</b> is set to <b>true</b>, then this parameter will determine whether the kernel size will be set automatically by SoFiA (<b>auto</b>) or scaled by the given factor (<b>&gt;&nbsp;0</b>). If <b>autoKernel</b> is <b>false</b>, this option will be ignored."));
	tabParametrisationButtonRelPlot->setWhatsThis(tr("<h3>reliability.makePlot</h3><p>If set to <b>true</b>, PDF files showing the distribution of positive and negative detections in parameter space will be created for diagnostic purposes.</p>"));
	tabParametrisationFieldRelMin->setWhatsThis(tr("<h3>reliability.threshold</h3><p>Discard sources whose reliability is below this threshold.</p>"));
	tabInFilterFieldEdgeX->setWhatsThis(tr("<h3>scaleNoise.edgeX</h3><p>Size of edge (in pixels) to be excluded in first coordinate. The setting will be ignored if <b>scaleNoise.method = local</b>.</p>"));
	tabInFilterFieldEdgeY->setWhatsThis(tr("<h3>scaleNoise.edgeY</h3><p>Size of edge (in pixels) to be excluded in second coordinate. The setting will be ignored if <b>scaleNoise.method = local</b>.</p>"));
	tabInFilterFieldEdgeZ->setWhatsThis(tr("<h3>scaleNoise.edgeZ</h3><p>Size of edge (in pixels) to be excluded in third coordinate. The setting will be ignored if <b>scaleNoise.method = local</b>.</p>"));
	tabInFilterFieldFluxRange->setWhatsThis(tr("<h3>scaleNoise.fluxRange</h3><p>Range of flux values to be used in noise measurement. Can be <b>negative</b>, <b>positive</b> or <b>all</b> to use only negative, only positive or all pixels, respectively.</p>"));
	tabInFilterFieldInterpolation->setWhatsThis(tr("<h3>scaleNoise.interpolation</h3><p>Select if the local noise measurement made by SoFiA should be interpolated in between grid points or not. The following settings are allowed:</p><ul><li><b>none</b> will simply fill each grid cell with the noise value. This is the default setting.</li><li><b>linear</b> will carry out a linear interpolation between grid points, which is slower, but results in smoother noise maps.</li><li><b>cubic</b> will carry out a cubic spline interpolation between grid points, which is by far the slowest method, but results in the smoothest noise maps.</li></ul><p>Note that this setting will only be relevant if <b>scaleNoise.method = local</b>.</p>"));
	//tabInFilterFieldGridSpatial->setWhatsThis(tr("<h3>scaleNoise.gridSpatial</h3><p>This defines the spatial grid size on which the local noise measurement takes place. The setting will be ignored if <b>scaleNoise.method = global</b>. It must be an even number of 2 or greater.</p>"));
	//tabInFilterFieldGridSpectral->setWhatsThis(tr("<h3>scaleNoise.gridSpectral</h3><p>This defines the spectral grid size on which the local noise measurement takes place. The setting will be ignored if <b>scaleNoise.method = global</b>. It must be an even number of 2 or greater.</p>"));
	tabInFilterFieldMethod->setWhatsThis(tr("<h3>scaleNoise.method</h3><p>If set to <b>global</b>, the noise measurement will be carried out on the entire projected image plane perpendicular to the axis along which the noise is to be scaled (default). If set to <b>local</b>, the noise measurement will occur in a running window of specified size on a specified grid. Note that the latter can be slow and memory-heavy and is only recommended for small cubes or 2D images that are affected by localised noise variations.</p>"));
	tabInFilterFieldScaleX->setWhatsThis(tr("<h3>scaleNoise.scaleX</h3><p>Noise normalisation in first (spatial) dimension. The setting will be ignored if <b>scaleNoise.method = local</b>.</p>"));
	tabInFilterFieldScaleY->setWhatsThis(tr("<h3>scaleNoise.scaleY</h3><p>Noise normalisation in second (spatial) dimension. The setting will be ignored if <b>scaleNoise.method = local</b>.</p>"));
	tabInFilterFieldScaleZ->setWhatsThis(tr("<h3>scaleNoise.scaleZ</h3><p>Noise normalisation in third (spectral) dimension. The setting will be ignored if <b>scaleNoise.method = local</b>.</p>"));
	tabInFilterFieldStatistic->setWhatsThis(tr("<h3>scaleNoise.statistic</h3><p>Statistic used to measure the noise. This can be median absolute deviation (<b>mad</b>), standard deviation (<b>std</b>), Gaussian fit to flux histogram (<b>gauss</b>) or Gaussian fit to negative fluxes (<b>negative</b>).</p>"));
	tabInFilterFieldWindowSpatial->setWhatsThis(tr("<h3>scaleNoise.windowSpatial</h3><p>This defines the spatial window size over which the local noise measurement is taken. The setting will only be relevant if <b>scaleNoise.method = local</b>. It must be an odd number of 1 or greater.</p>"));
	tabInFilterFieldWindowSpectral->setWhatsThis(tr("<h3>scaleNoise.windowSpectral</h3><p>This defines the spectral window size over which the local noise measurement is taken. The setting will only be relevant if <b>scaleNoise.method = local</b>. It must be an odd number of 1 or greater.</p>"));
	tabSourceFindingFieldEdgeMode->setWhatsThis(tr("<h3>SCfind.edgeMode</h3><p>Behaviour near the edge of the cube. The following values are possible:<p><ul><li><b>constant:</b> assume constant value of 0</li><li><b>nearest:</b> assume constant value equal to edge pixel</li><li><b>reflect:</b> mirror values at edge, thereby including the edge pixel itself</li><li><b>mirror:</b> mirror values at position of outermost pixel, thereby excluding the edge pixel itself</li><li><b>wrap:</b> copy values from opposite edge of the array</li></ul>"));
	tabSourceFindingFieldKernels->setWhatsThis(tr("<h3>SCfind.kernels</h3><p>List of kernels to be used for smoothing. The format is:</p><p style=\"font-family:monospace;\">[[dx, dy, dz, 'type'], ...]</p><p>where <b>dx</b>, <b>dy</b>, and <b>dz</b> are the spatial and spectral kernel sizes (FWHM), and <b>'type'</b> can be boxcar (<b>'b'</b>) or Gaussian (<b>'g'</b>). Note that 'type' only applies to the spectral axis, and the spatial kernel is always Gaussian.</p>"));
	tabSourceFindingFieldKunit->setWhatsThis(tr("<h3>SCfind.kernelUnit</h3><p>Are kernel parameters specified in <b>pixel</b> or <b>world</b> coordinates?</p>"));
	tabSourceFindingFieldRmsMode->setWhatsThis(tr("<h3>SCfind.rmsMode</h3><p>Noise determination method: Gaussian fit to flux histogram (<b>gauss</b>),Gaussian fit to negative flux histogram (<b>negative</b>), median absolute deviation (<b>mad</b>), or standard deviation (<b>std</b>).</p>"));
	tabSourceFindingFieldFluxRange->setWhatsThis(tr("<h3>SCfind.fluxRange</h3><p>Range of flux values to be used in noise measurement. Can be <b>negative</b>, <b>positive</b> or <b>all</b> to use only negative, only positive or all pixels, respectively.</p>"));
	tabSourceFindingFieldThreshold->setWhatsThis(tr("<h3>SCfind.threshold</h3><p>Flux threshold relative to the noise level.</p>"));
	tabInFilterFieldBorder->setWhatsThis(tr("<h3>smooth.edgeMode</h3><p>Behaviour near the edge of the cube. The following options are supported:</p><ul><li><b>constant:</b> assume constant value of 0</li><li><b>nearest:</b> assume constant value equal to edge pixel</li><li><b>reflect:</b> mirror values at edge, thereby including the edge pixel itself</li><li><b>mirror:</b> mirror values at position of outermost pixel, thereby excluding the edge pixel itself</li><li><b>wrap:</b> copy values from opposite edge of the array</li></ul>"));
	tabInFilterFieldKernel->setWhatsThis(tr("<h3>smooth.kernel</h3><p>Type of smoothing kernel used in both spatial and spectral smoothing. Can be <b>gaussian</b>, <b>boxcar</b> or <b>median</b>.</p>"));
	tabInFilterFieldSmoothingSpatialLon->setWhatsThis(tr("<h3>smooth.kernelX</h3><p>Kernel size in pixels for first coordinate. For Gaussian kernels the value refers to the FWHM.</p>"));
	tabInFilterFieldSmoothingSpatialLat->setWhatsThis(tr("<h3>smooth.kernelY</h3><p>Kernel size in pixels for second coordinate. For Gaussian kernels the value refers to the FWHM.</p>"));
	tabInFilterFieldSmoothingSpectral->setWhatsThis(tr("<h3>smooth.kernelZ</h3><p>Kernel size in pixels for third coordinate. For Gaussian kernels the value refers to the FWHM.</p>"));
	tabSourceFindingGroupBox3->setWhatsThis(tr("<h3>steps.doCNHI</h3><p>Run the Characterised Noise HI (CNHI) source finder (<a href=\"http://adsabs.harvard.edu/abs/2012PASA...29..251J\">Jurek 2012</a>).</p>"));
	tabOutputButtonCubelets->setWhatsThis(tr("<h3>steps.doCubelets</h3><p>Create and save data products for each individual source, including sub-cubes, moment 0, 1 and 2 maps, integrated spectra and position&ndash;velocity diagrams.</p>"));
	tabInputGroupBox3->setWhatsThis(tr("<h3>steps.doFlag</h3><p>Flag pixel and channel ranges before proceeding. Details are specified with the <b>flag.file</b> and <b>flag.regions</b> options.</p>"));
	tabMergingGroupBox1->setWhatsThis(tr("<h3>steps.doMerge</h3><p>Merge detected pixels into final sources.</p>"));
	tabOutputButtonMom0->setWhatsThis(tr("<h3>steps.doMom0</h3><p>Create and save moment-0 map of all detections.</p>"));
	tabOutputButtonMom1->setWhatsThis(tr("<h3>steps.doMom1</h3><p>Create and save moment-1 map of all detections.</p>"));
	tabInputGroupBox2->setWhatsThis(tr("<h3>steps.doOptical</h3><p>Run SoFiA on multiple, smaller sub-cubes centred on positions defined in an input source catalogue. A catalogue file will need to be specified (see parameter <b>optical.sourceCatalogue</b>). This could, e.g., be an optical galaxy catalogue with the aim to search for HI detections at the positions of all galaxies.</p>"));
	tabParametrisationGroupBox1->setWhatsThis(tr("<h3>steps.doParameterise</h3><p>Run the mask optimisation and source parameterisation module.</p>"));
	tabParametrisationGroupBox2->setWhatsThis(tr("<h3>steps.doReliability</h3><p>Use negative detections to calculate the reliability of each source by comparing the distribution of positive and negative detections in parameter space (<a href=\"http://adsabs.harvard.edu/abs/2012PASA...29..296S\">Serra et al. 2012</a>). Note that this requires a flux threshold low enough to ensure that a substantial number of positive and negative noise peaks are detected.</p>"));
	tabInFilterGroupBox2->setWhatsThis(tr("<h3>steps.doScaleNoise</h3><p>Automatically normalise noise levels prior to source finding.</p>"));
	tabSourceFindingGroupBox1->setWhatsThis(tr("<h3>steps.doSCfind</h3><p>Run the smooth + clip finder on the data cube.</p>"));
	tabInFilterGroupBox1->setWhatsThis(tr("<h3>steps.doSmooth</h3><p>Spatially and spectrally smooth cube prior to source finding.</p>"));
	tabInputGroupBox4->setWhatsThis(tr("<h3>steps.doSubcube</h3><p>If set to true, source finding will be carried out on a subcube to be defined by the <b>import.subcube</b> and <b>import.subcubeMode</b> options.</p>"));
	tabSourceFindingGroupBox2->setWhatsThis(tr("<h3>steps.doThreshold</h3><p>Run the threshold finder on the data cube.</p>"));
	tabInFilterGroupBox3->setWhatsThis(tr("<h3>steps.doWavelet</h3><p>Decompose the data cube into wavelet components using the 2D&ndash;1D wavelet decomposition algorithm of <a href=\"http://adsabs.harvard.edu/abs/2012PASA...29..244F\">Fl&ouml;er et al. (2012)</a>.</p>"));
	tabOutputButtonFilteredCube->setWhatsThis(tr("<h3>steps.doWriteFilteredCube</h3><p>Save a copy of the filtered data cube. Note that this will only make sense if at least one of the input filters was applied.</p>"));
	tabOutputButtonNoiseCube->setWhatsThis(tr("<h3>steps.doWriteNoiseCube</h3><p>If set to <b>true</b>, then a copy of the noise cube estimated by SoFiA will be written to disk as a FITS file. This will only work if noise scaling is enabled.</p>"));
	tabOutputButtonMask->setWhatsThis(tr("<h3>steps.doWriteMask</h3><p>Save source mask cube.</p>"));
	tabSourceFindingFieldClipMethod->setWhatsThis(tr("<h3>threshold.clipMethod</h3><p>Define whether the flux threshold is <b>relative</b> to the noise level or in <b>absolute</b> flux units.</p>"));
	tabSourceFindingFieldRmsMode2->setWhatsThis(tr("<h3>threshold.rmsMode</h3><p>Noise determination method: Gaussian fit to flux histogram (<b>gauss</b>), Gaussian fit to negative flux histogram (<b>negative</b>), median absolute deviation (<b>mad</b>), or standard deviation (<b>std</b>).</p>"));
	tabSourceFindingFieldFluxRange2->setWhatsThis(tr("<h3>threshold.fluxRange</h3><p>Range of flux values to be used in noise measurement. Can be <b>negative</b>, <b>positive</b> or <b>all</b> to use only negative, only positive or all pixels, respectively.</p>"));
	tabSourceFindingFieldThreshold2->setWhatsThis(tr("<h3>threshold.threshold</h3><p>Absolute or relative flux threshold for detections (see <b>threshold.clipMethod</b>).</p>"));
	tabInFilterField2d1dIterations->setWhatsThis(tr("<h3>wavelet.iterations</h3><p>Number of iterations in the reconstruction process.</p>"));
	tabInFilterField2d1dPositivity->setWhatsThis(tr("<h3>wavelet.positivity</h3><p>If <b>true</b>, include only positive wavelet components in the decomposition.</p>"));
	tabInFilterField2d1dScaleXY->setWhatsThis(tr("<h3>wavelet.scaleXY</h3><p>Number of spatial scales used in the decomposition. The default value of <b>-1</b> will automatically determine the appropriate number of scales based on the data cube.</p>"));
	tabInFilterField2d1dScaleZ->setWhatsThis(tr("<h3>wavelet.scaleZ</h3><p>Number of spectral scales used in the decomposition. The default value of <b>-1</b> will automatically determine the appropriate number of scales based on the data cube.</p>"));
	tabInFilterField2d1dThreshold->setWhatsThis(tr("<h3>wavelet.threshold</h3><p>Flux threshold used in the wavelet reconstruction in multiples of the rms noise. Note that this threshold only determines which wavelet components are added to the decomposed cube; any source finding will be done separately, using a different flux threshold.</p>"));
	tabOutputFieldBaseName->setWhatsThis(tr("<h3>writeCat.baseName</h3><p>By default, SoFiA will use the file name of the input data cube as the template for all output data file names, usually extended by a specific postfix based on the actual output data product. By setting the <b>writeCat.basename</b> parameter, a different template for all output file names can be specified.</p>"));
	tabOutputButtonOverwrite->setWhatsThis(tr("<h3>writeCat.overwrite</h3><p>If <b>true</b>, overwrite any existing output files. If <b>false</b>, print an error message for any output file that already exists and skip writing that file.</p>"));
	tabOutputButtonCompress->setWhatsThis(tr("<h3>writeCat.compress</h3><p>If set to <b>true</b>, use <a href=\"http://www.gzip.org/\">gzip</a> to compress all output files.</p>"));
	tabOutputFieldDirectory->setWhatsThis(tr("<h3>writeCat.outputDir</h3><p>Optional directory path to which all output files are written. If not specified, the directory of the input cube will be used by default.</p>"));
	tabOutputFieldParameters->setWhatsThis(tr("<h3>writeCat.parameters</h3><p>List of parameters to appear in source catalogue. Format:</p><p>['par1', 'par2', ...]</p><p>An asterisk, <span style=\"font-family:monospace;\">['*']</span>, means that all available parameters are written to the catalogue (default behaviour).</p><p> Note that, depending on the actual pipeline settings, some selected parameters may not appear in the source catalogue.</p>"));
	tabOutputButtonASCII->setWhatsThis(tr("<h3>writeCat.writeASCII</h3><p>Write catalogue in ASCII format.</p>"));
	tabOutputButtonSQL->setWhatsThis(tr("<h3>writeCat.writeSQL</h3><p>Write catalogue in SQL format.</p>"));
	tabOutputButtonXML->setWhatsThis(tr("<h3>writeCat.writeXML</h3><p>Write catalogue in VO table (XML) format.</p>"));
	
	return;
}



// ------------------------------
// Support drag and drop of files
// ------------------------------

void SoFiA::dragMoveEvent(QDragMoveEvent *event)
{
	event->accept();
	return;
}

void SoFiA::dragEnterEvent(QDragEnterEvent *event)
{
	event->acceptProposedAction();
	return;
}

void SoFiA::dropEvent(QDropEvent *event)
{
	const QMimeData* mimeData = event->mimeData();
	
	// Check for single file:
	if(mimeData->hasUrls())
	{
		QString path;
		QList<QUrl> urlList = mimeData->urls();
		
		if(urlList.size() == 1)
		{
			path.append(urlList.at(0).toLocalFile());
			
			this->loadFile(path);
		}
	}
	
	return;
}




// --------------------------------
// Reimplementation of closeEvent()
// --------------------------------

void SoFiA::closeEvent(QCloseEvent *event)
{
	if(pipelineProcess->state() == QProcess::NotRunning)
	{
		if(settingsChanged)
		{
			// Unsaved changes found, ask the user what to do:
			QMessageBox messageBox(this);
			messageBox.setWindowTitle(tr("SoFiA - Quit"));
			messageBox.setText(tr("<p>This action will terminate SoFiA. Any unsaved changes may be lost.</p><p>Do you wish to exit from SoFiA?</p>"));
			messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
			messageBox.setDefaultButton(QMessageBox::Ok);
			messageBox.setIcon(QMessageBox::Warning);
			int choice = messageBox.exec();
			
			if(choice == QMessageBox::Ok)
			{
				shutdownInitiated = true;
				event->accept();
			}
			else event->ignore();
		}
		else
		{
			// No unsaved changes found, quit SoFiA:
			shutdownInitiated = true;
			event->accept();
		}
	}
	else
	{
		// Pipeline still running, don’t quit:
		event->ignore();
		
		QString messageText = tr("<p>The pipeline is still running.</p><p>If you wish to exit from SoFiA, you will need to either manually abort the pipeline or wait until the pipeline run has finished.</p>");
		QString statusText = tr("Pipeline still running.");
		showMessage(MESSAGE_WARNING, messageText, statusText);
	}
	
	return;
}
