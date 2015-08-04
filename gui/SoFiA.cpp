/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 0.4.0 (SoFiA.cpp) - Source Finding Application                 ///
/// Copyright (C) 2013-2014 Tobias Westmeier                             ///
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
#include "SoFiA.h"



// -----------
// Constructor
// -----------

SoFiA::SoFiA(int argc, char *argv[])
{
    pipelineProcess = new QProcess(this);
    connect(pipelineProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(pipelineProcessReadStd()));
    connect(pipelineProcess, SIGNAL(readyReadStandardError()), this, SLOT(pipelineProcessReadErr()));
    connect(pipelineProcess, SIGNAL(started()), this, SLOT(pipelineProcessStarted()));
    connect(pipelineProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(pipelineProcessFinished(int, QProcess::ExitStatus)));
    connect(pipelineProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(pipelineProcessError(QProcess::ProcessError)));
    
    this->createInterface();
    this->setDefaults();
    this->setAcceptDrops(true);        // Enable drag and drop
    
    // Create catalogue viewer window, but don't show it yet:
    spreadsheet = new WidgetSpreadsheet(this);
    spreadsheet->hide();
    
    // Load parameter file if specified:
    if(argc > 1)
    {
        QString fileName = QString(argv[1]);
        
        if(!fileName.isEmpty() and fileName[0] != '/')
        {
            fileName.prepend('/');
            fileName.prepend(QDir::currentPath());
        }
        
        if(loadFile(fileName))
        {
            QString messageText = tr("<p>Failed to read input file %1.</p>").arg(fileName.section('/', -1));
            QString statusText = tr("Failed to read input file %1.").arg(fileName.section('/', -1));
            showMessage(MESSAGE_ERROR, messageText, statusText);
        }
        
        //std::cout << filename.toLocal8Bit().constData() << '\n';
    }
    
    return;
}



// ----------
// Destructor
// ----------

SoFiA::~SoFiA()
{
    if(QFile::exists(SOFIA_TEMP_FILE))
    {
        // Temporary parameter file present, needs to be deleted
        if(QFile::remove(SOFIA_TEMP_FILE) == false)
        {
            std::cerr << "Warning: Failed to remove temporary parameter file on exit.\n";
        }
    }
    
    return;
}



// ----------------------------
// Function to display messages
// ----------------------------

int SoFiA::showMessage(int severity, QString &messageText, QString &statusText)
{
    if(severity < MESSAGE_INFO or severity > MESSAGE_ERROR) severity = MESSAGE_INFO;
    
    if(!statusText.isEmpty())
    {
        QString statusTitle[3] = {tr("Information"), tr("Warning"), tr("Error")};
        this->statusBar()->showMessage(QString("%1: %2").arg(statusTitle[severity]).arg(statusText));
    }
    
    if(!messageText.isEmpty())
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

int SoFiA::setDefaults()
{
    QString fileName = SOFIA_DEFAULT_SETTINGS;
    
    if(loadFile(fileName))
    {
        QString messageText = tr("<p>Failed to load default parameters.</p><p>Please close the programme and check your installation. SoFiA will not function properly without the default parameters.</p>");
        QString statusText = tr("Failed to load default parameters.");
        showMessage(MESSAGE_ERROR, messageText, statusText);
        
        return 1;
    }
    
    currentFileName.clear();
    
    QString messageText = QString("");
    QString statusText = tr("Parameters reset to default.");
    showMessage(MESSAGE_INFO, messageText, statusText);
    
    this->setWindowTitle(tr("SoFiA"));
    
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
            parameters.insert(w->objectName(), w->text());
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
    
    // Treat list of output parameters separately and check respective check boxes:
    QString listOutputPar;
    
    if(not tabOutputGroupBox2->isChecked())
    {
        listOutputPar = QString("\'*\'");
    }
    else
    {
        if(tabOutputButtonParameter_id->isChecked())        listOutputPar.append(QString("\'id\',"));
        if(tabOutputButtonParameter_x_geo->isChecked())     listOutputPar.append(QString("\'x_geo\',"));
        if(tabOutputButtonParameter_y_geo->isChecked())     listOutputPar.append(QString("\'y_geo\',"));
        if(tabOutputButtonParameter_z_geo->isChecked())     listOutputPar.append(QString("\'z_geo\',"));
        if(tabOutputButtonParameter_x->isChecked())         listOutputPar.append(QString("\'x\',"));
        if(tabOutputButtonParameter_y->isChecked())         listOutputPar.append(QString("\'y\',"));
        if(tabOutputButtonParameter_z->isChecked())         listOutputPar.append(QString("\'z\',"));
        if(tabOutputButtonParameter_x_min->isChecked())     listOutputPar.append(QString("\'x_min\',"));
        if(tabOutputButtonParameter_x_max->isChecked())     listOutputPar.append(QString("\'x_max\',"));
        if(tabOutputButtonParameter_y_min->isChecked())     listOutputPar.append(QString("\'y_min\',"));
        if(tabOutputButtonParameter_y_max->isChecked())     listOutputPar.append(QString("\'y_max\',"));
        if(tabOutputButtonParameter_z_min->isChecked())     listOutputPar.append(QString("\'z_min\',"));
        if(tabOutputButtonParameter_z_max->isChecked())     listOutputPar.append(QString("\'z_max\',"));
        if(tabOutputButtonParameter_n_pix->isChecked())     listOutputPar.append(QString("\'n_pix\',"));
        if(tabOutputButtonParameter_snr_min->isChecked())   listOutputPar.append(QString("\'snr_min\',"));
        if(tabOutputButtonParameter_snr_max->isChecked())   listOutputPar.append(QString("\'snr_max\',"));
        if(tabOutputButtonParameter_snr_sum->isChecked())   listOutputPar.append(QString("\'snr_sum\',"));
        if(tabOutputButtonParameter_n_pos->isChecked())     listOutputPar.append(QString("\'n_pos\',"));
        if(tabOutputButtonParameter_n_neg->isChecked())     listOutputPar.append(QString("\'n_neg\',"));
        if(tabOutputButtonParameter_rel->isChecked())       listOutputPar.append(QString("\'rel\',"));
        if(tabOutputButtonParameter_bf_a->isChecked())      listOutputPar.append(QString("\'bf_a\',"));
        if(tabOutputButtonParameter_bf_b1->isChecked())     listOutputPar.append(QString("\'bf_b1\',"));
        if(tabOutputButtonParameter_bf_b2->isChecked())     listOutputPar.append(QString("\'bf_b2\',"));
        if(tabOutputButtonParameter_bf_c->isChecked())      listOutputPar.append(QString("\'bf_c\',"));
        if(tabOutputButtonParameter_bf_chi2->isChecked())   listOutputPar.append(QString("\'bf_chi2\',"));
        if(tabOutputButtonParameter_bf_flag->isChecked())   listOutputPar.append(QString("\'bf_flag\',"));
        if(tabOutputButtonParameter_bf_f_int->isChecked())  listOutputPar.append(QString("\'bf_f_int\',"));
        if(tabOutputButtonParameter_bf_f_peak->isChecked()) listOutputPar.append(QString("\'bf_f_peak\',"));
        if(tabOutputButtonParameter_bf_w->isChecked())      listOutputPar.append(QString("\'bf_w\',"));
        if(tabOutputButtonParameter_bf_w20->isChecked())    listOutputPar.append(QString("\'bf_w20\',"));
        if(tabOutputButtonParameter_bf_w50->isChecked())    listOutputPar.append(QString("\'bf_w50\',"));
        if(tabOutputButtonParameter_bf_xe->isChecked())     listOutputPar.append(QString("\'bf_xe\',"));
        if(tabOutputButtonParameter_bf_xp->isChecked())     listOutputPar.append(QString("\'bf_xp\',"));
        if(tabOutputButtonParameter_bf_z->isChecked())      listOutputPar.append(QString("\'bf_z\',"));
        if(tabOutputButtonParameter_ell_maj->isChecked())   listOutputPar.append(QString("\'ell_maj\',"));
        if(tabOutputButtonParameter_ell_min->isChecked())   listOutputPar.append(QString("\'ell_min\',"));
        if(tabOutputButtonParameter_ell_pa->isChecked())    listOutputPar.append(QString("\'ell_pa\',"));
        if(tabOutputButtonParameter_f_peak->isChecked())    listOutputPar.append(QString("\'f_peak\',"));
        if(tabOutputButtonParameter_f_int->isChecked())     listOutputPar.append(QString("\'f_int\',"));
        if(tabOutputButtonParameter_f_wm50->isChecked())    listOutputPar.append(QString("\'f_wm50\',"));
        if(tabOutputButtonParameter_rms->isChecked())       listOutputPar.append(QString("\'rms\',"));
        if(tabOutputButtonParameter_w20->isChecked())       listOutputPar.append(QString("\'w20\',"));
        if(tabOutputButtonParameter_w50->isChecked())       listOutputPar.append(QString("\'w50\',"));
        if(tabOutputButtonParameter_wm50->isChecked())      listOutputPar.append(QString("\'wm50\',"));
        if(tabOutputButtonParameter_ra->isChecked())        listOutputPar.append(QString("\'ra\',"));
        if(tabOutputButtonParameter_dec->isChecked())       listOutputPar.append(QString("\'dec\',"));
        if(tabOutputButtonParameter_lon->isChecked())       listOutputPar.append(QString("\'lon\',"));
        if(tabOutputButtonParameter_lat->isChecked())       listOutputPar.append(QString("\'lat\',"));
        if(tabOutputButtonParameter_freq->isChecked())      listOutputPar.append(QString("\'freq\',"));
        if(tabOutputButtonParameter_velo->isChecked())      listOutputPar.append(QString("\'velo\',"));
        
        if(!listOutputPar.isEmpty())
        {
            listOutputPar.truncate(listOutputPar.size() - 1);
        }
        else
        {
            listOutputPar = QString("\'*\'");
        }
    }
    
    listOutputPar.prepend(QString("["));
    listOutputPar.append(QString("]"));
    
    parameters.insert(QString("writeCat.parameters"), listOutputPar);
    
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
            //int index = w->findText(parameters.value(w->objectName()));
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
    
    // Treat list of output parameters separately and set respective check boxes:
    QString listOutputPar = parameters.value(QString("writeCat.parameters"));
    
    if(listOutputPar.contains(QString("\'*\'")))
    {
        tabOutputGroupBox2->setChecked(false);
    }
    else
    {
        tabOutputGroupBox2->setChecked(true);
        
        tabOutputButtonParameter_id->setChecked(false);
        tabOutputButtonParameter_x_geo->setChecked(false);
        tabOutputButtonParameter_y_geo->setChecked(false);
        tabOutputButtonParameter_z_geo->setChecked(false);
        tabOutputButtonParameter_x->setChecked(false);
        tabOutputButtonParameter_y->setChecked(false);
        tabOutputButtonParameter_z->setChecked(false);
        tabOutputButtonParameter_x_min->setChecked(false);
        tabOutputButtonParameter_x_max->setChecked(false);
        tabOutputButtonParameter_y_min->setChecked(false);
        tabOutputButtonParameter_y_max->setChecked(false);
        tabOutputButtonParameter_z_min->setChecked(false);
        tabOutputButtonParameter_z_max->setChecked(false);
        tabOutputButtonParameter_n_pix->setChecked(false);
        tabOutputButtonParameter_snr_min->setChecked(false);
        tabOutputButtonParameter_snr_max->setChecked(false);
        tabOutputButtonParameter_snr_sum->setChecked(false);
        tabOutputButtonParameter_n_pos->setChecked(false);
        tabOutputButtonParameter_n_neg->setChecked(false);
        tabOutputButtonParameter_rel->setChecked(false);
        tabOutputButtonParameter_bf_a->setChecked(false);
        tabOutputButtonParameter_bf_b1->setChecked(false);
        tabOutputButtonParameter_bf_b2->setChecked(false);
        tabOutputButtonParameter_bf_c->setChecked(false);
        tabOutputButtonParameter_bf_chi2->setChecked(false);
        tabOutputButtonParameter_bf_flag->setChecked(false);
        tabOutputButtonParameter_bf_f_int->setChecked(false);
        tabOutputButtonParameter_bf_f_peak->setChecked(false);
        tabOutputButtonParameter_bf_w->setChecked(false);
        tabOutputButtonParameter_bf_w20->setChecked(false);
        tabOutputButtonParameter_bf_w50->setChecked(false);
        tabOutputButtonParameter_bf_xe->setChecked(false);
        tabOutputButtonParameter_bf_xp->setChecked(false);
        tabOutputButtonParameter_bf_z->setChecked(false);
        tabOutputButtonParameter_ell_maj->setChecked(false);
        tabOutputButtonParameter_ell_min->setChecked(false);
        tabOutputButtonParameter_ell_pa->setChecked(false);
        tabOutputButtonParameter_f_peak->setChecked(false);
        tabOutputButtonParameter_f_int->setChecked(false);
        tabOutputButtonParameter_f_wm50->setChecked(false);
        tabOutputButtonParameter_rms->setChecked(false);
        tabOutputButtonParameter_w20->setChecked(false);
        tabOutputButtonParameter_w50->setChecked(false);
        tabOutputButtonParameter_wm50->setChecked(false);
        tabOutputButtonParameter_ra->setChecked(false);
        tabOutputButtonParameter_dec->setChecked(false);
        tabOutputButtonParameter_lon->setChecked(false);
        tabOutputButtonParameter_lat->setChecked(false);
        tabOutputButtonParameter_freq->setChecked(false);
        tabOutputButtonParameter_velo->setChecked(false);
        
        if(listOutputPar.contains(QString("\'id\'")))        tabOutputButtonParameter_id->setChecked(true);
        if(listOutputPar.contains(QString("\'x_geo\'")))     tabOutputButtonParameter_x_geo->setChecked(true);
        if(listOutputPar.contains(QString("\'y_geo\'")))     tabOutputButtonParameter_y_geo->setChecked(true);
        if(listOutputPar.contains(QString("\'z_geo\'")))     tabOutputButtonParameter_z_geo->setChecked(true);
        if(listOutputPar.contains(QString("\'x\'")))         tabOutputButtonParameter_x->setChecked(true);
        if(listOutputPar.contains(QString("\'y\'")))         tabOutputButtonParameter_y->setChecked(true);
        if(listOutputPar.contains(QString("\'z\'")))         tabOutputButtonParameter_z->setChecked(true);
        if(listOutputPar.contains(QString("\'x_min\'")))     tabOutputButtonParameter_x_min->setChecked(true);
        if(listOutputPar.contains(QString("\'x_max\'")))     tabOutputButtonParameter_x_max->setChecked(true);
        if(listOutputPar.contains(QString("\'y_min\'")))     tabOutputButtonParameter_y_min->setChecked(true);
        if(listOutputPar.contains(QString("\'y_max\'")))     tabOutputButtonParameter_y_max->setChecked(true);
        if(listOutputPar.contains(QString("\'z_min\'")))     tabOutputButtonParameter_z_min->setChecked(true);
        if(listOutputPar.contains(QString("\'z_max\'")))     tabOutputButtonParameter_z_max->setChecked(true);
        if(listOutputPar.contains(QString("\'n_pix\'")))     tabOutputButtonParameter_n_pix->setChecked(true);
        if(listOutputPar.contains(QString("\'snr_min\'")))   tabOutputButtonParameter_snr_min->setChecked(true);
        if(listOutputPar.contains(QString("\'snr_max\'")))   tabOutputButtonParameter_snr_max->setChecked(true);
        if(listOutputPar.contains(QString("\'snr_sum\'")))   tabOutputButtonParameter_snr_sum->setChecked(true);
        if(listOutputPar.contains(QString("\'n_pos\'")))     tabOutputButtonParameter_n_pos->setChecked(true);
        if(listOutputPar.contains(QString("\'n_neg\'")))     tabOutputButtonParameter_n_neg->setChecked(true);
        if(listOutputPar.contains(QString("\'rel\'")))       tabOutputButtonParameter_rel->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_a\'")))      tabOutputButtonParameter_bf_a->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_b1\'")))     tabOutputButtonParameter_bf_b1->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_b2\'")))     tabOutputButtonParameter_bf_b2->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_c\'")))      tabOutputButtonParameter_bf_c->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_chi2\'")))   tabOutputButtonParameter_bf_chi2->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_flag\'")))   tabOutputButtonParameter_bf_flag->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_f_int\'")))  tabOutputButtonParameter_bf_f_int->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_f_peak\'"))) tabOutputButtonParameter_bf_f_peak->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_w\'")))      tabOutputButtonParameter_bf_w->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_w20\'")))    tabOutputButtonParameter_bf_w20->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_w50\'")))    tabOutputButtonParameter_bf_w50->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_xe\'")))     tabOutputButtonParameter_bf_xe->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_xp\'")))     tabOutputButtonParameter_bf_xp->setChecked(true);
        if(listOutputPar.contains(QString("\'bf_z\'")))      tabOutputButtonParameter_bf_z->setChecked(true);
        if(listOutputPar.contains(QString("\'ell_maj\'")))   tabOutputButtonParameter_ell_maj->setChecked(true);
        if(listOutputPar.contains(QString("\'ell_min\'")))   tabOutputButtonParameter_ell_min->setChecked(true);
        if(listOutputPar.contains(QString("\'ell_pa\'")))    tabOutputButtonParameter_ell_pa->setChecked(true);
        if(listOutputPar.contains(QString("\'f_peak\'")))    tabOutputButtonParameter_f_peak->setChecked(true);
        if(listOutputPar.contains(QString("\'f_int\'")))     tabOutputButtonParameter_f_int->setChecked(true);
        if(listOutputPar.contains(QString("\'f_wm50\'")))    tabOutputButtonParameter_f_wm50->setChecked(true);
        if(listOutputPar.contains(QString("\'rms\'")))       tabOutputButtonParameter_rms->setChecked(true);
        if(listOutputPar.contains(QString("\'w20\'")))       tabOutputButtonParameter_w20->setChecked(true);
        if(listOutputPar.contains(QString("\'w50\'")))       tabOutputButtonParameter_w50->setChecked(true);
        if(listOutputPar.contains(QString("\'wm50\'")))      tabOutputButtonParameter_wm50->setChecked(true);
        if(listOutputPar.contains(QString("\'ra\'")))        tabOutputButtonParameter_ra->setChecked(true);
        if(listOutputPar.contains(QString("\'dec\'")))       tabOutputButtonParameter_dec->setChecked(true);
        if(listOutputPar.contains(QString("\'lon\'")))       tabOutputButtonParameter_lon->setChecked(true);
        if(listOutputPar.contains(QString("\'lat\'")))       tabOutputButtonParameter_lat->setChecked(true);
        if(listOutputPar.contains(QString("\'freq\'")))      tabOutputButtonParameter_freq->setChecked(true);
        if(listOutputPar.contains(QString("\'velo\'")))      tabOutputButtonParameter_velo->setChecked(true);
    }
    
    return 0;
}



// ---------------------
// Slot to load settings
// ---------------------

void SoFiA::loadSettings()
{
    QString newFileName = QFileDialog::getOpenFileName(this, tr("SoFiA - Load Parameters"), QDir::currentPath());
    
    if(!newFileName.isEmpty())
    {
        setDefaults();       // Load default settings first.
        
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

int SoFiA::loadFile(QString &fileName)
{
    if(!fileName.isEmpty())
    {
        /*if(!currentFileName.isEmpty())
        {
            QMessageBox messageBox(this);
            messageBox.setWindowTitle(tr("SoFiA - Load Parameters"));
            messageBox.setText(tr("<p>Opening a new file will override all current parameter settings. Unsaved changes will be lost.</p><p>Do you wish to open a new file?</p>"));
            messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
            messageBox.setDefaultButton(QMessageBox::Ok);
            messageBox.setIcon(QMessageBox::Warning);
            int choice = messageBox.exec();
            
            if(choice != QMessageBox::Ok)
            {
                return 0;
            }
        }*/
        
        QFile inFile(fileName);
        
        if(!inFile.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            return 1;        // Error message should be generated by calling function.
        }
        
        QTextStream inStream(&inFile);
        QString keyname;
        QString value;
        
        int counter = 0;
        
        while(!inStream.atEnd())
        {
            QString line = inStream.readLine().trimmed();
            
            if((!line.isEmpty()) and (!line.startsWith("#")))
            {
                keyname = line.section(QChar('='), 0, 0).trimmed();
                value   = line.section(QChar('='), 1).trimmed();
                
                // Surprisingly, the following actually works:
                //QWidget *widget = tabs->findChild<QWidget*>(keyname);
                
                if((fileName == SOFIA_DEFAULT_SETTINGS or parameters.contains(keyname)) and keyname.size() != 0) //and widget != 0
                {
                    //parameters.insert(widget->objectName(), value);
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
        }
        else
        {
            currentFileName = fileName;
            
            setFields();
            updateFields();
            
            QString messageText = QString("");
            QString statusText = tr("Parameters loaded from %1.").arg(currentFileName.section('/', -1));
            showMessage(MESSAGE_INFO, messageText, statusText);
            
            this->setWindowTitle(tr("SoFiA - %1").arg(currentFileName.section('/', -1)));
        }
    }
    
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
        
        if(!outFile.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            QString messageText = tr("<p>Failed to write to output file %1.</p>").arg(currentFileName.section('/', -1));
            QString statusText = tr("Failed to write to output file %1.").arg(currentFileName.section('/', -1));
            showMessage(MESSAGE_ERROR, messageText, statusText);
            
            currentFileName.clear();
            this->setWindowTitle(tr("SoFiA"));
            
            return;
        }
        
        updateVariables();   // This will update all parameters before saving.
        
        QTextStream outStream(&outFile);
        
        /*QList<QWidget*> widget = tabs->findChildren<QWidget*>();
        
        foreach(QWidget *w, widget)
        {
            if(parameter.contains(w))
            {
                outStream << w->objectName() << "\t=\t" << parameter.value(w) << endl;
            }
        }*/
        
        for(QMap<QString, QString>::iterator iter = parameters.begin(); iter != parameters.end(); iter++)
        {
            outStream << iter.key() << "\t=\t" << iter.value() << endl;
        }
        
        outFile.close();
        
        QString messageText = QString("");
        QString statusText = tr("Parameters saved to %1.").arg(currentFileName.section('/', -1));
        showMessage(MESSAGE_INFO, messageText, statusText);
        
        this->setWindowTitle(tr("SoFiA - %1").arg(currentFileName.section('/', -1)));
    }
    
    return;
}



// ---------------------------
// Slot to save settings as...
// ---------------------------

void SoFiA::saveSettingsAs()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("SoFiA - Save Parameters"), QDir::currentPath());
    
    if(!fileName.isEmpty())
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
    /*QMessageBox messageBox(this);
    messageBox.setWindowTitle(tr("SoFiA - Clear Pipeline Messages"));
    messageBox.setText(tr("<p>This action will clear the pipeline message interface and discard all output messages generated by previous pipeline runs.</p><p>Do you wish to clear all pipeline messages?</p>"));
    messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
    messageBox.setDefaultButton(QMessageBox::Ok);
    messageBox.setIcon(QMessageBox::Warning);
    int choice = messageBox.exec();*/
    
    //if(choice == QMessageBox::Ok)
    if(true)
    {
        outputText->clear();
        outputProgress->setValue(0);
        
        QString messageText = tr("");
        QString statusText  = tr("Pipeline messages cleared.");
        showMessage(MESSAGE_INFO, messageText, statusText);
        
        updateActions();
    }
    
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
    
    // Disable RMS mode field in threshold finder if clip mode is 'absolute':
    tabSourceFindingFieldRmsMode2->setEnabled(tabSourceFindingFieldClipMethod->currentIndex() == 0 and tabSourceFindingGroupBox2->isChecked());
    
    // Enable/disable writing of filtered cube when no filters are selected:
    tabOutputButtonFilteredCube->setEnabled(tabInFilterGroupBox1->isChecked() or tabInFilterGroupBox2->isChecked() or tabInFilterGroupBox3->isChecked());
    
    // Enable output filter fields when respective parameter is selected:
    tabOutFilterFieldW50Min->setEnabled(tabOutFilterButtonW50->isChecked());
    tabOutFilterFieldW50Max->setEnabled(tabOutFilterButtonW50->isChecked());
    tabOutFilterFieldW20Min->setEnabled(tabOutFilterButtonW20->isChecked());
    tabOutFilterFieldW20Max->setEnabled(tabOutFilterButtonW20->isChecked());
    tabOutFilterFieldFpeakMin->setEnabled(tabOutFilterButtonFpeak->isChecked());
    tabOutFilterFieldFpeakMax->setEnabled(tabOutFilterButtonFpeak->isChecked());
    tabOutFilterFieldFintMin->setEnabled(tabOutFilterButtonFint->isChecked());
    tabOutFilterFieldFintMax->setEnabled(tabOutFilterButtonFint->isChecked());
    
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
    
    // Disable output parameter buttons if "all" button is selected:
    tabOutputButtonParameter_id->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_x_geo->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_y_geo->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_z_geo->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_x->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_y->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_z->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_x_min->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_x_max->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_y_min->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_y_max->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_z_min->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_z_max->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_n_pix->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_snr_min->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_snr_max->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_snr_sum->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_n_pos->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_n_neg->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_rel->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_a->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_b1->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_b2->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_c->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_chi2->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_flag->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_f_int->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_f_peak->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_w->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_w20->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_w50->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_xe->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_xp->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_bf_z->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_ell_maj->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_ell_min->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_ell_pa->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_f_peak->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_f_int->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_f_wm50->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_rms->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_w20->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_w50->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_wm50->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_ra->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_dec->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_lon->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_lat->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_freq->setEnabled(tabOutputGroupBox2->isChecked());
    tabOutputButtonParameter_velo->setEnabled(tabOutputGroupBox2->isChecked());
    
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
    
    if(tabParametrisationGroupBox1->isChecked()) toolBoxPA->setItemIcon(0, iconTaskComplete);
    else                                         toolBoxPA->setItemIcon(0, iconTaskReject);
    if(tabParametrisationGroupBox2->isChecked()) toolBoxPA->setItemIcon(1, iconTaskComplete);
    else                                         toolBoxPA->setItemIcon(1, iconTaskReject);
    
    if(tabOutputButtonASCII->isChecked() or tabOutputButtonXML->isChecked() or tabOutputButtonSQL->isChecked() or (tabOutputButtonFilteredCube->isEnabled() and tabOutputButtonFilteredCube->isChecked()) or tabOutputButtonMask->isChecked() or tabOutputButtonMom0->isChecked() or tabOutputButtonMom1->isChecked() or tabOutputButtonCubelets->isChecked()) toolBoxOP->setItemIcon(0, iconTaskComplete);
    else toolBoxOP->setItemIcon(0, iconTaskReject);
    
    if(tabOutputGroupBox2->isChecked())          toolBoxOP->setItemIcon(1, iconTaskComplete);
    else                                         toolBoxOP->setItemIcon(1, iconTaskReject);
    
    updateActions();
    
    return;
}



// --------------------------------------
// Function to update actions and buttons
// --------------------------------------

void SoFiA::updateActions()
{
    // Activate or de-activate actions and buttons
    
    tabOutputButtonGo->setEnabled(!((tabInputFieldData->text()).isEmpty()) and (pipelineProcess->state() == QProcess::NotRunning));
    actionRun        ->setEnabled(!((tabInputFieldData->text()).isEmpty()) and (pipelineProcess->state() == QProcess::NotRunning));
    
    actionAbort->setEnabled(pipelineProcess->state() == QProcess::Running);
    //actionExit->setEnabled(pipelineProcess->state() == QProcess::NotRunning);
    
    actionSaveLogAs->setEnabled(outputText->toPlainText() != "" and pipelineProcess->state() == QProcess::NotRunning);
    actionClearLog->setEnabled(outputText->toPlainText() != "" and pipelineProcess->state() == QProcess::NotRunning);
    actionShowCatalogue->setEnabled(!((tabInputFieldData->text()).isEmpty()));
    
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



// -------------------------------------
// Slot to select optical catalogue file
// -------------------------------------

void SoFiA::selectOpticalCatalogFile()
{
    selectFile(tabInputFieldCatalog);
    
    return;
}



// ---------------------------------
// Slot to select input weights cube
// ---------------------------------

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
    QMessageBox messageBox(this);
    messageBox.setWindowTitle(tr("SoFiA - New Parameter File"));
    messageBox.setText(tr("<p>This action will reset all parameters to their default values and close the current file. All unsaved changes will be lost.</p><p>Do you wish to create a new parameter file?</p>"));
    messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
    messageBox.setDefaultButton(QMessageBox::Ok);
    messageBox.setIcon(QMessageBox::Warning);
    int choice = messageBox.exec();
    
    if(choice == QMessageBox::Ok)
    {
        setDefaults();
    }
    
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
    QString messageText = tr("<h3>About SoFiA</h3><p>Version 0.4.0 (using Qt %1)</p><p>SoFiA, the <b>Source Finding Application</b>, is a 3D source finding pipeline designed to detect and parameterise galaxies in HI data cubes. The acronym SoFiA is based on the Greek word %2, which means wisdom.</p><p>SoFiA is free software: you can redistribute it and/or modify it under the terms of the <b>GNU General Public License</b> as published by the Free Software Foundation, either version 3 of the licence, or (at your option) any later version.</p><p>SoFiA is distributed in the hope that it will be useful, but <b>without any warranty</b>; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.</p><p>You should have received a copy of the GNU General Public License along with SoFiA. If not, see <a href=\"http://www.gnu.org/licenses/\">http://www.gnu.org/licenses/</a>.</p><p>SoFiA uses the Oxygen icon set which is licensed under version&nbsp;3 of the <a href=\"http://www.gnu.org/licenses/lgpl-3.0.txt\">GNU Lesser General Public License</a>. For more details please visit the website of the <a href=\"http://www.oxygen-icons.org/\"><s>Oxygen project</s></a> (website no longer exists).</p><p>&copy; 2013&ndash;2014 The SoFiA Authors</p>").arg(QString(qVersion())).arg(QString::fromUtf8("σοφία"));
    QString statusText = QString("");
    showMessage(MESSAGE_INFO, messageText, statusText);
    
    return;
}



// --------------------
// Slot to run pipeline
// --------------------

void SoFiA::runPipeline()
{
    /*QMessageBox messageBox(this);
    messageBox.setWindowTitle(tr("SoFiA - Run Pipeline"));
    messageBox.setText(tr("<p>This action will run the pipeline with the parameters as currently set in the user interface. Depending on your hardware, the size of your data file, and the actual parameters, this can take up to several hours.</p><p>Do you wish to run the pipeline?</p>"));
    messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
    messageBox.setDefaultButton(QMessageBox::Ok);
    messageBox.setIcon(QMessageBox::Information);
    int choice = messageBox.exec();
    
    if(choice == QMessageBox::Ok)
    {
    }*/
    
    if(pipelineProcess->state() != QProcess::NotRunning)
    {
        QString messageText = tr("<p>The pipeline is already running.</p>");
        QString statusText = tr("Pipeline already running.");
        showMessage(MESSAGE_INFO, messageText, statusText);
    }
    else
    {
        bool flagTmpFile = false;
        if(dockWidgetOutput->isHidden()) dockWidgetOutput->show();
        
        QString command("python");
        QStringList arguments;
        
        if(currentFileName.isEmpty())
        {
            currentFileName = SOFIA_TEMP_FILE;
            flagTmpFile = true;
        }
        
        QString SOFIA_FULL_PATH = getenv("SOFIA_PIPELINE_PATH");
        
        // Replace name of pipeline if optically motivated source finding is requested:
        if(tabInputGroupBox2->isChecked()) SOFIA_FULL_PATH.replace("sofia_pipeline.py", "optical_find.py");
        
        arguments << SOFIA_FULL_PATH.toUtf8().data() << currentFileName;
        saveSettings();
        
        if(flagTmpFile == true)
        {
            currentFileName.clear();
            this->setWindowTitle(tr("SoFiA"));
        }
        
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
    
    unsigned int progress = outputProgress->value();
    
    if(outputStd.contains(QString("--- SoFiA: Reading default parameters ---")))                progress = 0;
    if(outputStd.contains(QString("--- SoFiA: Reading user parameters ---")))                   progress = 5;
    if(outputStd.contains(QString("--- SoFiA: Reading data cube(s) ---")))                      progress = 10;
    if(outputStd.contains(QString("--- SoFiA: Running input filters ---")))                     progress = 15;
    if(outputStd.contains(QString("--- SoFiA: Running source finder ---")))                     progress = 25;
    if(outputStd.contains(QString("--- SoFiA: Merging detections ---")))                        progress = 35;
    if(outputStd.contains(QString("--- SoFiA: Determining reliability ---")))                   progress = 45;
    if(outputStd.contains(QString("--- SoFiA: Removing unreliable sources ---")))               progress = 50;
    if(outputStd.contains(QString("--- SoFiA: Parameterising sources ---")))                    progress = 55;
    if(outputStd.contains(QString("--- SoFiA: Correcting parameters for sub-cube offset ---"))) progress = 65;
    if(outputStd.contains(QString("--- SoFiA: Writing mask cube ---")))                         progress = 70;
    if(outputStd.contains(QString("--- SoFiA: Writing moment-0 map ---")))                      progress = 75;
    if(outputStd.contains(QString("--- SoFiA: Writing moment-1 map ---")))                      progress = 80;
    if(outputStd.contains(QString("--- SoFiA: Writing cubelets ---")))                          progress = 85;
    if(outputStd.contains(QString("--- SoFiA: Adding WCS position to catalogue ---")))          progress = 90;
    if(outputStd.contains(QString("--- SoFiA: Writing output catalogue ---")))                  progress = 95;
    if(outputStd.contains(QString("--- SoFiA: Pipeline finished ---")))                         progress = 100;
    
    outputProgress->setValue(progress);
    
    if(!outputStd.isEmpty())
    {
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
    
    if(!outputErr.isEmpty())
    {
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
    
    outputProgress->setValue(0);
    
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
            
            if(!spreadsheet->isHidden()) showCatalogue();  // Reload catalogue if currently visible.
        }
        else
        {
            // Pipeline finished with error:
            QString messageText("");
            QString statusText = tr("Pipeline failed.");
            showMessage(MESSAGE_INFO, messageText, statusText);
            
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
        
        outputText->setTextColor(Qt::red);
	outputText->insertPlainText(QString("Pipeline aborted with exit code %1.\n").arg(exitCode));
	outputText->verticalScrollBar()->setValue(outputText->verticalScrollBar()->maximum());
    }
    
    updateActions();
    
    return;
}



// ---------------------------
// Slot to cancel pipeline run
// ---------------------------

void SoFiA::pipelineProcessCancel()
{
    actionAbort->setEnabled(false);
    
    if(pipelineProcess->state() == QProcess::Running)
    {
        pipelineProcess->terminate();                      // First try to terminate process.
    }
    
    if(!pipelineProcess->waitForFinished(30000))           // Give it 30 s to terminate.
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
// Slot to show source catalogue
// -----------------------------

void SoFiA::showCatalogue()
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
    filename.append("_cat.xml");
    
    // Determine directory name:
    if(dirname.isEmpty() and separatorPos > 0) dirname = fullFilePath.left(separatorPos);
    if(not (dirname.isEmpty() or dirname.endsWith("/"))) dirname.append("/");
    
    // Concatenate directory and file names:
    filename.prepend(dirname);
    
    // Load catalogue:
    if(spreadsheet->loadCatalog(filename))
    {
        QString messageText = tr("<p>Failed to load source catalogue:</p><p>\"%1\"</p>").arg(filename);
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
    tabOutFilter       = new QWidget(tabs);
    tabOutput          = new QWidget(tabs);
    
    tabs->addTab(tabInput,           tr("Input"));
    tabs->addTab(tabInFilter,        tr("Input Filter"));
    tabs->addTab(tabSourceFinding,   tr("Source Finding"));
    tabs->addTab(tabMerging,         tr("Merging"));
    tabs->addTab(tabParametrisation, tr("Parameterisation"));
    tabs->addTab(tabOutFilter,       tr("Output Filter"));
    tabs->addTab(tabOutput,          tr("Output"));
    
    tabs->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Minimum);
    tabs->setUsesScrollButtons(false);
    
    // Note that additional spaces at the end of the display text for checkboxes 
    // and radio buttons are a workaround for rendering issues under Mac OS X 
    // whereby the last character of the text was partially cut off.
    
    // Set up input tab
    // ----------------
    
    toolBoxIP = new QToolBox(tabInput);
    
    tabInputLayout = new QVBoxLayout;
    
    // input files
    tabInputGroupBox1 = new QGroupBox(tr("Files and settings"), toolBoxIP);
    tabInputForm1 = new QFormLayout;
    
    tabInputWidgetData = new QWidget(tabInputGroupBox1);
    tabInputLayoutData = new QHBoxLayout;
    tabInputFieldData  = new QLineEdit(tabInputWidgetData);
    tabInputFieldData->setObjectName("import.inFile");
    //tabInputFieldData->setToolTip("Name of input data cube (required)");
    connect(tabInputFieldData, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    tabInputButtonData = new QPushButton(tr("Select..."), tabInputWidgetData);
    connect(tabInputButtonData, SIGNAL(clicked()), this, SLOT(selectInputDataFile()));
    tabInputButtonData->setIcon(iconDocumentOpen);
    tabInputLayoutData->addWidget(tabInputFieldData);
    tabInputLayoutData->addWidget(tabInputButtonData);
    tabInputLayoutData->setContentsMargins(0, 0, 0, 0);
    tabInputWidgetData->setLayout(tabInputLayoutData);
    
    tabInputWidgetMask = new QWidget(tabInputGroupBox1);
    tabInputLayoutMask = new QHBoxLayout;
    tabInputFieldMask  = new QLineEdit(tabInputWidgetMask);
    tabInputFieldMask->setObjectName("import.maskFile");
    //tabInputFieldMask->setToolTip("Name of mask cube (optional)");
    tabInputButtonMask = new QPushButton(tr("Select..."), tabInputWidgetMask);
    connect(tabInputButtonMask, SIGNAL(clicked()), this, SLOT(selectInputMaskFile()));
    tabInputButtonMask->setIcon(iconDocumentOpen);
    tabInputLayoutMask->addWidget(tabInputFieldMask);
    tabInputLayoutMask->addWidget(tabInputButtonMask);
    tabInputLayoutMask->setContentsMargins(0, 0, 0, 0);
    tabInputWidgetMask->setLayout(tabInputLayoutMask);
    
    tabInputWidgetWeights = new QWidget(tabInputGroupBox1);
    tabInputLayoutWeights = new QHBoxLayout;
    tabInputFieldWeights  = new QLineEdit(tabInputWidgetWeights);
    tabInputFieldWeights->setObjectName("import.weightsFile");
    //tabInputFieldWeights->setToolTip("Name of data cube containing weights (optional)");
    tabInputButtonWeights = new QPushButton(tr("Select..."), tabInputWidgetWeights);
    connect(tabInputButtonWeights, SIGNAL(clicked()), this, SLOT(selectInputWeightsFile()));
    tabInputButtonWeights->setIcon(iconDocumentOpen);
    tabInputLayoutWeights->addWidget(tabInputFieldWeights);
    tabInputLayoutWeights->addWidget(tabInputButtonWeights);
    tabInputLayoutWeights->setContentsMargins(0, 0, 0, 0);
    tabInputWidgetWeights->setLayout(tabInputLayoutWeights);
    
    tabInputFieldWeightsFunction = new QLineEdit(tabInputGroupBox1);
    tabInputFieldWeightsFunction->setObjectName("import.weightsFunction");
    //tabInputFieldWeightsFunction->setToolTip("Analytic function describing data weights (optional)");
    
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
    tabInputForm4 = new QFormLayout;
    
    tabInputFieldSubcube = new QLineEdit(tabInputGroupBox4);
    tabInputFieldSubcube->setObjectName("import.subcube");
    //tabInputFieldSubcube->setToolTip(tr("Define subcube range to be used (optional)"));
    tabInputFieldSubcubeMode = new QComboBox(tabInputGroupBox4);
    tabInputFieldSubcubeMode->setObjectName("import.subcubeMode");
    //tabInputFieldSubcubeMode->setToolTip(tr("Range given in pixels or world coordinates"));
    tabInputFieldSubcubeMode->addItem(tr("Pixels"), QVariant(QString("pixel")));
    tabInputFieldSubcubeMode->addItem(tr("World coordinates"), QVariant(QString("world")));
    
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
    tabInputForm3 = new QFormLayout;
    
    tabInputFieldFlags = new QLineEdit(tabInputGroupBox3);
    tabInputFieldFlags->setObjectName("flag.regions");
    //tabInputFieldFlags->setToolTip("Pixel/channel ranges to be flagged (optional)");
    
    tabInputForm3->addRow(tr("Range:"), tabInputFieldFlags);
    tabInputForm3->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
    tabInputGroupBox3->setLayout(tabInputForm3);
    
    // optical source finding
    tabInputGroupBox2 = new QGroupBox(tr("Enable"), toolBoxIP);
    tabInputGroupBox2->setObjectName("steps.doOptical");
    tabInputGroupBox2->setCheckable(true);
    tabInputGroupBox2->setChecked(false);
    connect(tabInputGroupBox2, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabInputForm2 = new QFormLayout;
    
    tabInputWidgetCatalog = new QWidget(tabInputGroupBox2);
    tabInputLayoutCatalog = new QHBoxLayout;
    tabInputFieldCatalog = new QLineEdit(tabInputWidgetCatalog);
    tabInputFieldCatalog->setObjectName("optical.sourceCatalogue");
    //tabInputFieldCatalog->setToolTip("Source catalogue for catalogue-based source finding");
    tabInputButtonCatalog = new QPushButton(tr("Select..."), tabInputWidgetCatalog);
    connect(tabInputButtonCatalog, SIGNAL(clicked()), this, SLOT(selectOpticalCatalogFile()));
    tabInputButtonCatalog->setIcon(iconDocumentOpen);
    tabInputLayoutCatalog->addWidget(tabInputFieldCatalog);
    tabInputLayoutCatalog->addWidget(tabInputButtonCatalog);
    tabInputLayoutCatalog->setContentsMargins(0, 0, 0, 0);
    tabInputWidgetCatalog->setLayout(tabInputLayoutCatalog);
    
    tabInputFieldSpatialSize = new QLineEdit(tabInputGroupBox2);
    tabInputFieldSpatialSize->setObjectName("optical.spatSize");
    //tabInputFieldSpatialSize->setToolTip("Spatial size of subcube to be searched");
    connect(tabInputFieldSpatialSize, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    
    tabInputFieldSpectralSize = new QLineEdit(tabInputGroupBox2);
    tabInputFieldSpectralSize->setObjectName("optical.specSize");
    //tabInputFieldSpectralSize->setToolTip("Spectral size of subcube to be searched");
    connect(tabInputFieldSpectralSize, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    
    tabInputFieldMultiCat = new QCheckBox(tr("Create separate output catalogues "), tabInputGroupBox2);
    tabInputFieldMultiCat->setObjectName("optical.storeMultiCat");
    //tabInputFieldMultiCat->setToolTip("Create separate output catalogue file for each field?");
    tabInputFieldMultiCat->setChecked(false);
    
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
    
    tabInFilterLayout = new QVBoxLayout;
    
    tabInFilterGroupBox1 = new QGroupBox(tr("Enable"), toolBoxIF);
    tabInFilterGroupBox1->setObjectName("steps.doSmooth");
    tabInFilterGroupBox1->setCheckable(true);
    tabInFilterGroupBox1->setChecked(false);
    connect(tabInFilterGroupBox1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    
    tabInFilterForm1 = new QFormLayout;
    
    tabInFilterFieldKernel = new QComboBox(tabInFilterGroupBox1);
    tabInFilterFieldKernel->setObjectName("smooth.kernel");
    tabInFilterFieldKernel->addItem(tr("Gaussian"), QVariant(QString("gaussian")));
    tabInFilterFieldKernel->addItem(tr("Boxcar"), QVariant(QString("boxcar")));
    tabInFilterFieldKernel->addItem(tr("Median"), QVariant(QString("median")));
    
    tabInFilterFieldBorder = new QComboBox(tabInFilterGroupBox1);
    tabInFilterFieldBorder->setObjectName("smooth.edgeMode");
    tabInFilterFieldBorder->addItem(tr("Constant"), QVariant(QString("constant")));
    tabInFilterFieldBorder->addItem(tr("Reflect"), QVariant(QString("reflect")));
    tabInFilterFieldBorder->addItem(tr("Mirror"), QVariant(QString("mirror")));
    tabInFilterFieldBorder->addItem(tr("Nearest"), QVariant(QString("nearest")));
    tabInFilterFieldBorder->addItem(tr("Wrap"), QVariant(QString("wrap")));
    
    tabInFilterFieldSmoothingSpatialLon  = new QLineEdit(tabInFilterGroupBox1);
    tabInFilterFieldSmoothingSpatialLon->setObjectName("smooth.kernelX");
    //tabInFilterFieldSmoothingSpatialLon->setToolTip(tr("Longitude smoothing scale in pixels"));
    tabInFilterFieldSmoothingSpatialLon->setMaxLength(10);
    tabInFilterFieldSmoothingSpatialLon->setMaximumWidth(100);
    connect(tabInFilterFieldSmoothingSpatialLon, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    tabInFilterFieldSmoothingSpatialLat  = new QLineEdit(tabInFilterGroupBox1);
    tabInFilterFieldSmoothingSpatialLat->setObjectName("smooth.kernelY");
    //tabInFilterFieldSmoothingSpatialLat->setToolTip(tr("Latitude smoothing scale in pixels"));
    tabInFilterFieldSmoothingSpatialLat->setMaxLength(10);
    tabInFilterFieldSmoothingSpatialLat->setMaximumWidth(100);
    connect(tabInFilterFieldSmoothingSpatialLat, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    tabInFilterFieldSmoothingSpectral = new QLineEdit(tabInFilterGroupBox1);
    tabInFilterFieldSmoothingSpectral->setObjectName("smooth.kernelZ");
    //tabInFilterFieldSmoothingSpectral->setToolTip(tr("Spectral smoothing scale in channels"));
    tabInFilterFieldSmoothingSpectral->setMaxLength(10);
    tabInFilterFieldSmoothingSpectral->setMaximumWidth(100);
    connect(tabInFilterFieldSmoothingSpectral, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    
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
    
    tabInFilterForm2 = new QFormLayout;
    
    tabInFilterWidgetScaleXYZ = new QWidget(tabInFilterGroupBox2);
    tabInFilterLayoutScaleXYZ = new QHBoxLayout;
    tabInFilterFieldScaleX = new QCheckBox(tr("X "), tabInFilterWidgetScaleXYZ);
    tabInFilterFieldScaleX->setObjectName("scaleNoise.scaleX");
    tabInFilterFieldScaleX->setChecked(false);
    tabInFilterFieldScaleY = new QCheckBox(tr("Y "), tabInFilterWidgetScaleXYZ);
    tabInFilterFieldScaleY->setObjectName("scaleNoise.scaleY");
    tabInFilterFieldScaleY->setChecked(false);
    tabInFilterFieldScaleZ = new QCheckBox(tr("Z "), tabInFilterWidgetScaleXYZ);
    tabInFilterFieldScaleZ->setObjectName("scaleNoise.scaleZ");
    tabInFilterFieldScaleZ->setChecked(true);
    tabInFilterLayoutScaleXYZ->setContentsMargins(0, 0, 0, 0);
    //tabInFilterLayoutScaleXYZ->setSpacing(0);
    tabInFilterLayoutScaleXYZ->addWidget(tabInFilterFieldScaleX);
    tabInFilterLayoutScaleXYZ->addWidget(tabInFilterFieldScaleY);
    tabInFilterLayoutScaleXYZ->addWidget(tabInFilterFieldScaleZ);
    tabInFilterWidgetScaleXYZ->setLayout(tabInFilterLayoutScaleXYZ);
    
    tabInFilterFieldStatistic = new QComboBox(tabInFilterGroupBox2);
    tabInFilterFieldStatistic->setObjectName("scaleNoise.statistic");
    tabInFilterFieldStatistic->addItem(tr("Gaussian fit to negative fluxes"), QVariant(QString("negative")));
    tabInFilterFieldStatistic->addItem(tr("Median absolute deviation"), QVariant(QString("mad")));
    tabInFilterFieldStatistic->addItem(tr("Standard deviation"), QVariant(QString("std")));
    
    tabInFilterFieldEdgeX  = new QSpinBox(tabInFilterGroupBox2);
    tabInFilterFieldEdgeX->setObjectName("scaleNoise.edgeX");
    tabInFilterFieldEdgeX->setMaximumWidth(100);
    tabInFilterFieldEdgeX->setMinimum(0);
    tabInFilterFieldEdgeX->setMaximum(100);
    tabInFilterFieldEdgeY = new QSpinBox(tabInFilterGroupBox2);
    tabInFilterFieldEdgeY->setObjectName("scaleNoise.edgeY");
    tabInFilterFieldEdgeY->setMaximumWidth(100);
    tabInFilterFieldEdgeY->setMinimum(0);
    tabInFilterFieldEdgeY->setMaximum(100);
    tabInFilterFieldEdgeZ = new QSpinBox(tabInFilterGroupBox2);
    tabInFilterFieldEdgeZ->setObjectName("scaleNoise.edgeZ");
    tabInFilterFieldEdgeZ->setMaximumWidth(100);
    tabInFilterFieldEdgeZ->setMinimum(0);
    tabInFilterFieldEdgeZ->setMaximum(100);
    
    tabInFilterForm2->addRow(tr("Dimensions:"), tabInFilterWidgetScaleXYZ);
    tabInFilterForm2->addRow(tr("Statistic:"), tabInFilterFieldStatistic);
    tabInFilterForm2->addRow(tr("Edge X:"), tabInFilterFieldEdgeX);
    tabInFilterForm2->addRow(tr("Edge Y:"), tabInFilterFieldEdgeY);
    tabInFilterForm2->addRow(tr("Edge Z:"), tabInFilterFieldEdgeZ);
    
    
    
    tabInFilterGroupBox3 = new QGroupBox(tr("Enable"), toolBoxIF);
    tabInFilterGroupBox3->setObjectName("steps.doWavelet");
    tabInFilterGroupBox3->setCheckable(true);
    tabInFilterGroupBox3->setChecked(false);
    connect(tabInFilterGroupBox3, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    
    tabInFilterForm3 = new QFormLayout;
    
    tabInFilterField2d1dThreshold = new QLineEdit(tabInFilterGroupBox3);
    tabInFilterField2d1dThreshold->setObjectName("wavelet.threshold");
    //tabInFilterField2d1dThreshold->setToolTip(tr("Wavelet reconstruction threshold in multiples of the rms noise"));
    tabInFilterField2d1dThreshold->setMaximumWidth(100);
    tabInFilterField2d1dThreshold->setMaxLength(10);
    connect(tabInFilterField2d1dThreshold, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    
    tabInFilterField2d1dIterations = new QSpinBox(tabInFilterGroupBox3);
    tabInFilterField2d1dIterations->setObjectName("wavelet.iterations");
    //tabInFilterField2d1dIterations->setToolTip(tr("Number of iterations in the reconstruction process"));
    tabInFilterField2d1dIterations->setMaximumWidth(100);
    tabInFilterField2d1dIterations->setMinimum(1);
    tabInFilterField2d1dIterations->setMaximum(50);
    
    tabInFilterField2d1dScaleXY = new QSpinBox(tabInFilterGroupBox3);
    tabInFilterField2d1dScaleXY->setObjectName("wavelet.scaleXY");
    //tabInFilterField2d1dScaleXY->setToolTip(tr("Number of spatial scales used in decomposition"));
    tabInFilterField2d1dScaleXY->setMaximumWidth(100);
    tabInFilterField2d1dScaleXY->setMinimum(-1);
    tabInFilterField2d1dScaleXY->setMaximum(50);
    
    tabInFilterField2d1dScaleZ = new QSpinBox(tabInFilterGroupBox3);
    tabInFilterField2d1dScaleZ->setObjectName("wavelet.scaleZ");
    //tabInFilterField2d1dScaleZ->setToolTip(tr("Number of spectral scales used in decomposition"));
    tabInFilterField2d1dScaleZ->setMaximumWidth(100);
    tabInFilterField2d1dScaleZ->setMinimum(-1);
    tabInFilterField2d1dScaleZ->setMaximum(50);
    
    tabInFilterField2d1dPositivity = new QCheckBox(tr("Enable "), tabInFilterGroupBox3);
    tabInFilterField2d1dPositivity->setObjectName("wavelet.positivity");
    //tabInFilterField2d1dPositivity->setToolTip("Only include positive wavelet components");
    tabInFilterField2d1dPositivity->setChecked(false);
    
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
    
    tabSourceFindingLayout = new QVBoxLayout;
    
    //S+C
    tabSourceFindingGroupBox1 = new QGroupBox(tr("Enable"), toolBoxSF);
    tabSourceFindingGroupBox1->setObjectName("steps.doSCfind");
    tabSourceFindingGroupBox1->setCheckable(true);
    tabSourceFindingGroupBox1->setChecked(true);
    connect(tabSourceFindingGroupBox1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    
    tabSourceFindingWidget1Left = new QWidget(tabSourceFindingGroupBox1);
    tabSourceFindingWidget1Right = new QWidget(tabSourceFindingGroupBox1);
    tabSourceFindingForm1Left = new QFormLayout;
    tabSourceFindingForm1Right = new QFormLayout;
    tabSourceFindingForm1Layout = new QHBoxLayout;
    
    tabSourceFindingFieldThreshold  = new QLineEdit(tabSourceFindingGroupBox1);
    tabSourceFindingFieldThreshold->setObjectName("SCfind.threshold");
    tabSourceFindingFieldThreshold->setMaximumWidth(100);
    tabSourceFindingFieldThreshold->setMaxLength(10);
    connect(tabSourceFindingFieldThreshold, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    
    tabSourceFindingFieldEdgeMode = new QComboBox(tabSourceFindingGroupBox1);
    tabSourceFindingFieldEdgeMode->setObjectName("SCfind.edgeMode");
    tabSourceFindingFieldEdgeMode->addItem(tr("Constant"), QVariant(QString("constant")));
    tabSourceFindingFieldEdgeMode->addItem(tr("Reflect"), QVariant(QString("reflect")));
    tabSourceFindingFieldEdgeMode->addItem(tr("Mirror"), QVariant(QString("mirror")));
    tabSourceFindingFieldEdgeMode->addItem(tr("Nearest"), QVariant(QString("nearest")));
    tabSourceFindingFieldEdgeMode->addItem(tr("Wrap"), QVariant(QString("wrap")));
    
    tabSourceFindingFieldRmsMode = new QComboBox(tabSourceFindingGroupBox1);
    tabSourceFindingFieldRmsMode->setObjectName("SCfind.rmsMode");
    tabSourceFindingFieldRmsMode->addItem(tr("Gaussian fit to negative fluxes"), QVariant(QString("negative")));
    tabSourceFindingFieldRmsMode->addItem(tr("Median absolute deviation"), QVariant(QString("mad")));
    tabSourceFindingFieldRmsMode->addItem(tr("Standard deviation"), QVariant(QString("std")));
    
    tabSourceFindingFieldKunit = new QComboBox(tabSourceFindingGroupBox1);
    tabSourceFindingFieldKunit->setObjectName("SCfind.kernelUnit");
    tabSourceFindingFieldKunit->addItem(tr("Pixels"), QVariant(QString("pixel")));
    tabSourceFindingFieldKunit->addItem(tr("World coordinates"), QVariant(QString("world")));
    
    tabSourceFindingFieldKernels = new QTextEdit(tabSourceFindingGroupBox1);
    tabSourceFindingFieldKernels->setObjectName("SCfind.kernels");
    tabSourceFindingFieldKernels->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    tabSourceFindingFieldKernels->setMaximumHeight(120);
    
    tabSourceFindingForm1Left->addRow(tr("Threshold:"), tabSourceFindingFieldThreshold);
    tabSourceFindingForm1Left->addRow(tr("Edge mode:"), tabSourceFindingFieldEdgeMode);
    tabSourceFindingForm1Left->addRow(tr("RMS mode:"), tabSourceFindingFieldRmsMode);
    tabSourceFindingForm1Left->addRow(tr("Kernel units:"), tabSourceFindingFieldKunit);
    
    tabSourceFindingForm1Right->addRow(tr("Kernels:"), tabSourceFindingFieldKernels);
    tabSourceFindingForm1Right->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
    
    tabSourceFindingWidget1Left->setLayout(tabSourceFindingForm1Left);
    tabSourceFindingWidget1Right->setLayout(tabSourceFindingForm1Right);
    
    tabSourceFindingForm1Layout->addWidget(tabSourceFindingWidget1Left);
    tabSourceFindingForm1Layout->addWidget(tabSourceFindingWidget1Right);
    
    tabSourceFindingGroupBox1->setLayout(tabSourceFindingForm1Layout);
    
    
    
    //CNHI
    tabSourceFindingGroupBox3 = new QGroupBox(tr("Enable"), toolBoxSF);
    tabSourceFindingGroupBox3->setObjectName("steps.doCNHI");
    tabSourceFindingGroupBox3->setCheckable(true);
    tabSourceFindingGroupBox3->setChecked(false);
    connect(tabSourceFindingGroupBox3, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    
    tabSourceFindingForm3 = new QFormLayout;
    
    tabSourceFindingFieldPReq  = new QLineEdit(tabSourceFindingGroupBox3);
    tabSourceFindingFieldPReq->setObjectName("CNHI.pReq");
    tabSourceFindingFieldPReq->setMaximumWidth(100);
    tabSourceFindingFieldPReq->setMaxLength(10);
    connect(tabSourceFindingFieldPReq, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    
    tabSourceFindingFieldQReq  = new QLineEdit(tabSourceFindingGroupBox3);
    tabSourceFindingFieldQReq->setObjectName("CNHI.qReq");
    tabSourceFindingFieldQReq->setMaximumWidth(100);
    tabSourceFindingFieldQReq->setMaxLength(10);
    connect(tabSourceFindingFieldQReq, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    
    tabSourceFindingFieldMinScale = new QSpinBox(tabSourceFindingGroupBox3);
    tabSourceFindingFieldMinScale->setObjectName("CNHI.minScale");
    tabSourceFindingFieldMinScale->setMaximumWidth(100);
    tabSourceFindingFieldMinScale->setMinimum(1);
    tabSourceFindingFieldMinScale->setMaximum(50);
    
    tabSourceFindingFieldMaxScale = new QSpinBox(tabSourceFindingGroupBox3);
    tabSourceFindingFieldMaxScale->setObjectName("CNHI.maxScale");
    tabSourceFindingFieldMaxScale->setMaximumWidth(100);
    tabSourceFindingFieldMaxScale->setMinimum(-1);
    tabSourceFindingFieldMaxScale->setMaximum(50);
    
    tabSourceFindingMedianTest = new QCheckBox(tr("Enable "), tabSourceFindingGroupBox3);
    tabSourceFindingMedianTest->setObjectName("CNHI.medianTest");
    
    tabSourceFindingFieldVerbose = new QComboBox(tabSourceFindingGroupBox3);
    tabSourceFindingFieldVerbose->setObjectName("CNHI.verbose");
    tabSourceFindingFieldVerbose->addItem(tr("None"), QVariant(QString("0")));
    tabSourceFindingFieldVerbose->addItem(tr("Minimal"), QVariant(QString("1")));
    tabSourceFindingFieldVerbose->addItem(tr("Maximal"), QVariant(QString("2")));
    
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
    
    tabSourceFindingForm2 = new QFormLayout;
    
    tabSourceFindingFieldThreshold2 = new QLineEdit(tabSourceFindingGroupBox2);
    tabSourceFindingFieldThreshold2->setObjectName("threshold.threshold");
    tabSourceFindingFieldThreshold2->setMaximumWidth(100);
    tabSourceFindingFieldThreshold2->setMaxLength(10);
    connect(tabSourceFindingFieldThreshold2, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    
    tabSourceFindingFieldClipMethod = new QComboBox(tabSourceFindingGroupBox2);
    tabSourceFindingFieldClipMethod->setObjectName("threshold.clipMethod");
    tabSourceFindingFieldClipMethod->addItem(tr("Relative"), QVariant(QString("relative")));
    tabSourceFindingFieldClipMethod->addItem(tr("Absolute"), QVariant(QString("absolute")));
    connect(tabSourceFindingFieldClipMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFields()));
    
    tabSourceFindingFieldRmsMode2 = new QComboBox(tabSourceFindingGroupBox2);
    tabSourceFindingFieldRmsMode2->setObjectName("threshold.rmsMode");
    tabSourceFindingFieldRmsMode2->addItem(tr("Gaussian fit to negative fluxes"), QVariant(QString("negative")));
    tabSourceFindingFieldRmsMode2->addItem(tr("Median absolute deviation"), QVariant(QString("mad")));
    tabSourceFindingFieldRmsMode2->addItem(tr("Standard deviation"), QVariant(QString("std")));
    
    tabSourceFindingForm2->addRow(tr("Threshold:"), tabSourceFindingFieldThreshold2);
    tabSourceFindingForm2->addRow(tr("Clip mode:"), tabSourceFindingFieldClipMethod);
    tabSourceFindingForm2->addRow(tr("RMS mode:"), tabSourceFindingFieldRmsMode2);
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
    
    tabMergingLayout = new QVBoxLayout;
    
    tabMergingGroupBox1 = new QGroupBox(tr("Enable"), toolBoxME);
    tabMergingGroupBox1->setObjectName("steps.doMerge");
    tabMergingGroupBox1->setCheckable(true);
    tabMergingGroupBox1->setChecked(true);
    connect(tabMergingGroupBox1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabMergingForm1 = new QFormLayout;
    
    tabMergingFieldRadiusX = new QSpinBox(tabMergingGroupBox1);
    tabMergingFieldRadiusX->setObjectName("merge.radiusX");
    tabMergingFieldRadiusX->setMaximumWidth(100);
    tabMergingFieldRadiusX->setMinimum(0);
    tabMergingFieldRadiusX->setMaximum(50);
    tabMergingFieldRadiusY = new QSpinBox(tabMergingGroupBox1);
    tabMergingFieldRadiusY->setObjectName("merge.radiusY");
    tabMergingFieldRadiusY->setMaximumWidth(100);
    tabMergingFieldRadiusY->setMinimum(0);
    tabMergingFieldRadiusY->setMaximum(50);
    tabMergingFieldRadiusZ = new QSpinBox(tabMergingGroupBox1);
    tabMergingFieldRadiusZ->setObjectName("merge.radiusZ");
    tabMergingFieldRadiusZ->setMaximumWidth(100);
    tabMergingFieldRadiusZ->setMinimum(0);
    tabMergingFieldRadiusZ->setMaximum(50);
    tabMergingFieldMinSizeX = new QSpinBox(tabMergingGroupBox1);
    tabMergingFieldMinSizeX->setObjectName("merge.minSizeX");
    tabMergingFieldMinSizeX->setMaximumWidth(100);
    tabMergingFieldMinSizeX->setMinimum(1);
    tabMergingFieldMinSizeX->setMaximum(50);
    tabMergingFieldMinSizeY = new QSpinBox(tabMergingGroupBox1);
    tabMergingFieldMinSizeY->setObjectName("merge.minSizeY");
    tabMergingFieldMinSizeY->setMaximumWidth(100);
    tabMergingFieldMinSizeY->setMinimum(1);
    tabMergingFieldMinSizeY->setMaximum(50);
    tabMergingFieldMinSizeZ = new QSpinBox(tabMergingGroupBox1);
    tabMergingFieldMinSizeZ->setObjectName("merge.minSizeZ");
    tabMergingFieldMinSizeZ->setMaximumWidth(100);
    tabMergingFieldMinSizeZ->setMinimum(1);
    tabMergingFieldMinSizeZ->setMaximum(50);
    
    tabMergingForm1->addRow(tr("Radius X:"), tabMergingFieldRadiusX);
    tabMergingForm1->addRow(tr("Radius Y:"), tabMergingFieldRadiusY);
    tabMergingForm1->addRow(tr("Radius Z:"), tabMergingFieldRadiusZ);
    tabMergingForm1->addRow(tr("Min. size X:"), tabMergingFieldMinSizeX);
    tabMergingForm1->addRow(tr("Min. size Y:"), tabMergingFieldMinSizeY);
    tabMergingForm1->addRow(tr("Min. size Z:"), tabMergingFieldMinSizeZ);
    tabMergingGroupBox1->setLayout(tabMergingForm1);
    
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
    
    tabMergingLayout->addWidget(toolBoxME);
    tabMergingLayout->addStretch();
    tabMergingLayout->addWidget(tabMergingWidgetControls);
    tabMerging->setLayout(tabMergingLayout);
    
    
    
    
    // Set up parametrisation tab
    // --------------------------
    
    toolBoxPA = new QToolBox(tabParametrisation);
    
    tabParametrisationLayout = new QVBoxLayout;
    
    tabParametrisationGroupBox1 = new QGroupBox(tr("Enable"), toolBoxPA);
    tabParametrisationGroupBox1->setObjectName("steps.doParameterise");
    tabParametrisationGroupBox1->setCheckable(true);
    tabParametrisationGroupBox1->setChecked(true);
    connect(tabParametrisationGroupBox1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabParametrisationForm1 = new QFormLayout;
    
    tabParametrisationButtonMaskOpt = new QCheckBox(tr("Optimise mask (ellipse) "), tabParametrisationGroupBox1);
    tabParametrisationButtonMaskOpt->setObjectName("parameters.optimiseMask");
    //tabParametrisationButtonMaskOpt->setToolTip("Run mask optimisation algorithm to improve flux measurement");
    tabParametrisationButtonMaskOpt->setEnabled(true);
    tabParametrisationButtonMaskOpt->setChecked(false);
    tabParametrisationButtonDilateMask = new QCheckBox(tr("Optimise mask (dilation) "), tabParametrisationGroupBox1);
    tabParametrisationButtonDilateMask->setObjectName("parameters.dilateMask");
    //tabParametrisationButtonDilateMask->setToolTip("Dilate source mask to improve flux measurement");
    tabParametrisationButtonDilateMask->setEnabled(true);
    tabParametrisationButtonDilateMask->setChecked(false);
    tabParametrisationButtonBusyFunction = new QCheckBox(tr("Fit Busy Function "), tabParametrisationGroupBox1);
    tabParametrisationButtonBusyFunction->setObjectName("parameters.fitBusyFunction");
    //tabParametrisationButtonBusyFunction->setToolTip("Parametrise integrated spectrum by fitting Busy Function");
    tabParametrisationButtonBusyFunction->setEnabled(true);
    tabParametrisationButtonBusyFunction->setChecked(false);
    
    tabParametrisationForm1->addRow(tr(""), tabParametrisationButtonMaskOpt);
    tabParametrisationForm1->addRow(tr(""), tabParametrisationButtonDilateMask);
    tabParametrisationForm1->addRow(tr(""), tabParametrisationButtonBusyFunction);
    tabParametrisationGroupBox1->setLayout(tabParametrisationForm1);
    
    tabParametrisationGroupBox2 = new QGroupBox(tr("Enable"), toolBoxPA);
    tabParametrisationGroupBox2->setObjectName("steps.doReliability");
    tabParametrisationGroupBox2->setCheckable(true);
    tabParametrisationGroupBox2->setChecked(true);
    connect(tabParametrisationGroupBox2, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabParametrisationForm2 = new QFormLayout;
    
    tabParametrisationFieldRelMin = new QLineEdit(tabParametrisationGroupBox2);
    tabParametrisationFieldRelMin->setObjectName("reliability.threshold");
    tabParametrisationFieldRelMin->setMaximumWidth(100);
    tabParametrisationFieldRelMin->setMaxLength(10);
    connect(tabParametrisationFieldRelMin, SIGNAL(editingFinished()), this, SLOT(updateFields()));
    //tabParametrisationFieldRelMax = new QLineEdit(tabParametrisationGroupBox2);
    //tabParametrisationFieldRelMax->setObjectName("reliability.thresholdMax");
    //tabParametrisationFieldRelMax->setMaximumWidth(100);
    //tabParametrisationFieldRelMax->setMaxLength(10);
    //tabParametrisationFieldRelMax->setEnabled(false);
    //tabParametrisationFieldRelMax->setText("1.0");
    //QLabel *labelRel   = new QLabel(QString::fromUtf8("–"), tabParametrisationGroupBox2);
    
    //tabParametrisationWidgetRel = new QWidget(tabParametrisationGroupBox2);
    //tabParametrisationWidgetRel->setToolTip("Reliability cutoff for output catalogue");
    //tabParametrisationLayoutRel = new QHBoxLayout;
    //tabParametrisationLayoutRel->addWidget(tabParametrisationFieldRelMin);
    //tabParametrisationLayoutRel->addWidget(labelRel);
    //tabParametrisationLayoutRel->addWidget(tabParametrisationFieldRelMax);
    //tabParametrisationLayoutRel->addStretch();
    //tabParametrisationLayoutRel->setContentsMargins(0, 0, 0, 0);
    //tabParametrisationWidgetRel->setLayout(tabParametrisationLayoutRel);
    
    tabParametrisationFieldRelKernel = new QLineEdit(tabParametrisationGroupBox2);
    tabParametrisationFieldRelKernel->setObjectName("reliability.kernel");
    
    tabParametrisationFieldRelPlot = new QCheckBox(tr("Enable "), tabParametrisationGroupBox1);
    tabParametrisationFieldRelPlot->setObjectName("reliability.makePlot");
    //tabParametrisationFieldRelPlot->setToolTip("PDF file showing distribution of positive/negative sources in parameter space");
    tabParametrisationFieldRelPlot->setEnabled(true);
    tabParametrisationFieldRelPlot->setChecked(false);
    
    tabParametrisationForm2->addRow(tr("Threshold:"), tabParametrisationFieldRelMin);
    tabParametrisationForm2->addRow(tr("Kernel:"), tabParametrisationFieldRelKernel);
    tabParametrisationForm2->addRow(tr("Diagnostic plot:"), tabParametrisationFieldRelPlot);
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
    
    
    
    // Set up output filter tab
    // ------------------------
    
    toolBoxOF = new QToolBox(tabOutFilter);
    
    tabOutFilterLayout = new QVBoxLayout;
    
    tabOutFilterGroupBox1 = new QGroupBox(tr("Parameter range"), toolBoxOF);
    tabOutFilterGroupBox1->setEnabled(false);
    
    tabOutFilterForm1 = new QFormLayout;
    
    tabOutFilterFieldW50Min   = new QLineEdit(tabOutFilterGroupBox1);
    tabOutFilterFieldW50Min->setObjectName("widthW50Min");
    tabOutFilterFieldW50Max   = new QLineEdit(tabOutFilterGroupBox1);
    tabOutFilterFieldW50Max->setObjectName("widthW50Max");
    tabOutFilterFieldW20Min   = new QLineEdit(tabOutFilterGroupBox1);
    tabOutFilterFieldW20Min->setObjectName("widthW20Min");
    tabOutFilterFieldW20Max   = new QLineEdit(tabOutFilterGroupBox1);
    tabOutFilterFieldW20Max->setObjectName("widthW20Max");
    tabOutFilterFieldFpeakMin = new QLineEdit(tabOutFilterGroupBox1);
    tabOutFilterFieldFpeakMin->setObjectName("fpeakMin");
    tabOutFilterFieldFpeakMax = new QLineEdit(tabOutFilterGroupBox1);
    tabOutFilterFieldFpeakMax->setObjectName("fpeakMax");
    tabOutFilterFieldFintMin  = new QLineEdit(tabOutFilterGroupBox1);
    tabOutFilterFieldFintMin->setObjectName("fintegrMin");
    tabOutFilterFieldFintMax  = new QLineEdit(tabOutFilterGroupBox1);
    tabOutFilterFieldFintMax->setObjectName("fintegrMax");
    
    tabOutFilterFieldW50Min->setMaximumWidth(100);
    tabOutFilterFieldW50Max->setMaximumWidth(100);
    tabOutFilterFieldW20Min->setMaximumWidth(100);
    tabOutFilterFieldW20Max->setMaximumWidth(100);
    tabOutFilterFieldFpeakMin->setMaximumWidth(100);
    tabOutFilterFieldFpeakMax->setMaximumWidth(100);
    tabOutFilterFieldFintMin->setMaximumWidth(100);
    tabOutFilterFieldFintMax->setMaximumWidth(100);
    
    tabOutFilterFieldW50Min->setMaxLength(10);
    tabOutFilterFieldW50Max->setMaxLength(10);
    tabOutFilterFieldW20Min->setMaxLength(10);
    tabOutFilterFieldW20Max->setMaxLength(10);
    tabOutFilterFieldFpeakMin->setMaxLength(10);
    tabOutFilterFieldFpeakMax->setMaxLength(10);
    tabOutFilterFieldFintMin->setMaxLength(10);
    tabOutFilterFieldFintMax->setMaxLength(10);
    
    QLabel *labelW50   = new QLabel(QString::fromUtf8("–"), tabOutFilterGroupBox1);
    QLabel *labelW20   = new QLabel(QString::fromUtf8("–"), tabOutFilterGroupBox1);
    QLabel *labelFpeak = new QLabel(QString::fromUtf8("–"), tabOutFilterGroupBox1);
    QLabel *labelFint  = new QLabel(QString::fromUtf8("–"), tabOutFilterGroupBox1);
    
    tabOutFilterButtonW50   = new QCheckBox(tr("Apply "), tabOutFilterGroupBox1);
    tabOutFilterButtonW50->setObjectName("applyW50Filt");
    tabOutFilterButtonW50->setEnabled(true);
    connect(tabOutFilterButtonW50, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutFilterButtonW20   = new QCheckBox(tr("Apply "), tabOutFilterGroupBox1);
    tabOutFilterButtonW20->setObjectName("applyW20Filt");
    tabOutFilterButtonW20->setEnabled(true);
    connect(tabOutFilterButtonW20, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutFilterButtonFpeak = new QCheckBox(tr("Apply "), tabOutFilterGroupBox1);
    tabOutFilterButtonFpeak->setObjectName("applyFpeakFilt");
    tabOutFilterButtonFpeak->setEnabled(true);
    connect(tabOutFilterButtonFpeak, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutFilterButtonFint  = new QCheckBox(tr("Apply "), tabOutFilterGroupBox1);
    tabOutFilterButtonFint->setObjectName("applyFintFilt");
    tabOutFilterButtonFint->setEnabled(true);
    connect(tabOutFilterButtonFint, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    
    tabOutFilterWidgetW50 = new QWidget(tabOutFilterGroupBox1);
    tabOutFilterLayoutW50 = new QHBoxLayout;
    tabOutFilterLayoutW50->addWidget(tabOutFilterFieldW50Min);
    tabOutFilterLayoutW50->addWidget(labelW50);
    tabOutFilterLayoutW50->addWidget(tabOutFilterFieldW50Max);
    tabOutFilterLayoutW50->addWidget(tabOutFilterButtonW50);
    tabOutFilterLayoutW50->addStretch();
    tabOutFilterLayoutW50->setContentsMargins(0, 0, 0, 0);
    tabOutFilterWidgetW50->setLayout(tabOutFilterLayoutW50);
    
    tabOutFilterWidgetW20 = new QWidget(tabOutFilterGroupBox1);
    tabOutFilterLayoutW20 = new QHBoxLayout;
    tabOutFilterLayoutW20->addWidget(tabOutFilterFieldW20Min);
    tabOutFilterLayoutW20->addWidget(labelW20);
    tabOutFilterLayoutW20->addWidget(tabOutFilterFieldW20Max);
    tabOutFilterLayoutW20->addWidget(tabOutFilterButtonW20);
    tabOutFilterLayoutW20->addStretch();
    tabOutFilterLayoutW20->setContentsMargins(0, 0, 0, 0);
    tabOutFilterWidgetW20->setLayout(tabOutFilterLayoutW20);
    
    tabOutFilterWidgetFpeak = new QWidget(tabOutFilterGroupBox1);
    tabOutFilterLayoutFpeak = new QHBoxLayout;
    tabOutFilterLayoutFpeak->addWidget(tabOutFilterFieldFpeakMin);
    tabOutFilterLayoutFpeak->addWidget(labelFpeak);
    tabOutFilterLayoutFpeak->addWidget(tabOutFilterFieldFpeakMax);
    tabOutFilterLayoutFpeak->addWidget(tabOutFilterButtonFpeak);
    tabOutFilterLayoutFpeak->addStretch();
    tabOutFilterLayoutFpeak->setContentsMargins(0, 0, 0, 0);
    tabOutFilterWidgetFpeak->setLayout(tabOutFilterLayoutFpeak);
    
    tabOutFilterWidgetFint = new QWidget(tabOutFilterGroupBox1);
    tabOutFilterLayoutFint = new QHBoxLayout;
    tabOutFilterLayoutFint->addWidget(tabOutFilterFieldFintMin);
    tabOutFilterLayoutFint->addWidget(labelFint);
    tabOutFilterLayoutFint->addWidget(tabOutFilterFieldFintMax);
    tabOutFilterLayoutFint->addWidget(tabOutFilterButtonFint);
    tabOutFilterLayoutFint->addStretch();
    tabOutFilterLayoutFint->setContentsMargins(0, 0, 0, 0);
    tabOutFilterWidgetFint->setLayout(tabOutFilterLayoutFint);
    
    tabOutFilterForm1->addRow(tr("Line width (w50):"), tabOutFilterWidgetW50);
    tabOutFilterForm1->addRow(tr("Line width (w20):"), tabOutFilterWidgetW20);
    tabOutFilterForm1->addRow(tr("Peak flux:"), tabOutFilterWidgetFpeak);
    tabOutFilterForm1->addRow(tr("Integrated flux:"), tabOutFilterWidgetFint);
    tabOutFilterGroupBox1->setLayout(tabOutFilterForm1);
    
    tabOutFilterButtonPrev = new QPushButton(tr("Previous"), tabOutFilter);
    tabOutFilterButtonPrev->setIcon(iconGoPreviousView);
    connect(tabOutFilterButtonPrev, SIGNAL(clicked()), this, SLOT(displayPrevTab()));
    tabOutFilterButtonNext = new QPushButton(tr("Next"), tabOutFilter);
    tabOutFilterButtonNext->setIcon(iconGoNextView);
    connect(tabOutFilterButtonNext, SIGNAL(clicked()), this, SLOT(displayNextTab()));
    tabOutFilterLayoutControls = new QHBoxLayout();
    tabOutFilterLayoutControls->setContentsMargins(0, 0, 0, 0);
    tabOutFilterLayoutControls->setSpacing(0);
    tabOutFilterLayoutControls->addWidget(tabOutFilterButtonPrev);
    tabOutFilterLayoutControls->addStretch();
    tabOutFilterLayoutControls->addWidget(tabOutFilterButtonNext);
    tabOutFilterWidgetControls = new QWidget(tabOutFilter);
    tabOutFilterWidgetControls->setLayout(tabOutFilterLayoutControls);
    
    toolBoxOF->addItem(tabOutFilterGroupBox1, iconTaskReject,   tr("Output Parameter Filtering (not yet available)"));
    
    tabOutFilterLayout->addWidget(toolBoxOF);
    tabOutFilterLayout->addStretch();
    tabOutFilterLayout->addWidget(tabOutFilterWidgetControls);
    tabOutFilter->setLayout(tabOutFilterLayout);
    
    // Set up output tab
    // -----------------
    
    toolBoxOP = new QToolBox(tabOutput);
    
    tabOutputLayout = new QVBoxLayout;
    
    tabOutputGroupBox1 = new QGroupBox(tr("Files and settings"), toolBoxOP);
    
    tabOutputForm1 = new QFormLayout;
    
    tabOutputFieldBaseName = new QLineEdit(tabOutputGroupBox1);
    tabOutputFieldBaseName->setObjectName("writeCat.basename");
    //tabOutputFieldBaseName->setToolTip("Base name to be used for all output files (optional). Defaults to input file name.");
    tabOutputFieldBaseName->setEnabled(true);
    
    tabOutputWidgetDirectory = new QWidget(tabOutputGroupBox1);
    tabOutputLayoutDirectory = new QHBoxLayout;
    tabOutputFieldDirectory  = new QLineEdit(tabOutputWidgetDirectory);
    tabOutputFieldDirectory->setObjectName("writeCat.outputDir");
    //tabOutputFieldDirectory->setToolTip("Path to output directory (optional). Defaults to input cube directory.");
    tabOutputButtonDirectory = new QPushButton(tr("Select..."), tabOutputWidgetDirectory);
    connect(tabOutputButtonDirectory, SIGNAL(clicked()), this, SLOT(selectOutputDirectory()));
    tabOutputButtonDirectory->setIcon(iconDocumentOpen);
    tabOutputLayoutDirectory->addWidget(tabOutputFieldDirectory);
    tabOutputLayoutDirectory->addWidget(tabOutputButtonDirectory);
    tabOutputLayoutDirectory->setContentsMargins(0, 0, 0, 0);
    tabOutputWidgetDirectory->setLayout(tabOutputLayoutDirectory);
    
    tabOutputButtonASCII = new QCheckBox(tr("ASCII "), tabOutputGroupBox1);
    tabOutputButtonASCII->setObjectName("writeCat.writeASCII");
    //tabOutputButtonASCII->setToolTip(tr("Human-readable ASCII file"));
    tabOutputButtonASCII->setChecked(true);
    tabOutputButtonASCII->setEnabled(true);
    connect(tabOutputButtonASCII, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutputButtonXML = new QCheckBox(tr("VO table "), tabOutputGroupBox1);
    tabOutputButtonXML->setObjectName("writeCat.writeXML");
    tabOutputButtonXML->setChecked(false);
    tabOutputButtonXML->setEnabled(true);
    //tabOutputButtonXML->setToolTip(tr("Virtual Observatory XML table"));
    connect(tabOutputButtonXML, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutputButtonSQL = new QCheckBox(tr("SQL "), tabOutputGroupBox1);
    tabOutputButtonSQL->setObjectName("writeCat.writeSQL");
    tabOutputButtonSQL->setChecked(false);
    tabOutputButtonSQL->setEnabled(false);
    //tabOutputButtonSQL->setToolTip(tr("Structured Query Language"));
    connect(tabOutputButtonSQL, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    
    tabOutputWidgetFormat = new QWidget(tabOutputGroupBox1);
    tabOutputLayoutFormat = new QHBoxLayout();
    tabOutputLayoutFormat->setContentsMargins(0, 0, 0, 0);
    tabOutputLayoutFormat->setSpacing(10);
    tabOutputLayoutFormat->addWidget(tabOutputButtonASCII);
    tabOutputLayoutFormat->addWidget(tabOutputButtonXML);
    tabOutputLayoutFormat->addWidget(tabOutputButtonSQL);
    tabOutputLayoutFormat->addStretch();
    tabOutputWidgetFormat->setLayout(tabOutputLayoutFormat);
    
    tabOutputButtonFilteredCube = new QCheckBox(tr("Filtered cube "), tabOutputGroupBox1);
    tabOutputButtonFilteredCube->setObjectName("steps.doWriteFilteredCube");
    //tabOutputButtonFilteredCube->setToolTip(tr("Data cube with input filters applied"));
    tabOutputButtonFilteredCube->setChecked(false);
    connect(tabOutputButtonFilteredCube, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutputButtonMask = new QCheckBox(tr("Mask "), tabOutputGroupBox1);
    tabOutputButtonMask->setObjectName("steps.doWriteMask");
    //tabOutputButtonMask->setToolTip(tr("Source mask cube"));
    tabOutputButtonMask->setChecked(false);
    connect(tabOutputButtonMask, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutputButtonMom0 = new QCheckBox(tr("Mom. 0 "), tabOutputGroupBox1);
    tabOutputButtonMom0->setObjectName("steps.doMom0");
    tabOutputButtonMom0->setChecked(false);
    //tabOutputButtonMom0->setToolTip(tr("Integrated flux map (moment 0)"));
    connect(tabOutputButtonMom0, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutputButtonMom1 = new QCheckBox(tr("Mom. 1 "), tabOutputGroupBox1);
    tabOutputButtonMom1->setObjectName("steps.doMom1");
    tabOutputButtonMom1->setChecked(false);
    //tabOutputButtonMom1->setToolTip(tr("Velocity field map (moment 1)"));
    connect(tabOutputButtonMom1, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    tabOutputButtonCubelets = new QCheckBox(tr("Cubelets "), tabOutputGroupBox1);
    tabOutputButtonCubelets->setObjectName("steps.doCubelets");
    tabOutputButtonCubelets->setChecked(false);
    //tabOutputButtonCubelets->setToolTip(tr("Individual cubelets, moment maps and spectra for each source"));
    connect(tabOutputButtonCubelets, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    
    tabOutputButtonCompress = new QCheckBox(tr("Gzip "), tabOutputGroupBox1);
    tabOutputButtonCompress->setObjectName("writeCat.compress");
    tabOutputButtonCompress->setChecked(false);
    //tabOutputButtonCompress->setToolTip(tr("Use gzip to compress all output files"));
    
    tabOutputWidgetProducts = new QWidget(tabOutputGroupBox1);
    tabOutputLayoutProducts = new QHBoxLayout();
    tabOutputLayoutProducts->setContentsMargins(0, 0, 0, 0);
    tabOutputLayoutProducts->setSpacing(10);
    tabOutputLayoutProducts->addWidget(tabOutputButtonFilteredCube);
    tabOutputLayoutProducts->addWidget(tabOutputButtonMask);
    tabOutputLayoutProducts->addWidget(tabOutputButtonMom0);
    tabOutputLayoutProducts->addWidget(tabOutputButtonMom1);
    tabOutputLayoutProducts->addWidget(tabOutputButtonCubelets);
    tabOutputLayoutProducts->addStretch();
    tabOutputWidgetProducts->setLayout(tabOutputLayoutProducts);
    
    tabOutputForm1->addRow(tr("Base name:"), tabOutputFieldBaseName);
    tabOutputForm1->addRow(tr("Output directory:"), tabOutputWidgetDirectory);
    tabOutputForm1->addRow(tr("Source catalogue:"), tabOutputWidgetFormat);
    tabOutputForm1->addRow(tr("Data products:"), tabOutputWidgetProducts);
    tabOutputForm1->addRow(tr("Compression:"), tabOutputButtonCompress);
    tabOutputForm1->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
    tabOutputGroupBox1->setLayout(tabOutputForm1);
    
    tabOutputGroupBox2 = new QGroupBox(tr("Enable"), toolBoxOP);
    tabOutputGroupBox2->setEnabled(true);
    tabOutputGroupBox2->setCheckable(true);
    tabOutputGroupBox2->setChecked(false);
    connect(tabOutputGroupBox2, SIGNAL(toggled(bool)), this, SLOT(updateFields()));
    
    tabOutputForm2 = new QFormLayout;
    
    tabOutputButtonParameter_id        = new QCheckBox(tr("id"), tabOutputGroupBox2);
    tabOutputButtonParameter_x_geo     = new QCheckBox(tr("x_geo"), tabOutputGroupBox2);
    tabOutputButtonParameter_y_geo     = new QCheckBox(tr("y_geo"), tabOutputGroupBox2);
    tabOutputButtonParameter_z_geo     = new QCheckBox(tr("z_geo"), tabOutputGroupBox2);
    tabOutputButtonParameter_x         = new QCheckBox(tr("x"), tabOutputGroupBox2);
    tabOutputButtonParameter_y         = new QCheckBox(tr("y"), tabOutputGroupBox2);
    tabOutputButtonParameter_z         = new QCheckBox(tr("z"), tabOutputGroupBox2);
    tabOutputButtonParameter_x_min     = new QCheckBox(tr("x_min"), tabOutputGroupBox2);
    tabOutputButtonParameter_x_max     = new QCheckBox(tr("x_max"), tabOutputGroupBox2);
    tabOutputButtonParameter_y_min     = new QCheckBox(tr("y_min"), tabOutputGroupBox2);
    tabOutputButtonParameter_y_max     = new QCheckBox(tr("y_max"), tabOutputGroupBox2);
    tabOutputButtonParameter_z_min     = new QCheckBox(tr("z_min"), tabOutputGroupBox2);
    tabOutputButtonParameter_z_max     = new QCheckBox(tr("z_max"), tabOutputGroupBox2);
    tabOutputButtonParameter_n_pix     = new QCheckBox(tr("n_pix"), tabOutputGroupBox2);
    tabOutputButtonParameter_snr_min   = new QCheckBox(tr("snr_min"), tabOutputGroupBox2);
    tabOutputButtonParameter_snr_max   = new QCheckBox(tr("snr_max"), tabOutputGroupBox2);
    tabOutputButtonParameter_snr_sum   = new QCheckBox(tr("snr_sum"), tabOutputGroupBox2);
    tabOutputButtonParameter_n_pos     = new QCheckBox(tr("n_pos"), tabOutputGroupBox2);
    tabOutputButtonParameter_n_neg     = new QCheckBox(tr("n_neg"), tabOutputGroupBox2);
    tabOutputButtonParameter_rel       = new QCheckBox(tr("rel"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_a      = new QCheckBox(tr("bf_a"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_b1     = new QCheckBox(tr("bf_b1"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_b2     = new QCheckBox(tr("bf_b2"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_c      = new QCheckBox(tr("bf_c"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_chi2   = new QCheckBox(tr("bf_chi2"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_flag   = new QCheckBox(tr("bf_flag"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_f_int  = new QCheckBox(tr("bf_f_int"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_f_peak = new QCheckBox(tr("bf_f_peak"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_w      = new QCheckBox(tr("bf_w"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_w20    = new QCheckBox(tr("bf_w20"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_w50    = new QCheckBox(tr("bf_w50"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_xe     = new QCheckBox(tr("bf_xe"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_xp     = new QCheckBox(tr("bf_xp"), tabOutputGroupBox2);
    tabOutputButtonParameter_bf_z      = new QCheckBox(tr("bf_z"), tabOutputGroupBox2);
    tabOutputButtonParameter_ell_maj   = new QCheckBox(tr("ell_maj"), tabOutputGroupBox2);
    tabOutputButtonParameter_ell_min   = new QCheckBox(tr("ell_min"), tabOutputGroupBox2);
    tabOutputButtonParameter_ell_pa    = new QCheckBox(tr("ell_pa"), tabOutputGroupBox2);
    tabOutputButtonParameter_f_peak    = new QCheckBox(tr("f_peak"), tabOutputGroupBox2);
    tabOutputButtonParameter_f_int     = new QCheckBox(tr("f_int"), tabOutputGroupBox2);
    tabOutputButtonParameter_f_wm50    = new QCheckBox(tr("f_wm50"), tabOutputGroupBox2);
    tabOutputButtonParameter_rms       = new QCheckBox(tr("rms"), tabOutputGroupBox2);
    tabOutputButtonParameter_w20       = new QCheckBox(tr("w20"), tabOutputGroupBox2);
    tabOutputButtonParameter_w50       = new QCheckBox(tr("w50"), tabOutputGroupBox2);
    tabOutputButtonParameter_wm50      = new QCheckBox(tr("wm50"), tabOutputGroupBox2);
    tabOutputButtonParameter_ra        = new QCheckBox(tr("ra"), tabOutputGroupBox2);
    tabOutputButtonParameter_dec       = new QCheckBox(tr("dec"), tabOutputGroupBox2);
    tabOutputButtonParameter_lon       = new QCheckBox(tr("lon"), tabOutputGroupBox2);
    tabOutputButtonParameter_lat       = new QCheckBox(tr("lat"), tabOutputGroupBox2);
    tabOutputButtonParameter_freq      = new QCheckBox(tr("freq"), tabOutputGroupBox2);
    tabOutputButtonParameter_velo      = new QCheckBox(tr("velo"), tabOutputGroupBox2);
    
    tabOutputButtonParameter_id        -> setObjectName("parameter_id");
    tabOutputButtonParameter_x_geo     -> setObjectName("parameter_x_geo");
    tabOutputButtonParameter_y_geo     -> setObjectName("parameter_y_geo");
    tabOutputButtonParameter_z_geo     -> setObjectName("parameter_z_geo");
    tabOutputButtonParameter_x         -> setObjectName("parameter_x");
    tabOutputButtonParameter_y         -> setObjectName("parameter_y");
    tabOutputButtonParameter_z         -> setObjectName("parameter_z");
    tabOutputButtonParameter_x_min     -> setObjectName("parameter_x_min");
    tabOutputButtonParameter_x_max     -> setObjectName("parameter_x_max");
    tabOutputButtonParameter_y_min     -> setObjectName("parameter_y_min");
    tabOutputButtonParameter_y_max     -> setObjectName("parameter_y_max");
    tabOutputButtonParameter_z_min     -> setObjectName("parameter_z_min");
    tabOutputButtonParameter_z_max     -> setObjectName("parameter_z_max");
    tabOutputButtonParameter_n_pix     -> setObjectName("parameter_n_pix");
    tabOutputButtonParameter_snr_min   -> setObjectName("parameter_snr_min");
    tabOutputButtonParameter_snr_max   -> setObjectName("parameter_snr_max");
    tabOutputButtonParameter_snr_sum   -> setObjectName("parameter_snr_sum");
    tabOutputButtonParameter_n_pos     -> setObjectName("parameter_n_pos");
    tabOutputButtonParameter_n_neg     -> setObjectName("parameter_n_neg");
    tabOutputButtonParameter_rel       -> setObjectName("parameter_rel");
    tabOutputButtonParameter_bf_a      -> setObjectName("parameter_bf_a");
    tabOutputButtonParameter_bf_b1     -> setObjectName("parameter_bf_b1");
    tabOutputButtonParameter_bf_b2     -> setObjectName("parameter_bf_b2");
    tabOutputButtonParameter_bf_c      -> setObjectName("parameter_bf_c");
    tabOutputButtonParameter_bf_chi2   -> setObjectName("parameter_bf_chi2");
    tabOutputButtonParameter_bf_flag   -> setObjectName("parameter_bf_flag");
    tabOutputButtonParameter_bf_f_int  -> setObjectName("parameter_bf_f_int");
    tabOutputButtonParameter_bf_f_peak -> setObjectName("parameter_bf_f_peak");
    tabOutputButtonParameter_bf_w      -> setObjectName("parameter_bf_w");
    tabOutputButtonParameter_bf_w20    -> setObjectName("parameter_bf_w20");
    tabOutputButtonParameter_bf_w50    -> setObjectName("parameter_bf_w50");
    tabOutputButtonParameter_bf_xe     -> setObjectName("parameter_bf_xe");
    tabOutputButtonParameter_bf_xp     -> setObjectName("parameter_bf_xp");
    tabOutputButtonParameter_bf_z      -> setObjectName("parameter_bf_z");
    tabOutputButtonParameter_ell_maj   -> setObjectName("parameter_ell_maj");
    tabOutputButtonParameter_ell_min   -> setObjectName("parameter_ell_min");
    tabOutputButtonParameter_ell_pa    -> setObjectName("parameter_ell_pa");
    tabOutputButtonParameter_f_peak    -> setObjectName("parameter_f_peak");
    tabOutputButtonParameter_f_int     -> setObjectName("parameter_f_int");
    tabOutputButtonParameter_f_wm50    -> setObjectName("parameter_f_wm50");
    tabOutputButtonParameter_rms       -> setObjectName("parameter_rms");
    tabOutputButtonParameter_w20       -> setObjectName("parameter_w20");
    tabOutputButtonParameter_w50       -> setObjectName("parameter_w50");
    tabOutputButtonParameter_wm50      -> setObjectName("parameter_wm50");
    tabOutputButtonParameter_ra        -> setObjectName("parameter_ra");
    tabOutputButtonParameter_dec       -> setObjectName("parameter_dec");
    tabOutputButtonParameter_lon       -> setObjectName("parameter_lon");
    tabOutputButtonParameter_lat       -> setObjectName("parameter_lat");
    tabOutputButtonParameter_freq      -> setObjectName("parameter_freq");
    tabOutputButtonParameter_velo      -> setObjectName("parameter_velo");
    
    tabOutputButtonParameter_id        -> setToolTip(tr("Unique identifier"));
    tabOutputButtonParameter_x_geo     -> setToolTip(tr("X-coordinate of geometric centre"));
    tabOutputButtonParameter_y_geo     -> setToolTip(tr("Y-coordinate of geometric centre"));
    tabOutputButtonParameter_z_geo     -> setToolTip(tr("Z-coordinate of geometric centre"));
    tabOutputButtonParameter_x         -> setToolTip(tr("X-coordinate of flux-weighted centre"));
    tabOutputButtonParameter_y         -> setToolTip(tr("Y-coordinate of flux-weighted centre"));
    tabOutputButtonParameter_z         -> setToolTip(tr("Z-coordinate of flux-weighted centre"));
    tabOutputButtonParameter_x_min     -> setToolTip(tr("X-coordinate of lower boundary"));
    tabOutputButtonParameter_x_max     -> setToolTip(tr("X-coordinate of upper boundary"));
    tabOutputButtonParameter_y_min     -> setToolTip(tr("Y-coordinate of lower boundary"));
    tabOutputButtonParameter_y_max     -> setToolTip(tr("Y-coordinate of upper boundary"));
    tabOutputButtonParameter_z_min     -> setToolTip(tr("Z-coordinate of lower boundary"));
    tabOutputButtonParameter_z_max     -> setToolTip(tr("Z-coordinate of upper boundary"));
    tabOutputButtonParameter_n_pix     -> setToolTip(tr("Total number of spatial and spectral pixels"));
    tabOutputButtonParameter_snr_min   -> setToolTip(tr("Signal-to-noise ratio of negative peak"));
    tabOutputButtonParameter_snr_max   -> setToolTip(tr("Signal-to-noise ratio of positive peak"));
    tabOutputButtonParameter_snr_sum   -> setToolTip(tr("Integrated signal-to-noise ratio"));
    tabOutputButtonParameter_n_pos     -> setToolTip(tr("Number of spatial and spectral pixels with positive flux"));
    tabOutputButtonParameter_n_neg     -> setToolTip(tr("Number of spatial and spectral pixels with negative flux"));
    tabOutputButtonParameter_rel       -> setToolTip(tr("Reliability"));
    tabOutputButtonParameter_bf_a      -> setToolTip(tr("Busy Function free parameter A"));
    tabOutputButtonParameter_bf_b1     -> setToolTip(tr("Busy Function free parameter B<sub>1</sub>"));
    tabOutputButtonParameter_bf_b2     -> setToolTip(tr("Busy Function free parameter B<sub>2</sub>"));
    tabOutputButtonParameter_bf_c      -> setToolTip(tr("Busy Function free parameter C"));
    tabOutputButtonParameter_bf_chi2   -> setToolTip(tr("Reduced &chi;<sup>2</sup> of Busy Function fit"));
    tabOutputButtonParameter_bf_flag   -> setToolTip(tr("Success flag of Busy Function fit"));
    tabOutputButtonParameter_bf_f_int  -> setToolTip(tr("Integrated flux derived from Busy Function fit"));
    tabOutputButtonParameter_bf_f_peak -> setToolTip(tr("Peak flux density derived from Busy Function fit"));
    tabOutputButtonParameter_bf_w      -> setToolTip(tr("Busy Function free parameter W"));
    tabOutputButtonParameter_bf_w20    -> setToolTip(tr("Line width (w<sub>20</sub>) derived from Busy Function fit"));
    tabOutputButtonParameter_bf_w50    -> setToolTip(tr("Line width (w<sub>50</sub>) derived from Busy Function fit"));
    tabOutputButtonParameter_bf_xe     -> setToolTip(tr("Busy Function free parameter X<sub>e</sub>"));
    tabOutputButtonParameter_bf_xp     -> setToolTip(tr("Busy Function free parameter X<sub>p</sub>"));
    tabOutputButtonParameter_bf_z      -> setToolTip(tr("Line centroid derived from Busy Function fit"));
    tabOutputButtonParameter_ell_maj   -> setToolTip(tr("Major axis of ellipse fitted to source"));
    tabOutputButtonParameter_ell_min   -> setToolTip(tr("Minor axis of ellipse fitted to source"));
    tabOutputButtonParameter_ell_pa    -> setToolTip(tr("Position angle of ellipse fitted to source"));
    tabOutputButtonParameter_f_peak    -> setToolTip(tr("Peak flux density"));
    tabOutputButtonParameter_f_int     -> setToolTip(tr("Integrated flux"));
    tabOutputButtonParameter_f_wm50    -> setToolTip(tr("Integrated flux within w<sub>m50</sub> line width"));
    tabOutputButtonParameter_rms       -> setToolTip(tr("RMS noise level"));
    tabOutputButtonParameter_w20       -> setToolTip(tr("Line width (w<sub>20</sub>)"));
    tabOutputButtonParameter_w50       -> setToolTip(tr("Line width (w<sub>50</sub>)"));
    tabOutputButtonParameter_wm50      -> setToolTip(tr("Line width (w<sub>m50</sub>)"));
    tabOutputButtonParameter_ra        -> setToolTip(tr("Right ascension of flux-weighted centre"));
    tabOutputButtonParameter_dec       -> setToolTip(tr("Declination of flux-weighted centre"));
    tabOutputButtonParameter_lon       -> setToolTip(tr("Longitude of flux-weighted centre"));
    tabOutputButtonParameter_lat       -> setToolTip(tr("Latitude of flux-weighted centre"));
    tabOutputButtonParameter_freq      -> setToolTip(tr("Frequency of flux-weighted centre"));
    tabOutputButtonParameter_velo      -> setToolTip(tr("Radial velocity of flux-weighted centre"));
    
    tabOutputWidgetParameters = new QWidget(tabOutputGroupBox2);
    tabOutputLayoutParameters = new QGridLayout();
    tabOutputLayoutParameters->setContentsMargins(0, 0, 0, 0);
    tabOutputLayoutParameters->setHorizontalSpacing(20);
    tabOutputLayoutParameters->setVerticalSpacing(5);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_id,        0, 0);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_x,         1, 0);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_y,         2, 0);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_z,         3, 0);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_x_geo,     4, 0);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_y_geo,     5, 0);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_z_geo,     6, 0);
    
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_x_min,     0, 1);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_x_max,     1, 1);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_y_min,     2, 1);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_y_max,     3, 1);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_z_min,     4, 1);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_z_max,     5, 1);
    
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_n_pix,     0, 2);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_n_pos,     1, 2);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_n_neg,     2, 2);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_snr_min,   3, 2);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_snr_max,   4, 2);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_snr_sum,   5, 2);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_rms,       6, 2);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_rel,       7, 2);
    
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_ra,        0, 3);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_dec,       1, 3);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_lon,       2, 3);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_lat,       3, 3);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_freq,      4, 3);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_velo,      5, 3);
    
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_w20,       0, 4);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_w50,       1, 4);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_wm50,      2, 4);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_f_peak,    3, 4);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_f_int,     4, 4);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_f_wm50,    5, 4);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_ell_maj,   6, 4);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_ell_min,   7, 4);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_ell_pa,    8, 4);
    
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_a,      0, 5);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_b1,     1, 5);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_b2,     2, 5);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_c,      3, 5);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_xe,     4, 5);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_xp,     5, 5);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_w,      6, 5);
    
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_z,      0, 6);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_w20,    1, 6);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_w50,    2, 6);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_f_peak, 3, 6);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_f_int,  4, 6);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_chi2,   5, 6);
    tabOutputLayoutParameters->addWidget(tabOutputButtonParameter_bf_flag,   6, 6);
    
    tabOutputWidgetParameters->setLayout(tabOutputLayoutParameters);
    
    tabOutputForm2->addRow(tr(""), tabOutputWidgetParameters);
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
    toolBoxOP->addItem(tabOutputGroupBox2, iconTaskReject, tr("Output Parameters"));
    
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
    //outputText->setToolTip(tr("Pipeline messages"));
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
    //outputProgress->setToolTip(tr("Pipeline progress"));
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
    this->addDockWidget(Qt::TopDockWidgetArea, dockWidgetOutput);
    
    // Set up actions
    // --------------
    
    actionOpen = new QAction(tr("Open..."), this);
    actionOpen->setShortcuts(QKeySequence::Open);
    actionOpen->setIcon(iconDocumentOpen);
    connect(actionOpen, SIGNAL(triggered()), this, SLOT(loadSettings()));
    
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
    
    // Set up menu
    // -----------
    
    menuFile = new QMenu(tr("&File"), this);
    menuFile->addAction(actionDefault);
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
    
    menuView = new QMenu(tr("&Analysis"), this);
    menuView->addAction(actionShowCatalogue);
    
    menuSettings = new QMenu(tr("&Settings"), this);
    menuSettings->addAction(dockWidgetOutput->toggleViewAction());
    menuSettings->addAction(toolBar->toggleViewAction());
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
    this->setWindowTitle(tr("SoFiA"));
    this->setCentralWidget(widgetMain);
    this->resize(600, 300);
    this->setWindowIcon(iconSoFiA);
    
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
    tabInputFieldCatalog->setWhatsThis(tr("<h3>optical.sourceCatalogue</h3><p>This defines the full path and file name of the input catalogue required for catalogue-based source finding (see parameter <b>steps.doOptical</b>). There is no default.</p>"));
    tabInputFieldSpatialSize->setWhatsThis(tr("<h3>optical.spatSize</h3><p>This defines the <b>spatial</b> size of the sub-cube to be searched around each catalogue position. The size must be specified in the <b>native units</b> of the data cube, e.g. in degrees.</p>"));
    tabInputFieldSpectralSize->setWhatsThis(tr("<h3>optical.specSize</h3><p>This defines the <b>spectral</b> size of the sub-cube to be searched around each catalogue position. The size must be specified in the <b>native units</b> of the data cube, e.g. in km/s or Hz.</p>"));
    tabInputFieldMultiCat->setWhatsThis(tr("<h3>optical.storeMultiCat</h3><p>If set to <b>true</b>, a separate output catalogue will be created for each input position, containing only the sources found in that subcube. In addition, a single, merged catalogue will also be created. By default this parameter is set to <b>false</b>, in which case only a single output catalogue file is generated.</p>"));
    tabParametrisationButtonDilateMask->setWhatsThis(tr("<h3>parameters.dilateMask</h3><p>Run the mask optimisation algorithm based on spatially <b>dilating</b> the initial mask to achieve more accurate flux measurements.</p>"));
    tabParametrisationButtonBusyFunction->setWhatsThis(tr("<h3>parameters.fitBusyFunction</h3><p>Fit the Busy Function (<a href=\"http://adsabs.harvard.edu/abs/2014MNRAS.438.1176W\">Westmeier et al. 2014</a>) to the integrated spectrum for more accurate parameterisation.</p>"));
    tabParametrisationButtonMaskOpt->setWhatsThis(tr("<h3>parameters.optimiseMask</h3><p>Run the mask optimisation algorithm based on fitting and <b>growing</b> ellipses to achieve more accurate flux measurements.</p>"));
    tabParametrisationFieldRelKernel->setWhatsThis(tr("<h3>reliability.kernel</h3><p>Size of 3D smoothing kernel in log(parameter) space (see <b>reliability.parSpace</b>).</p>"));
    tabParametrisationFieldRelPlot->setWhatsThis(tr("<h3>reliability.makePlot</h3><p>If set to <b>true</b>, a PDF file showing the distribution of positive and negative detections in parameter space will be created for diagnostic purposes.</p>"));
    tabParametrisationFieldRelMin->setWhatsThis(tr("<h3>reliability.threshold</h3><p>Discard sources whose reliability is below this threshold.</p>"));
    tabInFilterFieldEdgeX->setWhatsThis(tr("<h3>scaleNoise.edgeX</h3><p>Size of edge (in pixels) to be excluded in first coordinate.</p>"));
    tabInFilterFieldEdgeY->setWhatsThis(tr("<h3>scaleNoise.edgeY</h3><p>Size of edge (in pixels) to be excluded in second coordinate.</p>"));
    tabInFilterFieldEdgeZ->setWhatsThis(tr("<h3>scaleNoise.edgeZ</h3><p>Size of edge (in pixels) to be excluded in third coordinate.</p>"));
    tabInFilterFieldScaleX->setWhatsThis(tr("<h3>scaleNoise.scaleX</h3><p>Noise normalisation in first (spatial) dimension.</p>"));
    tabInFilterFieldScaleY->setWhatsThis(tr("<h3>scaleNoise.scaleY</h3><p>Noise normalisation in second (spatial) dimension.</p>"));
    tabInFilterFieldScaleZ->setWhatsThis(tr("<h3>scaleNoise.scaleZ</h3><p>Noise normalisation in third (spectral) dimension.</p>"));
    tabInFilterFieldStatistic->setWhatsThis(tr("<h3>scaleNoise.statistic</h3><p>Statistic used to measure the noise. This can be median absolute deviation (<b>mad</b>), standard deviation (<b>std</b>) or Gaussian fit to negative fluxes (<b>negative</b>).</p>"));
    tabSourceFindingFieldEdgeMode->setWhatsThis(tr("<h3>SCfind.edgeMode</h3><p>Behaviour near the edge of the cube. The following values are possible:<p><ul><li><b>constant:</b> assume constant value of 0</li><li><b>nearest:</b> assume constant value equal to edge pixel</li><li><b>reflect:</b> mirror values at edge, thereby including the edge pixel itself</li><li><b>mirror:</b> mirror values at position of outermost pixel, thereby excluding the edge pixel itself</li><li><b>wrap:</b> copy values from opposite edge of the array</li></ul>"));
    tabSourceFindingFieldKernels->setWhatsThis(tr("<h3>SCfind.kernels</h3><p>List of kernels to be used for smoothing. The format is:</p><p style=\"font-family:monospace;\">[[dx, dy, dz, 'type'], ...]</p><p>where <b>dx</b>, <b>dy</b>, and <b>dz</b> are the spatial and spectral kernel sizes (FWHM), and <b>'type'</b> can be boxcar (<b>'b'</b>) or Gaussian (<b>'g'</b>). Note that 'type' only applies to the spectral axis, and the spatial kernel is always Gaussian.</p>"));
    tabSourceFindingFieldKunit->setWhatsThis(tr("<h3>SCfind.kernelUnit</h3><p>Are kernel parameters specified in <b>pixel</b> or <b>world</b> coordinates?</p>"));
    tabSourceFindingFieldRmsMode->setWhatsThis(tr("<h3>SCfind.rmsMode</h3><p>Noise determination method: Gaussian fit to negative flux histogram (<b>negative</b>), median absolute deviation (<b>mad</b>), or standard deviation (<b>std</b>).</p>"));
    tabSourceFindingFieldThreshold->setWhatsThis(tr("<h3>SCfind.threshold</h3><p>Flux threshold relative to the noise level.</p>"));
    tabInFilterFieldBorder->setWhatsThis(tr("<h3>smooth.edgeMode</h3><p>Behaviour near the edge of the cube. The following options are supported:</p><ul><li><b>constant:</b> assume constant value of 0</li><li><b>nearest:</b> assume constant value equal to edge pixel</li><li><b>reflect:</b> mirror values at edge, thereby including the edge pixel itself</li><li><b>mirror:</b> mirror values at position of outermost pixel, thereby excluding the edge pixel itself</li><li><b>wrap:</b> copy values from opposite edge of the array</li></ul>"));
    tabInFilterFieldKernel->setWhatsThis(tr("<h3>smooth.kernel</h3><p>Type of smoothing kernel used in both spatial and spectral smoothing. Can be <b>gaussian</b>, <b>boxcar</b> or <b>median</b>.</p>"));
    tabInFilterFieldSmoothingSpatialLon->setWhatsThis(tr("<h3>smooth.kernelX</h3><p>Kernel size in pixels for first coordinate. For Gaussian kernels the value refers to the FWHM.</p>"));
    tabInFilterFieldSmoothingSpatialLat->setWhatsThis(tr("<h3>smooth.kernelY</h3><p>Kernel size in pixels for second coordinate. For Gaussian kernels the value refers to the FWHM.</p>"));
    tabInFilterFieldSmoothingSpectral->setWhatsThis(tr("<h3>smooth.kernelZ</h3><p>Kernel size in pixels for third coordinate. For Gaussian kernels the value refers to the FWHM.</p>"));
    tabSourceFindingGroupBox3->setWhatsThis(tr("<h3>steps.doCNHI</h3><p>Run the Characterised Noise HI (CNHI) source finder (<a href=\"http://adsabs.harvard.edu/abs/2012PASA...29..251J\">Jurek 2012</a>).</p>"));
    tabOutputButtonCubelets->setWhatsThis(tr("<h3>steps.doCubelets</h3><p>Create and save data products for each individual source, including sub-cubes, moment 0, 1 and 2 maps, integrated spectra and position&ndash;velocity diagrams.</p>"));
    tabInputGroupBox3->setWhatsThis(tr("<h3>steps.doFlag</h3><p>Flag pixel and channel ranges before proceeding. Details are specified with the <b>flag.regions</b> option.</p>"));
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
    tabOutputButtonFilteredCube->setWhatsThis(tr("<h3>steps.doWriteFilteredCube</h3><p>Save a copy of the filtered data cube. Note that this will only work if at least one of the input filters was applied.</p>"));
    tabOutputButtonMask->setWhatsThis(tr("<h3>steps.doWriteMask</h3><p>Save source mask cube.</p>"));
    tabSourceFindingFieldClipMethod->setWhatsThis(tr("<h3>threshold.clipMethod</h3><p>Define whether the flux threshold is <b>relative</b> to the noise level or in <b>absolute</b> flux units.</p>"));
    tabSourceFindingFieldRmsMode2->setWhatsThis(tr("<h3>threshold.rmsMode</h3><p>Noise determination method: Gaussian fit to negative flux histogram (<b>negative</b>), median absolute deviation (<b>mad</b>), or standard deviation (<b>std</b>).</p>"));
    tabSourceFindingFieldThreshold2->setWhatsThis(tr("<h3>threshold.threshold</h3><p>Absolute or relative flux threshold for detections (see <b>threshold.clipMethod</b>).</p>"));
    tabInFilterField2d1dIterations->setWhatsThis(tr("<h3>wavelet.iterations</h3><p>Number of iterations in the reconstruction process.</p>"));
    tabInFilterField2d1dPositivity->setWhatsThis(tr("<h3>wavelet.positivity</h3><p>If <b>true</b>, include only positive wavelet components in the decomposition.</p>"));
    tabInFilterField2d1dScaleXY->setWhatsThis(tr("<h3>wavelet.scaleXY</h3><p>Number of spatial scales used in the decomposition. The default value of <b>-1</b> will automatically determine the appropriate number of scales based on the data cube.</p>"));
    tabInFilterField2d1dScaleZ->setWhatsThis(tr("<h3>wavelet.scaleZ</h3><p>Number of spectral scales used in the decomposition. The default value of <b>-1</b> will automatically determine the appropriate number of scales based on the data cube.</p>"));
    tabInFilterField2d1dThreshold->setWhatsThis(tr("<h3>wavelet.threshold</h3><p>Flux threshold used in the wavelet reconstruction in multiples of the rms noise. Note that this threshold only determines which wavelet components are added to the decomposed cube; any source finding will be done separately, using a different flux threshold.</p>"));
    tabOutputFieldBaseName->setWhatsThis(tr("<h3>writeCat.baseName</h3><p>By default, SoFiA will use the file name of the input data cube as the template for all output data file names, usually extended by a specific postfix based on the actual output data product. By setting the <b>writeCat.basename</b> parameter, a different template for all output file names can be specified.</p>"));
    tabOutputButtonCompress->setWhatsThis(tr("<h3>writeCat.compress</h3><p>If set to true, use <a href=\"http://www.gzip.org/\">gzip</a> to compress all output files.</p>"));
    tabOutputFieldDirectory->setWhatsThis(tr("<h3>writeCat.outputDir</h3><p>Optional directory path to which all output files are written. If not specified, the directory of the input cube will be used by default.</p>"));
    tabOutputGroupBox2->setWhatsThis(tr("<h3>writeCat.parameters</h3><p>List of parameters to appear in source catalogue. Format:</p><p>['par1', 'par2', ...]</p><p>An asterisk, <span style=\"font-family:monospace;\">['*']</span>, means that all parameters are written to the catalogue (default behaviour).</p>"));
    tabOutputButtonASCII->setWhatsThis(tr("<h3>writeCat.writeASCII</h3><p>Write catalogue in ASCII format.</p>"));
    tabOutputButtonSQL->setWhatsThis(tr("<h3>writeCat.writeSQL</h3><p>Write catalogue in SQL format (not yet implemented).</p>"));
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
        QMessageBox messageBox(this);
        messageBox.setWindowTitle(tr("SoFiA - Quit"));
        messageBox.setText(tr("<p>Do you wish to exit from SoFiA?</p><p>Note that any unsaved changes will be discarded.</p>"));
        messageBox.setStandardButtons(QMessageBox::Cancel | QMessageBox::Ok);
        messageBox.setDefaultButton(QMessageBox::Ok);
        messageBox.setIcon(QMessageBox::Warning);
        int choice = messageBox.exec();
        
        if(choice == QMessageBox::Ok) event->accept();
        else event->ignore();
    }
    else
    {
        event->ignore();
        
        QString messageText = tr("<p>The pipeline is still running.</p><p>If you wish to exit from SoFiA, you will need to either manually abort the pipeline or wait until the pipeline run has finished.</p>");
        QString statusText = tr("Pipeline still running.");
        showMessage(MESSAGE_WARNING, messageText, statusText);
    }
    
    return;
}
