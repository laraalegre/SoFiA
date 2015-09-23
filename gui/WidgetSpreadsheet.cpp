/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 0.5.0 (WidgetSpreadsheet.cpp) - Source Finding Application     ///
/// Copyright (C) 2013-2015 Tobias Westmeier                             ///
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

#include "TableWidget.h"
#include "WidgetSpreadsheet.h"

// ----------- //
// CONSTRUCTOR //
// ----------- //

WidgetSpreadsheet::WidgetSpreadsheet(QWidget *parent)
{
	this->setParent(parent);
	
	// Create button icons:
	iconDialogClose.addFile(QString(":/icons/22/dialog-close.png"), QSize(22, 22));
	iconDialogClose.addFile(QString(":/icons/16/dialog-close.png"), QSize(16, 16));
	iconDialogClose = QIcon::fromTheme("dialog-close", iconDialogClose);
	
	iconViewRefresh.addFile(QString(":/icons/22/view-refresh.png"), QSize(22, 22));
	iconViewRefresh.addFile(QString(":/icons/16/view-refresh.png"), QSize(16, 16));
	iconViewRefresh = QIcon::fromTheme("view-refresh", iconViewRefresh);
	
	// Create table widget:
	tableWidget = new TableWidget(this);
	tableWidth  = 0;
	tableHeight = 0;
	tableWidget->setColumnCount(tableWidth);
	tableWidget->setRowCount(tableHeight);
	tableWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	tableWidget->setStyleSheet("QTableWidget::item {padding:5px 10px;}");
	
	// Create control widget:
	widgetControls = new QWidget(this);
	
	widgetSort = new QWidget(widgetControls);
	buttonSort = new QComboBox(widgetSort);
	buttonSort->setMinimumWidth(100);
	buttonSort->setEnabled(false);
	connect(buttonSort, SIGNAL(currentIndexChanged(int)), this, SLOT(sortTable(int)));
	layoutSort = new QFormLayout;
	layoutSort->addRow(tr("Sort by"), buttonSort);
	layoutSort->setContentsMargins(0, 0, 0, 0);
	layoutSort->setSpacing(5);
	widgetSort->setLayout(layoutSort);
	
	buttonOrder = new QCheckBox(tr("descending "), widgetControls);
	buttonOrder->setChecked(false);
	connect(buttonOrder, SIGNAL(toggled(bool)), this, SLOT(changeSortOrder()));
	
	buttonReload = new QPushButton(tr("Reload"), widgetControls);
	buttonReload->setIcon(iconViewRefresh);
	buttonReload->setEnabled(false);
	connect(buttonReload, SIGNAL(clicked()), this, SLOT(reloadCatalog()));
	
	buttonClose = new QPushButton(tr("Close"), widgetControls);
	buttonClose->setIcon(iconDialogClose);
	connect(buttonClose, SIGNAL(clicked()), this, SLOT(close()));
	
	layoutControls = new QHBoxLayout;
	layoutControls->setContentsMargins(5, 5, 5, 5);
	layoutControls->setSpacing(10);
	layoutControls->addWidget(widgetSort);
	layoutControls->addWidget(buttonOrder);
	layoutControls->addStretch();
	layoutControls->addWidget(buttonReload);
	layoutControls->addWidget(buttonClose);
	widgetControls->setLayout(layoutControls);
	widgetControls->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);
	
	// Set up main layout:
	mainLayout = new QVBoxLayout;
	mainLayout->setContentsMargins(0, 0, 0, 0);
	mainLayout->setSpacing(0);
	mainLayout->addWidget(tableWidget);
	mainLayout->addWidget(widgetControls);
	
	// Set up main window:
	this->setLayout(mainLayout);
	this->setWindowFlags(Qt::Window);
	this->setWindowTitle(tr("SoFiA - Source Catalogue"));
	this->resize(720, 480);
	
	return;
}



// --------------------- //
// LOAD SOURCE CATALOGUE //
//---------------------- //

int WidgetSpreadsheet::loadCatalog(QString &filename)
{
	if(filename.isEmpty()) return 1;
	
	currentFileName = filename;
	tableWidget->clear();              // Clear the existing table. This will automatically delete all QTableWidgetItem objects
	                                   // created before, so no memory will be leaked.
	tableWidth  = 0;
	tableHeight = 0;
	buttonSort->clear();               // Clear sort button as well.
	buttonSort->setEnabled(false);
	buttonOrder->setEnabled(false);
	buttonReload->setEnabled(false);
	
	QDomDocument catalogue(currentFileName);
	QFile file(currentFileName);
	
	if(!file.open(QIODevice::ReadOnly)) return 1;
	
	if(!catalogue.setContent(&file))
	{
		file.close();
		return 1;
	}
	
	file.close();
	
	// Get all tags named FIELD (contains header information):
	QDomNodeList headerTags = catalogue.elementsByTagName("FIELD");
	if(headerTags.isEmpty()) return 1;
	
	tableWidth = headerTags.size();
	tableWidget->setColumnCount(tableWidth);
	
	QStringList datatype;    // This will store the data type of each column for later use.
	
	// Extract all header names:
	for(int i = 0; i < headerTags.size(); i++)
	{
		QDomNode field = headerTags.item(i);
		QDomNamedNodeMap fieldAttributes = field.attributes();
		QDomNode attribute = fieldAttributes.namedItem("name");
		QDomNode unit      = fieldAttributes.namedItem("unit");
		datatype.append((fieldAttributes.namedItem("datatype")).nodeValue());
		
		QString headerText = attribute.nodeValue().trimmed();
		
		buttonSort->addItem(headerText);
		
		if(not unit.isNull() and unit.nodeValue() != "-")
		{
			headerText.append("\n(");
			headerText.append(unit.nodeValue());
			headerText.append(")");
		}
		//else headerText.append("\n");
		
		QTableWidgetItem* headerItem = new QTableWidgetItem(headerText, QTableWidgetItem::Type);
		tableWidget->setHorizontalHeaderItem(i, headerItem);
	}
	
	// Get all tags named TR (data rows):
	QDomNodeList rowTags = catalogue.elementsByTagName("TR");
	if(rowTags.isEmpty()) return 1;
	
	tableHeight = rowTags.size();
	tableWidget->setRowCount(tableHeight);
	
	// Extract all fields (TD) within each row:
	for(int i = 0; i < rowTags.size(); i++)
	{
		QDomNode row = rowTags.item(i);
		QDomNodeList cellTags = row.childNodes();
		
		for(int j = 0; j < tableWidth; j++)
		{
			if(j < cellTags.size())
			{
				QDomNode cell = cellTags.item(j);
				
				if(not cell.isNull() and cell.hasChildNodes())
				{
					QDomNode entry = cell.firstChild();
					QString  text  = entry.nodeValue().trimmed();
					double   value = text.toDouble();
					QTableWidgetItem* cellItem = new QTableWidgetItem();
					
					if(datatype[j] == "char") cellItem->setData(Qt::DisplayRole, text);
					else                      cellItem->setData(Qt::DisplayRole, value);
					
					cellItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
					cellItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
					tableWidget->setItem(i, j, cellItem);
				}
			}
		}
	}
	
	tableWidget->resizeColumnsToContents();
	tableWidget->resizeRowsToContents();
	
	buttonSort->setEnabled(true);
	buttonOrder->setEnabled(true);
	buttonOrder->setChecked(false);
	buttonReload->setEnabled(true);
	
	return 0;
}



// -----------------  //
// SLOT TO SORT TABLE //
// -----------------  //

void WidgetSpreadsheet::sortTable(int column)
{
	if(buttonOrder->isChecked()) tableWidget->sortItems(column, Qt::DescendingOrder);
	else                         tableWidget->sortItems(column, Qt::AscendingOrder);
	
	return;
}



// ------------------------  //
// SLOT TO CHANGE SORT ORDER //
// ------------------------  //

void WidgetSpreadsheet::changeSortOrder()
{
	sortTable(buttonSort->currentIndex());
	
	return;
}



// -----------------------  //
// SLOT TO RELOAD CATALOGUE //
// -----------------------  //

void WidgetSpreadsheet::reloadCatalog()
{
	if(!currentFileName.isEmpty() and loadCatalog(currentFileName) != 0) std::cerr << "Failed to reload catalogue.\n";
	
	return;
}



// ----------------- //
// CLEAN UP ON CLOSE //
// ----------------- //

void WidgetSpreadsheet::closeEvent(QCloseEvent *event)
{
	emit widgetClosed();
	event->accept();
}
