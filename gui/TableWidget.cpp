/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 1.3.2 (TableWidget.cpp) - Source Finding Application           ///
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

#include "TableWidget.h"

// ----------- //
// Constructor //
// ----------- //

TableWidget::TableWidget(QWidget *parent)
{
	this->setParent(parent);
	
	// Create icons:
	iconCopy.addFile(QString(":/icons/16/edit-copy.png"), QSize(16, 16));
	iconCopy.addFile(QString(":/icons/22/edit-copy.png"), QSize(22, 22));
	iconCopy = QIcon::fromTheme("edit-copy", iconCopy);
	
	iconSelectAll.addFile(QString(":/icons/16/edit-select-all.png"), QSize(16, 16));
	iconSelectAll.addFile(QString(":/icons/22/edit-select-all.png"), QSize(22, 22));
	iconSelectAll = QIcon::fromTheme("edit-select-all", iconSelectAll);
	
	iconClearSelection.addFile(QString(":/icons/16/edit-clear.png"), QSize(16, 16));
	iconClearSelection.addFile(QString(":/icons/22/edit-clear.png"), QSize(22, 22));
	iconClearSelection = QIcon::fromTheme("edit-clear", iconClearSelection);
	
	// Set up actions:
	actionCopy = new QAction(tr("Copy"), this);
	actionCopy->setShortcut(QKeySequence::Copy);
	actionCopy->setEnabled(false);
	actionCopy->setIcon(iconCopy);
	connect(actionCopy, SIGNAL(triggered()), this, SLOT(copy()));
	
	actionSelectAll = new QAction(tr("Select All"), this);
	actionSelectAll->setShortcut(QKeySequence::SelectAll);
	actionSelectAll->setEnabled(true);
	actionSelectAll->setIcon(iconSelectAll);
	connect(actionSelectAll, SIGNAL(triggered()), this, SLOT(selectEverything()));
	
	actionClearSelection = new QAction(tr("Clear Selection"), this);
	actionClearSelection->setShortcut(Qt::Key_Escape);
	actionClearSelection->setEnabled(false);
	actionClearSelection->setIcon(iconClearSelection);
	connect(actionClearSelection, SIGNAL(triggered()), this, SLOT(selectNothing()));
	
	// React to change in item selection:
	connect(this, SIGNAL(itemSelectionChanged()), this, SLOT(updateActions()));
	
	return;
}



// ---------------------------------------- //
// Slot to copy selected cells to clipboard //
// ---------------------------------------- //

void TableWidget::copy()
{
	// Get a list of all selected ranges:
	QList<QTableWidgetSelectionRange> selectedRanges = this->selectedRanges();
	if(selectedRanges.isEmpty()) return;
	
	// Establish the outer boundary of all selections:
	int leftColumn  = this->columnCount() - 1;
	int rightColumn = 0;
	int topRow      = this->rowCount() - 1;
	int bottomRow   = 0;
	
	for(int i = 0; i < selectedRanges.size(); i++)
	{
		QTableWidgetSelectionRange range = selectedRanges.at(i);
		
		if(range.leftColumn()  < leftColumn)  leftColumn  = range.leftColumn();
		if(range.rightColumn() > rightColumn) rightColumn = range.rightColumn();
		if(range.topRow()      < topRow)      topRow      = range.topRow();
		if(range.bottomRow()   > bottomRow)   bottomRow   = range.bottomRow();
	}
	
	if(bottomRow < topRow or rightColumn < leftColumn) return;
	
	// Loop through selection range and extract data:
	QString outputText;
	
	for(int i = topRow; i <= bottomRow; i++)
	{
		for(int j = leftColumn; j <= rightColumn; j++)
		{
			if(this->item(i, j)->isSelected())
			{
				outputText += this->item(i, j)->text();
			}
			
			if (j < rightColumn) outputText += "\t";
			else                 outputText += "\n";
		}
	}
	
	// Copy data to clipboard:
	QClipboard *clipboard = QApplication::clipboard();
	clipboard->setText(outputText);
	
	return;
}



// ------------------------ //
// Slot to select all cells //
// ------------------------ //

void TableWidget::selectEverything()
{
	this->selectAll();
	
	return;
}



// -------------------------- //
// Slot to deselect all cells //
// -------------------------- //

void TableWidget::selectNothing()
{
	this->selectionModel()->clearSelection();
	
	return;
}



// ---------------------------------- //
// Slot to update actions and buttons //
// ---------------------------------- //

void TableWidget::updateActions()
{
	QList<QTableWidgetSelectionRange> selectedRanges = this->selectedRanges();
	
	actionCopy->setEnabled(not selectedRanges.isEmpty());
	actionClearSelection->setEnabled(not selectedRanges.isEmpty());
	
	return;
}



// --------------------------------------------------- //
// Function to create context menu when right-clicking //
// --------------------------------------------------- //

void TableWidget::contextMenuEvent(QContextMenuEvent *event)
{
	QMenu menuContext(this);
	
	menuContext.addAction(actionCopy);
	menuContext.addSeparator();
	menuContext.addAction(actionSelectAll);
	menuContext.addAction(actionClearSelection);
	
	menuContext.exec(event->globalPos());
	
	return;
}



// ------------------------------------------------------------ //
// Function to override default key press event of QTableWidget //
// ------------------------------------------------------------ //

void TableWidget::keyPressEvent(QKeyEvent *event)
{
	if(event->matches(QKeySequence::Copy)) this->copy();
	if(event->matches(QKeySequence::SelectAll)) this->selectEverything();
	else if(event->key() == Qt::Key_Escape) this->selectNothing();
	else QTableWidget::keyPressEvent(event);
	
	return;
}
