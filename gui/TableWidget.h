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
