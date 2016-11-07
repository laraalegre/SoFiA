#########################################
# Created by hand from a qmake template #
#########################################

TEMPLATE = app
TARGET = SoFiA
DEPENDPATH += .
INCLUDEPATH += .
QT += xml
LIBS += -lz

# Include module 'widgets' for Qt 5 or greater:
greaterThan(QT_MAJOR_VERSION, 4)
{
    QT += widgets
}

# Input
HEADERS   += HelpBrowser.h \
             TableWidget.h \
             WidgetSpreadsheet.h \
             SoFiA.h
SOURCES   += HelpBrowser.cpp \
             TableWidget.cpp \
             WidgetSpreadsheet.cpp \
             SoFiA.cpp \
             main.cpp
RESOURCES  = SoFiA.qrc
