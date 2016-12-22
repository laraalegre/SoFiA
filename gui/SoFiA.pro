#########################################
# Created by hand from a qmake template #
#########################################

TEMPLATE = app
TARGET = SoFiA
DEPENDPATH += .
INCLUDEPATH += .
QT += xml
LIBS += -lz
QMAKE_CXXFLAGS += -std=c++11 # for Qt 4

# Include module 'widgets' for Qt 5 or greater:
greaterThan(QT_MAJOR_VERSION, 4)
{
    QT += widgets
    CONFIG += c++11 # for Qt 5
}

# Input
HEADERS   += HelpBrowser.h \
             TableWidget.h \
             WidgetSpreadsheet.h \
             Fips.hpp \
             WidgetDataViewer.h \
             SoFiA.h
SOURCES   += HelpBrowser.cpp \
             TableWidget.cpp \
             WidgetSpreadsheet.cpp \
             Fips.cpp \
             WidgetDataViewer.cpp \
             SoFiA.cpp \
             main.cpp
RESOURCES  = SoFiA.qrc
