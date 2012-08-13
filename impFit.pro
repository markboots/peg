QT     -= gui core

CONFIG -= app_bundle
TARGET = impFit

QMAKE_CXXFLAGS += -c -g -Wall -fopenmp
QMAKE_LFLAGS += -fopenmp

INCLUDEPATH += /Users/mboots/dev/gsl-install/include

LIBS += -L/Users/mboots/dev/gsl-install/lib -lgsl -lgslcblas

HEADERS += src/PEG.h \
	src/PESolver.h \
	src/PEMainSupport.h

SOURCES += src/PEG.cpp\
	src/PESolver.cpp \
	src/PEMainSupport.cpp \
    src/impFit.cpp
