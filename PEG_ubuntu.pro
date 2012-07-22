QT     -= gui core

TARGET = pegSerial

QMAKE_CXXFLAGS += -c -g -Wall -fopenmp
QMAKE_LFLAGS += -fopenmp

#INCLUDEPATH += /Users/mboots/dev/gsl-install/include

LIBS += -lgsl -lgslcblas
GSL_LIB_DIR = /usr/local/lib

QMAKE_LFLAGS_DEBUG += "-Wl,-rpath,$$GSL_LIB_DIR"
QMAKE_LFLAGS_RELEASE += "-Wl,-rpath,$$GSL_LIB_DIR"

HEADERS += src/PEG.h \
	src/PESolver.h \
	src/PEMainSupport.h

SOURCES += src/PEG.cpp\
	src/PESolver.cpp \
	src/PEMainSupport.cpp \
	src/mainSerial.cpp
