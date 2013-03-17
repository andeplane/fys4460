TEMPLATE = app
CONFIG += console
CONFIG -= qt

release {
    DEFINES += ARMA_NO_DEBUG
}

DEFINES += MPI_ENABLED

SOURCES += \
    System.cpp \
    StatisticsSampler.cpp \
    md.cpp \
    random.cpp \
    thermostat.cpp \
    CVector.cpp \
    CUtil.cpp \
    CMath.cpp \
    unitconverter.cpp \
    settings.cpp \
    mdio.cpp

HEADERS += \
    System.h \
    StatisticsSampler.h \
    inlines.h \
    random.h \
    thermostat.h \
    CVector.h \
    CUtil.h \
    CIniFile.h \
    CMatrix.h \
    CMath.h \
    Stdafx.h \
    unitconverter.h \
    settings.h \
    mdio.h

mac {
    QMAKE_CXX = icc
    CONFIG -= app_bundle
    LIBS   += -larmadillo -llapack -lblas
    INCLUDEPATH += /usr/local/Cellar/boost/1.49.0/include/boost
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
}

unix:!mac {
    QMAKE_CXX = mpic++
    LIBS   += -larmadillo -llapack -lblas
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
}

# MPI Settings
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
#QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
