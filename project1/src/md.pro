TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    System.cpp \
    StatisticsSampler.cpp \
    md.cpp \
    InitialConditions.cpp \
    Atom.cpp \
    random.cpp \
    Cell.cpp \
    thermostat.cpp \
    threadcontrol.cpp

HEADERS += \
    System.h \
    StatisticsSampler.h \
    inlines.h \
    Atom.h \
    random.h \
    Cell.h \
    thermostat.h \
    threadcontrol.h

mac {
    CONFIG -= app_bundle
    LIBS   += -larmadillo
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
}

unix:!mac {
    LIBS   += -larmadillo
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
}

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
