TEMPLATE = app
CONFIG += console
CONFIG -= qt
CONFIG -= app_bundle

SOURCES += \
    System.cpp \
    StatisticsSampler.cpp \
    md.cpp \
    lib.cpp \
    InitialConditions.cpp \
    Atom.cpp

HEADERS += \
    System.h \
    StatisticsSampler.h \
    lib.h \
    inlines.h \
    initialConditions.h \
    Atom.h

mac {
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
