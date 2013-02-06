TEMPLATE = app
CONFIG += console
CONFIG -= qt
CONFIG -= app_bundle

SOURCES += \
    System.cpp \
    StatisticsSampler.cpp \
    md.cpp \
    InitialConditions.cpp \
    Atom.cpp \
    random.cpp \
    Cell.cpp

HEADERS += \
    System.h \
    StatisticsSampler.h \
    inlines.h \
    Atom.h \
    random.h \
    Cell.h

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
