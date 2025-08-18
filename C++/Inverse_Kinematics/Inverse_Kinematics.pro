TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        InverseKinematics.cpp \
        test.cpp


# eigen
INCLUDEPATH += $$PWD/dependence/eigen-3.4.0

HEADERS += \
    InverseKinematics.h
