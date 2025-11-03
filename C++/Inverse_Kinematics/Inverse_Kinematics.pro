/*
 * Copyright (c) 2025 Yuhan Chen (chenyuhan19930920@163.com)
 * Licensed under the MIT License.
 * See LICENSE file in the project root for full license information.
*/


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
