#-------------------------------------------------
#
# Project created by QtCreator 2013-10-02T22:33:16
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Motion
TEMPLATE = app


SOURCES += main.cpp \
    Math/cubespline.cpp \
    Math/vector3d.cpp \
    Math/vector2d.cpp \
    Math/matrix4x4.cpp \
    Math/matrix3x3.cpp \
    Math/matrix2x2.cpp \
    Graph/graphicwidget.cpp \
    Graph/graphdata.cpp \
    Math/rkf45.cpp \
    RightFuncs.cpp

HEADERS  += \
    Math/cubespline.h \
    Math/vector3d.h \
    Math/vector2d.h \
    Math/matrix4x4.h \
    Math/matrix3x3.h \
    Math/matrix2x2.h \
    Math/allmath.h \
    Graph/graphicwidget.h \
    Graph/graphdata.h \
    Math/rkf45.h \
    RightFuncs.h

LIBS += \
    -LC:\Qt\Qt5.1.0\5.1.0\mingw48_32\lib\libboost -lboost_math_tr1

