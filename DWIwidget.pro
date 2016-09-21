#-------------------------------------------------
#
# Project created by QtCreator 2016-08-09T16:55:43
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = DWIwidget
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h

FORMS    += mainwindow.ui

#INCLUDEPATH += $$E:\Tech318\MITK\MITK\Build\ep\src\CTK-build\CTK-build\Libs\DICOM\Widgets


#win32:CONFIG(release, debug|release): LIBS += -L$$C:\Tech318\VTK-7.0.0\Build\Install\lib\ -lvtkGUISupportQt-7.0
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$C:\Tech318\VTK-7.0.0\Build\Install\lib\ -lvtkGUISupportQt-7.0d

#win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\libvtkGUISupportQt-7.0.a
#else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\libvtkGUISupportQt-7.0d.a
#else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\vtkGUISupportQt-7.0.lib
#else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\vtkGUISupportQt-7.0d.lib

#win32:CONFIG(release, debug|release): LIBS += -L$$C:\Tech318\VTK-7.0.0\Build\Install\lib\ -lvtkIOImage-7.0
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$C:\Tech318\VTK-7.0.0\Build\Install\lib\ -lvtkIOImage-7.0d

#win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\libvtkIOImage-7.0.a
#else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\libvtkIOImage-7.0d.a
#else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\vtkIOImage-7.0.lib
#else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\vtkIOImage-7.0d.lib

#win32:CONFIG(release, debug|release): LIBS += -L$$C:\Tech318\VTK-7.0.0\Build\Install\lib\ -lvtkInteractionImage-7.0
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$C:\Tech318\VTK-7.0.0\Build\Install\lib\ -lvtkInteractionImage-7.0d

#win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\libvtkInteractionImage-7.0.a
#else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\libvtkInteractionImage-7.0d.a
#else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\vtkInteractionImage-7.0.lib
#else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$C:\Tech318\VTK-7.0.0\Build\Install\lib\vtkInteractionImage-7.0d.lib

win32: LIBS += -L$$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/ -lvtkGUISupportQt-7.0

INCLUDEPATH += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/include
DEPENDPATH += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/include

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/vtkGUISupportQt-7.0.lib
else:win32-g++: PRE_TARGETDEPS += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/libvtkGUISupportQt-7.0.a

#win32: LIBS += -L$$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/ -lvtkIOImage-7.0

#INCLUDEPATH += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/include
#DEPENDPATH += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/include

#win32:!win32-g++: PRE_TARGETDEPS += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/vtkIOImage-7.0.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/libvtkIOImage-7.0.a

#win32: LIBS += -L$$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/ -lvtkInteractionImage-7.0

#INCLUDEPATH += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/include
#DEPENDPATH += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/include

#win32:!win32-g++: PRE_TARGETDEPS += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/vtkInteractionImage-7.0.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/../../../../Tech318/VTK-7.0.0/Build/Install/lib/libvtkInteractionImage-7.0.a


#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/E:/Tech318/MITK/MITK/Build/ep/src/CTK-build/CTK-build/bin/release/ -lCTKDICOMWidgets
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/E:/Tech318/MITK/MITK/Build/ep/src/CTK-build/CTK-build/bin/debug/ -lCTKDICOMWidgets

#INCLUDEPATH += $$PWD/E:/Tech318/MITK/MITK/Build/ep/src/CTK-build/CTK-build/bin/Release
#DEPENDPATH += $$PWD/E:/Tech318/MITK/MITK/Build/ep/src/CTK-build/CTK-build/bin/Release

#win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/E:/Tech318/MITK/MITK/Build/ep/src/CTK-build/CTK-build/bin/release/libCTKDICOMWidgets.a
#else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/E:/Tech318/MITK/MITK/Build/ep/src/CTK-build/CTK-build/bin/debug/libCTKDICOMWidgets.a
#else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/E:/Tech318/MITK/MITK/Build/ep/src/CTK-build/CTK-build/bin/release/CTKDICOMWidgets.lib
#else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/E:/Tech318/MITK/MITK/Build/ep/src/CTK-build/CTK-build/bin/debug/CTKDICOMWidgets.lib
