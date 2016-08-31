# include_directories(${CTK_INCLUDE_DIRS})

# MITK_CREATE_MODULE(
  # DEPENDS MitkCore
  # PACKAGE_DEPENDS
    # PUBLIC CTK|CTKDICOMWidgets Qt5|Sql
    # PRIVATE Qt5|OpenGL+Widgets+Xml
  # )
cmake_minimum_required(VERSION 3.0)
if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
endif()
if(POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW) # CMake 3.0
endif()
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 NEW) # CMake 3.0
endif()
if(POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW) # CMake 3.1
endif()

project(DicomUI)

set(CTK_DIR "E:/Tech318/MITK/MITK/Build/ep/src/CTK-build")
find_package(CTK REQUIRED)

include_directories(${CTK_INCLUDE_DIRS})

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

set(Hdrs
  QmitkDicomExternalDataWidget.h
  QmitkDicomLocalStorageWidget.h
)

set(Srcs
  QmitkDicomLocalStorageWidget.cpp
  QmitkDicomExternalDataWidget.cpp
)

set(UIs
  QmitkDicomLocalStorageWidgetControls.ui
  QmitkDicomExternalDataWidgetControls.ui
)

set(MOC_Hdrs
  QmitkDicomLocalStorageWidget.h
  QmitkDicomExternalDataWidget.h
)


# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets Core Gui Sql OpenGL Xml REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})

# CMAKE_AUTOMOC in ON so the MocHdrs will be automatically wrapped.
add_library(DicomUI ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs})
qt5_use_modules(DicomUI Core Gui Widgets Sql OpenGL Xml)
target_link_libraries(DicomUI CTKDICOMWidgets)

