cmake_minimum_required(VERSION 2.8.8)
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
project(MR_Diffusion)

add_subdirectory(DicomWidget)
add_subdirectory(DiffusionCore)

find_package(VTK REQUIRED)
find_package(CTK REQUIRED)
find_package(ITK REQUIRED)

include(${VTK_USE_FILE})
include(${CTK_USE_FILE})
include(${ITK_USE_FILE})

if("${VTK_QT_VERSION}" STREQUAL "")
  message(FATAL_ERROR "VTK was not built with Qt")
endif()

# Set your files and resources here
set( Srcs main.cpp mainwindow.cpp)

set( Hdrs mainwindow.h )

set( MOC_Hdrs mainwindow.h )

set( UIs mainwindow.ui )

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  "DicomWidget"
  "DiffusionCore"
)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

if(VTK_QT_VERSION VERSION_GREATER "4")
  find_package(Qt5 COMPONENTS Widgets REQUIRED QUIET)
  qt5_wrap_ui(UI_Srcs ${UIs})

  # CMAKE_AUTOMOC in ON so the MocHdrs will be automatically wrapped.
  add_executable(MR_Diffusion ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs})
  qt5_use_modules(MR_Diffusion Core Gui Widgets Sql Xml UiTools XmlPatterns Concurrent)
  target_link_libraries(MR_Diffusion ${VTK_LIBRARIES} ${QT_LIBRARIES}  DicomWidget DiffusionCore CTKWidgets CTKDICOMWidgets)
else()
  find_package(Qt4 REQUIRED)
  include(${QT_USE_FILE})
  # Use what VTK built with
  set(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
  set(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
  set(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
  qt4_wrap_ui(UI_Srcs ${UIs})

  add_executable(MR_Diffusion ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs})
  target_link_libraries(MR_Diffusion ${QT_LIBRARIES} ${VTK_LIBRARIES})
endif()
