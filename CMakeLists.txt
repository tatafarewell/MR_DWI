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

project(MR_diffusion)

add_subdirectory(DicomWidget)
add_subdirectory(DiffusionCore)

set(CTK_DIR "C:/Tech318/CTK/CTK/CTK-Superbuild")
find_package(CTK REQUIRED)
include(${CTK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set(CMAKE_AUTOMOC ON)

# Set your files and resources here
set( Srcs mainapp.cpp mainhub.cpp)
set( Hdrs mainhub.h)
set( MOC_Hdrs mainhub.h)
set( UIs mainhub.ui)

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  "DicomWidget"
)

find_package(Qt5 COMPONENTS Widgets REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})
qt5_add_resources(UI_qrc "resources/dicom.qrc")

# CMAKE_AUTOMOC in ON so the MocHdrs will be automatically wrapped.
add_executable(DWIAPP ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} ${UI_qrc})
qt5_use_modules(DWIAPP Core Gui Widgets Sql Xml UiTools XmlPatterns Concurrent)
target_link_libraries(DWIAPP ${QT_LIBRARIES} ${VTK_LIBRARIES} DicomWidget CTKDICOMWidgets)

# install(TARGETS DWIAPP DESTINATION bin)
