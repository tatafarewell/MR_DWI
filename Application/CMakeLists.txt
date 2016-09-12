project(DWIAPP)

set(CTK_DIR "C:/Tech318/CTK/CTK/CTK-Superbuild")
find_package(CTK REQUIRED)
include(${CTK_USE_FILE})

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)


# Set your files and resources here
set( Srcs mainapp.cpp mainhub.cpp)
set( Hdrs mainhub.h)
set( MOC_Hdrs mainhub.h)
set( UIs mainhub.ui)

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

find_package(Qt5 COMPONENTS Widgets REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})

# CMAKE_AUTOMOC in ON so the MocHdrs will be automatically wrapped.
add_library(DWIAPP SHARED ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs})
qt5_use_modules(DWIAPP Core Gui Widgets Sql Xml UiTools XmlPatterns Concurrent)
target_link_libraries(DWIAPP ${QT_LIBRARIES} DicomWidget CTKDICOMWidgets)

install(TARGETS DWIAPP DESTINATION bin)