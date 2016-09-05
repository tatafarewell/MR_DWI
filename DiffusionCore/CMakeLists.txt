
cmake_minimum_required(VERSION 2.8)

PROJECT(cDWI)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(cDWI MACOSX_BUNDLE cDWI)

if(VTK_LIBRARIES)
  target_link_libraries(cDWI ${VTK_LIBRARIES})
else()
  target_link_libraries(cDWI vtkHybrid vtkWidgets)
endif()
target_link_libraries(cDWI ${ITK_LIBRARIES} ${Glue})
