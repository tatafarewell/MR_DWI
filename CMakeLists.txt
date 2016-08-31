cmake_minimum_required(VERSION 3.0)

project(MR_DIFFUSION)

set(DICOMUI )
include(${CTK_USE_FILE})

find_package(Qt4 REQUIRED)
include(${QT_USE_FILE})

add_executable(example1 main.cpp)
target_link_libraries(example1 CTKWidgets)