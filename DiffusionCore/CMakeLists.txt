
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

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
find_package(CTK REQUIRED)
include(${CTK_USE_FILE})

set( Srcs 
    DicomHelper.cpp
    DiffusionCore.cpp
    )

set( Hdrs 
     itkAdcMapFilter.hxx 
     itkAdcMapFilter.h 
     itkComputedDwiFilter.h 
     itkComputedDwiFilter.hxx 
     itkDwiIVIMCostFunction.h
     itkDwiIVIMFilter.hxx
     itkDwiIVIMFilter.h
     itkComputedEadcFilter.h
     itkComputedEadcFilter.hxx
     itkDwiIVIMCostFunction2.h
     itkDwiIVIMFilter2.hxx
     itkDwiIVIMFilter2.h
     itkDisplayOptimizer.h
     itkDisplayOptimizer.hxx
     itkMaskVectorImageFilter.h
     itkMaskVectorImageFilter.hxx
     itkTensor.h
     itkTensor.hxx
     itkGetDiffusionImageFilter.h
     itkGetDiffusionImageFilter.hxx
     itkVtkAllHeaders.h
     DicomHelper.h
     DiffusionCore.h
   )

set( MOC_Hdrs DiffusionCore.h)
set( UIs DiffusionCore.ui)

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)


find_package(Qt5 COMPONENTS Widgets REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})
add_library(DiffusionCore STATIC ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} )
qt5_use_modules(DiffusionCore Core Gui Widgets Sql Xml UiTools XmlPatterns Concurrent)
target_link_libraries(DiffusionCore ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${Glue} CTKWidgets)
# if(VTK_LIBRARIES)
  # target_link_libraries(DiffusionCore ${VTK_LIBRARIES})
# else()
  # target_link_libraries(DiffusionCore vtkHybrid vtkWidgets)
# endif()


