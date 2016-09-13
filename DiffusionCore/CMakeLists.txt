
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

# find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})
# find_package(ITK REQUIRED)
# include(${ITK_USE_FILE})

set( Srcs 
    cDWiFilter.cpp 
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
	 cDWI.h
     DicomHelper.h
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
qt5_use_modules(DicomWidget Core Gui Widgets opengl)

add_library(DiffusionCore STATIC ${Srcs} ${Hdrs} ${MOC_Hdrs} ${UIs}})
target_link_libraries(DiffusionCore ${QT_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${Glue})
# if(VTK_LIBRARIES)
  # target_link_libraries(DiffusionCore ${VTK_LIBRARIES})
# else()
  #target_link_libraries(DiffusionCore vtkHybrid vtkWidgets)
# endif()


