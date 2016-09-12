
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

set( Srcs cDWiFilter.cpp )

set( Hdrs 
     itkAdcMapFilter.hxx 
	 itkAdcMapFilter.h 
	 itkComputedDwiFilter.h 
	 itkComputedDwiFilter.hxx 
	 itkDwiIVIMCostFunction.h
	 itkDwiIVIMFilter.hxx
	 itkDwiIVIMFilter.h
	 cDWI.h
   )

# set( MOC_Hdrs QmitkDicomExternalDataWidget.h QmitkDicomLocalStorageWidget.h)
# set( UIs QmitkDicomExternalDataWidgetControls.ui QmitkDicomLocalStorageWidgetControls.ui)

add_library(DiffusionCore STATIC ${Srcs} ${Hdrs})
target_link_libraries(DiffusionCore ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${Glue})
# if(VTK_LIBRARIES)
  # target_link_libraries(DiffusionCore ${VTK_LIBRARIES})
# else()
  #target_link_libraries(DiffusionCore vtkHybrid vtkWidgets)
# endif()


