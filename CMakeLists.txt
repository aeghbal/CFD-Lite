cmake_minimum_required(VERSION 3.9)
set(CMAKE_Fortran_COMPILER "/opt/rh/devtoolset-7/root/usr/bin/gfortran")
set(CMAKE_CXX_COMPILER "/opt/rh/devtoolset-7/root/usr/bin/g++")
set(CMAKE_CC_COMPILER "/opt/rh/devtoolset-7/root/usr/bin/gcc")
project(runcfd LANGUAGES Fortran CXX C)

set(USE_CATALYST OFF CACHE BOOL "Link the simulator with Catalyst")
set(USE_VTK ON CACHE BOOL "Link the simulator with VTK")
set(USE_HDF5 ON CACHE BOOL "Link the simulator with HDF5")
set(USE_CGNS ON CACHE BOOL "Link the simulator with CGNS")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_BUILD_FILES_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_BUILD_DIRECTORY ${CMAKE_BINARY_DIR}/build)
#set(CMAKE_BINARY_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CACHEFILE_DIR ${CMAKE_CACHEFILE_DIR}/cash)

set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)


set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_C_FLAGS  "${CMAKE_CXX_FLAGS} -w")
#set(CMAKE_CUDA_FLAGS  "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -w -g -cpp -fdefault-real-8 -ffree-line-length-512")

  set(Adaptor_SRCS )
  set(VTK_SRCS )

if(USE_CATALYST)
  # we only need C++ compilers if we're building the Adaptor
#  enable_language(CUDA)
  enable_language(CXX)
  enable_language(C)
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -DCatalyst ")

  find_package(ParaView 5.6  REQUIRED COMPONENTS vtkPVPythonCatalyst )
 
  include("${PARAVIEW_USE_FILE}")
  set(Adaptor_SRCS ${Adaptor_SRCS}
	${PROJECT_SOURCE_DIR}/src/modules/cfd_catalyst.cpp
	${PROJECT_SOURCE_DIR}/src/modules/cfd_adaptor.cpp
    )
  add_library(CatalystAdaptor ${Adaptor_SRCS})
  target_link_libraries(CatalystAdaptor vtkPVPythonCatalyst )
  add_definitions("-DCatalyst")

endif()

if(USE_CGNS)
  # we only need C++ compilers if we're building the CGNS
  enable_language(Fortran)
  enable_language(CXX)
  enable_language(C)

  include_directories("/opt/cgnslib_3.2.1_modified/include")
  add_library(CGNS_l STATIC IMPORTED)
  set_property(TARGET CGNS_l PROPERTY IMPORTED_LOCATION "/opt/cgnslib_3.2.1_modified/lib/libcgns.a")
  add_compile_options(-lcgns)

	if(USE_HDF5)
	  find_package(HDF5)
	  if(HDF5_FOUND)
		message("Find: HDF5 ...")
		set(HDF5_SRCS
			${PROJECT_SOURCE_DIR}/src/modules/mod_cgns.f90
			)
		add_library(HDF5_l ${HDF5_SRCS})
		target_link_libraries(HDF5_l ${HDF5_LIBRARIES})
		add_definitions("-DUSE_HDF5")
	  else()
		message("Not Find: HDF5 ...")
	  endif()
	endif()

endif()


if(USE_VTK)
#  enable_language(CUDA)
  enable_language(CXX)
  enable_language(C)
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -DVTK ")

	find_package(VTK REQUIRED COMPONENTS vtkCommonCore vtkIOXML vtkCommonDataModel)
	include(${VTK_USE_FILE})

	set(VTK_SRC
		${PROJECT_SOURCE_DIR}/src/modules/cfd_writer.cpp
	    )
	set(VTK_SRCS ${VTK_SRCS}
		${PROJECT_SOURCE_DIR}/src/modules/cfd_writer.cpp
		${PROJECT_SOURCE_DIR}/src/modules/vtk_data.f90
	    )
	add_library(VTK_l ${VTK_SRC})
	target_link_libraries(VTK_l ${VTK_LIBRARIES})
  add_definitions("-DVTK")
endif()

add_executable(runcfd 	${PROJECT_SOURCE_DIR}/src/main.f90
							${PROJECT_SOURCE_DIR}/src/cfile.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_cgns.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_subdomains.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_eqn_setup.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_physics.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_solver.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_util.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_properties.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_multiphase.f90
							${PROJECT_SOURCE_DIR}/src/modules/mod_vtu_output.f90
							${PROJECT_SOURCE_DIR}/src/equations/mod_energy.f90
							${PROJECT_SOURCE_DIR}/src/equations/mod_mfr.f90
							${PROJECT_SOURCE_DIR}/src/equations/mod_scalar.f90
							${PROJECT_SOURCE_DIR}/src/equations/mod_uvwp.f90
							${PROJECT_SOURCE_DIR}/src/equations/mod_vfr.f90
							${PROJECT_SOURCE_DIR}/src/setup/calc_aip_xyzip.f90
							${PROJECT_SOURCE_DIR}/src/setup/calc_vol_cv_centers.f90
							${PROJECT_SOURCE_DIR}/src/setup/cell_input.f90
							${PROJECT_SOURCE_DIR}/src/setup/mod_agglomeration.f90
							${PROJECT_SOURCE_DIR}/src/setup/mod_mg_lvl_uns.f90
							${PROJECT_SOURCE_DIR}/src/setup/mod_meshds_uns.f90

							${Adaptor_SRCS}
							${VTK_SRCS}
				)

set_property(TARGET runcfd PROPERTY LINKER_LANGUAGE Fortran)

if(USE_CGNS AND USE_HDF5 AND USE_VTK AND USE_CATALYST)
  target_link_libraries(runcfd 		LINK_PRIVATE CGNS_l
												LINK_PRIVATE HDF5_l
												LINK_PRIVATE VTK_l
												LINK_PRIVATE CatalystAdaptor
						)
  include(vtkModuleMacros)
elseif(USE_CGNS AND USE_HDF5 AND USE_VTK)
  target_link_libraries(runcfd 		LINK_PRIVATE CGNS_l
												LINK_PRIVATE HDF5_l
												LINK_PRIVATE VTK_l
						)
  include(vtkModuleMacros)

elseif(USE_CGNS AND USE_HDF5)
  target_link_libraries(runcfd 		LINK_PRIVATE CGNS_l
												LINK_PRIVATE HDF5_l
						)
endif()

option(BUILD_TESTING "Build Testing" OFF)
# Setup testing.
if (BUILD_TESTING)
  include(CTest)
  add_test(NAME runTest COMMAND runcfd ${CMAKE_CURRENT_SOURCE_DIR}/SampleScripts/feslicescript.py)
  set_tests_properties(runTest PROPERTIES LABELS "CFD_PV;Simple")
endif()
