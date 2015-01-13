cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(SingleBrillanceCyrstal)

#if user does not specify prefix we assign it to the exec directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	message(STATUS "Install Prefix not specified.")
	file(MAKE_DIRECTORY exec)
	get_filename_component(INSTALL_DIR ${CMAKE_SOURCE_DIR}/exec REALPATH)
	set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH
			"Install Prefix" FORCE)
endif ()
message(STATUS "Installing to ${CMAKE_INSTALL_PREFIX}")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by 
# default. You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui to build a batch mode only executable.
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#--------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

install(DIRECTORY share DESTINATION ${CMAKE_INSTALL_PREFIX})

add_subdirectory(src)
