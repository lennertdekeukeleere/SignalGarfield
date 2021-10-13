#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(SignalGarfield)
set(CMAKE_C_FLAGS "-std=c++14")
#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT QUIET REQUIRED COMPONENTS RIO Net)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT package not found. --> program disabled")  
  return()
endif()

FIND_PACKAGE(Boost)

include_directories(${BOOST_INCLUDE_DIRS})

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

find_package(Garfield)
if(NOT Garfield_FOUND)
  MESSAGE(STATUS "Could not find cmake configuration for packgage Garfield. Assuming it was build using the make file.")
else()
  MESSAGE(STATUS "Found package Garfield in ${GARFIELD_INCLUDE_DIRS}")
endif()
# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)
message(${ROOT_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS} ${Geant4_INCLUDE_DIR})
# Setup GARFIELD++ include directories
if(NOT Garfield_FOUND)
include_directories($ENV{GARFIELD_HOME}/Include $ENV{GARFIELD_HOME}/Heed)
link_directories($ENV{GARFIELD_HOME}/Library)
else()
include_directories(${GARFIELD_INCLUDE_DIRS})
endif()
#---------------------------------------------------------------------------
#include the macro for generating dictionaries


#--------------------------------------------------------------------------
#Create a shared library with the generated dictionary



#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
#file(GLOB headers $ENV{GARFIELD_HOME}/include/*.hh)
#file(GLOB headers $ENV{GARFIELD_HOME}/heed/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 and GARFIELD++ libraries
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g")
add_executable(SignalGarfield SignalGarfield.cc ${sources} ${headers})
if(NOT Garfield_FOUND)
target_link_libraries(SignalGarfield -lGarfield -lgfortran ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${Boost_LIBRARIES} /usr/lib64/libboost_program_options.so.1.66.0)
else()
target_link_libraries(SignalGarfield ${GARFIELD_LIBRARIES} -lgfortran ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${Boost_LIBRARIES} /usr/lib64/libboost_program_options.so.1.66.0)
endif()
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(SCRIPTS
  input.txt
  config.txt
)

foreach(_script ${SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SignalGarfield DESTINATION bin)
