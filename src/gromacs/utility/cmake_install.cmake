# Install script for directory: /home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local/gromacs")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gromacs/utility" TYPE FILE FILES
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/arrayref.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/arraysize.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/basedefinitions.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/baseversion.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/classhelpers.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/cstringutil.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/current_function.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/datafilefinder.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/errorcodes.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/exceptions.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/fatalerror.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/flags.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/futil.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/gmxassert.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/init.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/programcontext.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/real.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/smalloc.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/stringutil.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/tests/cmake_install.cmake")

endif()

