# Install script for directory: /home/fanli/workspace/gromacs_fh_debug/src/gromacs

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gromacs" TYPE FILE FILES
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/analysisdata.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/commandline.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/options.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/random.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/selection.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/trajectoryanalysis.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility.h"
    "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/version.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/fanli/workspace/gromacs_fh_debug/lib/libgromacs_mpi.so.2.0.0"
    "/home/fanli/workspace/gromacs_fh_debug/lib/libgromacs_mpi.so.2"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgromacs_mpi.so.2.0.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgromacs_mpi.so.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/fanli/workspace/gromacs_fh_debug/lib/libgromacs_mpi.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgromacs_mpi.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgromacs_mpi.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgromacs_mpi.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE RENAME "libgromacs_mpi.pc" FILES "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/libgromacs.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/gromacs_mpi/libgromacs.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/gromacs_mpi/libgromacs.cmake"
         "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/CMakeFiles/Export/share/cmake/gromacs_mpi/libgromacs.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/gromacs_mpi/libgromacs-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/gromacs_mpi/libgromacs.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/gromacs_mpi" TYPE FILE FILES "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/CMakeFiles/Export/share/cmake/gromacs_mpi/libgromacs.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/gromacs_mpi" TYPE FILE FILES "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/CMakeFiles/Export/share/cmake/gromacs_mpi/libgromacs-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/gromacs_mpi" TYPE FILE RENAME "gromacs_mpi-config.cmake" FILES "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/gromacs-config.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/gromacs_mpi" TYPE FILE RENAME "gromacs_mpi-config-version.cmake" FILES "/home/fanli/workspace/gromacs_fh_debug/src/gromacs/gromacs-config-version.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/gmxlib/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/mdlib/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/listed-forces/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/commandline/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/domdec/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/ewald/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/gpu_utils/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/hardware/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/linearalgebra/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/math/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/mdrunutility/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/mdtypes/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/onlinehelp/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/options/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/pbcutil/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/random/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/tables/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/timing/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/topology/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/trajectory/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/utility/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/fileio/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/swap/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/essentialdynamics/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/pulling/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/imd/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/gmxana/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/gmxpreprocess/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/correlationfunctions/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/statistics/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/analysisdata/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/selection/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/trajectoryanalysis/cmake_install.cmake")
  include("/home/fanli/workspace/gromacs_fh_debug/src/gromacs/tools/cmake_install.cmake")

endif()

