# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if("/home/fanli/workspace/gromacs_fh_debug/src/contrib/fftw/fftw.tar.gz" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if(NOT EXISTS "/home/fanli/workspace/gromacs_fh_debug/src/contrib/fftw/fftw.tar.gz")
  message(FATAL_ERROR "File not found: /home/fanli/workspace/gromacs_fh_debug/src/contrib/fftw/fftw.tar.gz")
endif()

if("MD5" STREQUAL "")
  message(WARNING "File will not be verified since no URL_HASH specified")
  return()
endif()

if("6cc08a3b9c7ee06fdd5b9eb02e06f569" STREQUAL "")
  message(FATAL_ERROR "EXPECT_VALUE can't be empty")
endif()

message(STATUS "verifying file...
     file='/home/fanli/workspace/gromacs_fh_debug/src/contrib/fftw/fftw.tar.gz'")

file("MD5" "/home/fanli/workspace/gromacs_fh_debug/src/contrib/fftw/fftw.tar.gz" actual_value)

if(NOT "${actual_value}" STREQUAL "6cc08a3b9c7ee06fdd5b9eb02e06f569")
  message(FATAL_ERROR "error: MD5 hash of
  /home/fanli/workspace/gromacs_fh_debug/src/contrib/fftw/fftw.tar.gz
does not match expected value
  expected: '6cc08a3b9c7ee06fdd5b9eb02e06f569'
    actual: '${actual_value}'
")
endif()

message(STATUS "verifying file... done")
