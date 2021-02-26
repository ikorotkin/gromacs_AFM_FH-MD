#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libgromacs" for configuration "Debug"
set_property(TARGET libgromacs APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(libgromacs PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so;dl;rt;m;gmxfftw;-lpthread;-fopenmp"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libgromacs_mpi.so.2.0.0"
  IMPORTED_SONAME_DEBUG "libgromacs_mpi.so.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS libgromacs )
list(APPEND _IMPORT_CHECK_FILES_FOR_libgromacs "${_IMPORT_PREFIX}/lib/libgromacs_mpi.so.2.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
