# CMake generated Testfile for 
# Source directory: /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fileio/tests
# Build directory: /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fileio/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(FileIOTests "/home/fanli/workspace/gromacs_fh_debug/bin/fileio-test" "--gtest_output=xml:/home/fanli/workspace/gromacs_fh_debug/Testing/Temporary/FileIOTests.xml")
set_tests_properties(FileIOTests PROPERTIES  LABELS "GTest;UnitTest" _BACKTRACE_TRIPLES "/home/fanli/workspace/gromacs_fh_debug/src/testutils/TestMacros.cmake;86;add_test;/home/fanli/workspace/gromacs_fh_debug/src/testutils/TestMacros.cmake;166;gmx_register_unit_test;/home/fanli/workspace/gromacs_fh_debug/src/gromacs/fileio/tests/CMakeLists.txt;42;gmx_add_unit_test;/home/fanli/workspace/gromacs_fh_debug/src/gromacs/fileio/tests/CMakeLists.txt;0;")
