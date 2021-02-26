# CMake generated Testfile for 
# Source directory: /home/fanli/workspace/gromacs_fh_debug/src/gromacs/gmxana/legacytests
# Build directory: /home/fanli/workspace/gromacs_fh_debug/src/gromacs/gmxana/legacytests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(LegacyToolsTests "/home/fanli/workspace/gromacs_fh_debug/bin/legacy-tools-test" "--gtest_output=xml:/home/fanli/workspace/gromacs_fh_debug/Testing/Temporary/LegacyToolsTests.xml")
set_tests_properties(LegacyToolsTests PROPERTIES  LABELS "IntegrationTest" _BACKTRACE_TRIPLES "/home/fanli/workspace/gromacs_fh_debug/src/testutils/TestMacros.cmake;96;add_test;/home/fanli/workspace/gromacs_fh_debug/src/gromacs/gmxana/legacytests/CMakeLists.txt;43;gmx_register_integration_test;/home/fanli/workspace/gromacs_fh_debug/src/gromacs/gmxana/legacytests/CMakeLists.txt;0;")
