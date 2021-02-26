# CMake generated Testfile for 
# Source directory: /home/fanli/workspace/gromacs_fh_debug/src/programs/mdrun/tests
# Build directory: /home/fanli/workspace/gromacs_fh_debug/src/programs/mdrun/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(MdrunTests "/home/fanli/workspace/gromacs_fh_debug/bin/mdrun-test" "--gtest_output=xml:/home/fanli/workspace/gromacs_fh_debug/Testing/Temporary/MdrunTests.xml")
set_tests_properties(MdrunTests PROPERTIES  LABELS "IntegrationTest" _BACKTRACE_TRIPLES "/home/fanli/workspace/gromacs_fh_debug/src/testutils/TestMacros.cmake;96;add_test;/home/fanli/workspace/gromacs_fh_debug/src/programs/mdrun/tests/CMakeLists.txt;63;gmx_register_integration_test;/home/fanli/workspace/gromacs_fh_debug/src/programs/mdrun/tests/CMakeLists.txt;0;")
add_test(MdrunMpiTests "/usr/bin/mpiexec" "-n" "2" "/home/fanli/workspace/gromacs_fh_debug/bin/mdrun-mpi-test" "--gtest_output=xml:/home/fanli/workspace/gromacs_fh_debug/Testing/Temporary/MdrunMpiTests.xml")
set_tests_properties(MdrunMpiTests PROPERTIES  LABELS "MpiIntegrationTest" _BACKTRACE_TRIPLES "/home/fanli/workspace/gromacs_fh_debug/src/testutils/TestMacros.cmake;135;add_test;/home/fanli/workspace/gromacs_fh_debug/src/programs/mdrun/tests/CMakeLists.txt;83;gmx_register_mpi_integration_test;/home/fanli/workspace/gromacs_fh_debug/src/programs/mdrun/tests/CMakeLists.txt;0;")
