add_executable(simulation simulation_boilerplate.cc)
dune_symlink_to_source_files(FILES ecolimodel.conf ecoli_dfz.conf ecoli_dfz_inflow.conf)

# create vtk dir, as otherwise the executable just falls over
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/vtk")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()
