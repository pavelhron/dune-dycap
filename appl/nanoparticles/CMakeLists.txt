add_executable(estimation estimation.cc)
dune_symlink_to_source_files(FILES Readme.txt tracer.conf)


add_subdirectory("data")
add_subdirectory("parameters")
add_subdirectory("utilities")
