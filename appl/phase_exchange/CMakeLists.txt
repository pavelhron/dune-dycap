add_executable(exchange exchange.cc)
dune_symlink_to_source_files(FILES exchange.ini data_to_fit.dat exchange.gpl exchange.eps create_figures.sh)

add_subdirectory(doc)
