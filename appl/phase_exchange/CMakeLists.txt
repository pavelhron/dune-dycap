add_executable(exchange exchange.cc)
dune_symlink_to_source_files(FILES exchange.ini data_to_fit.dat exchange.gpl exchange.eps)
