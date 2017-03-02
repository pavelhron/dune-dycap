# module providing convenience methods for compiling binaries with gsl support
#
function(add_dune_gsl_flags _targets)
  if(GSL_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} ${GSL_LIBRARY} -lgsl -lgslcblas -lm)
      set_property(TARGET ${_target} APPEND PROPERTY INCLUDE_DIRECTORIES ${GSL_INCLUDE_DIRS} )
      set_property(TARGET ${_target} APPEND_STRING PROPERTY COMPILE_FLAGS "-lgsl -lgsllcblas")
    endforeach(_target ${_targets})
  endif(GSL_FOUND)
endfunction(add_dune_gsl_flags)
