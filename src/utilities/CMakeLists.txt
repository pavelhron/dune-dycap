install(FILES adhesion.hh
              cfltransportcontroller.hh
	      concentration_utilities.hh
	      controller.hh
	      gridfunction_utilities.hh
	      gridinfo.hh
	      inlets_utilities.hh
	      inlets_utilities_old.hh
	      interpolate.hh
	      multispan.hh
	      multispantwophase.hh
	      parser.hh
	      residuum.hh
	      rownorm_preconditioner.hh
	      rt0qfem.hh
	      timer.hh
	      timer_step_manager.hh
	      imex.hh
	      onestepmethods.hh
	      utilities.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/src/utilities)
