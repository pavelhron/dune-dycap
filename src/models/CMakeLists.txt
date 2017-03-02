install(FILES adroperator.hh
              componenttransportop.hh
	      componenttransportop_nocflcontrolle.hh
	      convectiondiffusiondg.hh
	      defaultimp.hh
	      fluxreconstruction.hh
	      limiter.hh
	      multicomponenttransportop.hh
	      pcpl_twophaseccfv.hh
	      pcpl_twophaseccfv_incompressible.hh
	      pcpl_twophaseccfv_totalflux.hh
	      plsg_twophaseccfv.hh
	      reaction_ode_models.hh
	      reactionop.hh
	      reactiontransportop.hh
	      transporttraits.hh
	      twophaseccfv.hh
	      twophasetraits.hh
	      twophasevelocity.hh
	      twosideimp.hh
	      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/src/models)
