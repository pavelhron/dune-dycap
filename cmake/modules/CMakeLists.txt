set(modules
  AddGSLFlags.cmake
  DuneDycapMacros.cmake
  FindGSL.cmake)


install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
