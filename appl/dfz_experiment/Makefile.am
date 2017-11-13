SUBDIRS =

noinst_PROGRAMS = simulation_boilerplate simulation_boilerplate_second


simulation_SOURCES = simulation.cc
simulation_CPPFLAGS = $(AM_CPPFLAGS)\
	$(SUPERLU_CPPFLAGS)\
	$(DUNEMPICPPFLAGS)\
	$(PARMETIS_CPPFLAGS)
simulation_LDADD = $(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) \
	$(DUNE_LIBS) \
	$(DUNEMPILIBS) \
	$(LDADD) \
	$(PARMETIS_LIBS)
simulation_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(DUNE_LDFLAGS) \
	$(PARMETIS_LDFLAGS)

simulation_boilerplate_SOURCES = simulation_boilerplate.cc
simulation_boilerplate_CPPFLAGS = $(AM_CPPFLAGS)\
	$(SUPERLU_CPPFLAGS)\
	$(DUNEMPICPPFLAGS)\
	$(PARMETIS_CPPFLAGS)\
	'-D SECONDORDER=false'
simulation_boilerplate_LDADD = $(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) \
	$(DUNE_LIBS) \
	$(DUNEMPILIBS) \
	$(LDADD) \
	$(PARMETIS_LIBS)
simulation_boilerplate_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(DUNE_LDFLAGS) \
	$(PARMETIS_LDFLAGS)

simulation_boilerplate_second_SOURCES = simulation_boilerplate.cc
simulation_boilerplate_second_CPPFLAGS = $(AM_CPPFLAGS)\
	$(SUPERLU_CPPFLAGS)\
	$(DUNEMPICPPFLAGS)\
	$(PARMETIS_CPPFLAGS)\
	'-D SECONDORDER=true'
simulation_boilerplate_second_LDADD = $(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) \
	$(DUNE_LIBS) \
	$(DUNEMPILIBS) \
	$(LDADD) \
	$(PARMETIS_LIBS)
simulation_boilerplate_second_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(DUNE_LDFLAGS) \
	$(PARMETIS_LDFLAGS)



# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

CLEANFILES = *.pvd *.txt ./vtk/* ./cut/* log* log *.pdf *.eps ./boundary/*

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
