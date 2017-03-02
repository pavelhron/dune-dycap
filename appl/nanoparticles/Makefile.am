SUBDIRS =
noinst_PROGRAMS = estimation

estimation_SOURCES = estimation.cc
estimation_CPPFLAGS = $(AM_CPPFLAGS)\
	$(SUPERLU_CPPFLAGS)\
	$(DUNEMPICPPFLAGS)\
	$(PARMETIS_CPPFLAGS)
estimation_LDADD = $(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) \
	$(DUNE_LIBS) \
	$(DUNEMPILIBS) \
	$(LDADD) \
	$(PARMETIS_LIBS)
estimation_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(DUNE_LDFLAGS) \
	$(PARMETIS_LDFLAGS)


# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

CLEANFILES = *.pvd ./vtk/* ./cut/* log ./boundary/* *~

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
