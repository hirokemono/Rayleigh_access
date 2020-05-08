#! /bin/csh
#
#
TARGET_PICK_RAYLEIGH_SPH =        pickup_rayleigh_spectr

MOD_PICK_RAYLEIGH_SPH = \
main_pickup_rayleigh_spectr.o \
\
m_precision.o \
calypso_mpi_4_rayleigh.o \
skip_comment_f.o \
t_field_data_IO.o \
t_read_control_elements.o \
t_control_array_integer2.o \
\
t_ctl_pick_rayleigh_spectr.o \
t_picked_rayleigh_parameter.o \
t_picked_rayleigh_spectr.o \
\
t_rayleigh_restart_IO.o \
rayleigh_restart_IO.o \
rayleigh99_rst_param_IO.o \
MPI_read_rayleigh_restart.o


# -------  do not edit folloing lines
#
# MPI settings
#
MPICHDIR = /usr/local/mpich
MPILIBDIR = $(MPICHDIR)/lib
MPILIBS = -mpich -fmpich
MPIBINDIR = $(MPICHDIR)/bin
MPIINCDIR = $(MPICHDIR)/include
#
#
MAKEDIR= 
BUILDDIR= 
SRCDIR = $(SRCDIR)

MAKENAME= $(MAKEDIR)/Makefile

F_INCLUDE = -I.

ifdef MPILIBS
  F90LIBS+= -L$(MPILIBDIR) $(MPILIBS)
endif

DEPFLAGS= -MM -w -DDEPENDENCY_CHECK

#
#  -------------------------------------------------------------------------
#

.SUFFIXES: .o .f90

.f90.o:
	$(F90) -o $@ @<



all: fluid

fluid:
	echo "# Compiling start"; \
		cd work;  \
		$(GMAKE) \
		GMAKE='$(GMAKE)' \
		fluid

makemake:
	if [ ! -d $(MAKEDIR) ]; then \
		mkdir $(MAKEDIR); \
	fi
	if [ ! -d $(BUILDDIR) ]; then \
		mkdir $(BUILDDIR); \
	fi
	@echo # create module list
	echo '#' > $(MAKENAME)
	echo '# Makefile for Kemono Dynamo simulation' >> $(MAKENAME)
	echo '#' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo 'SHELL = /bin/sh' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '#  all:         compile programs for simulations' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '# Compiler commands' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo F90=$(MPIF90)          >> $(MAKENAME)
	@echo AR=$(AR)               >> $(MAKENAME)
	@echo RANLIB=$(RANLIB) >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '# Compiler options' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo 'F90OPTFLAGS = $(F90OPTFLAGS)  $$(F_INCLUDE)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '# -----------------------------------------------------------' >> $(MAKENAME)
	@echo '# -----------------------------------------------------------' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '# Directory lists' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo MAKEDIR =  $(MAKEDIR)  >> $(MAKENAME)
	@echo BUILDDIR = $(BUILDDIR) >> $(MAKENAME)
	@echo SRCDIR =   $(SRCDIR)   >> $(MAKENAME)
	@echo INSTDIR =  $(INSTDIR)  >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	>> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '# -----------------------------------------------------------' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '# target lists ' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo 'PROGRAMS= rayleigh_link' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo 'create excecute file list...'
		$(GMAKE) \
		BUILDDIR='$(BUILDDIR)' \
		MAKENAME='$(MAKENAME)' \
		target_list
	@echo '#' >> $(MAKENAME)
	@echo '# -----------------------------------------------------------' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '# objects list' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo 'create object module list...'
		$(GMAKE) \
		MAKENAME='$(MAKENAME)' \
		mod_list
	@echo '' >> $(MAKENAME)
	@echo 'create excecution construct list...'
	@echo '#' >> $(MAKENAME)
	@echo 'all: $$(PROGRAMS)' >> $(MAKENAME)
		$(GMAKE) \
		GMAKE='$(GMAKE)' \
		MAKENAME='$(MAKENAME)' \
		target_task
	@echo '' >> $(MAKENAME)
	@echo 'create objects construct list...'
	@echo '#' >> $(MAKENAME)
	@echo '# Rules for routines' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
		$(GMAKE) \
		GMAKE='$(GMAKE)' \
		MAKENAME='$(MAKENAME)' \
		module
	@echo '#' >> $(MAKENAME)
	@echo '# Rules for install' >> $(MAKENAME)
	@echo '#' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo 'install:' >> $(MAKENAME)
	@echo '	if [ ! -d $$(INSTDIR) ]; then \'       >> $(MAKENAME)
	@echo '		mkdir $$(INSTDIR); \'              >> $(MAKENAME)
	@echo '	fi'                                   >> $(MAKENAME)
	@echo '	if [ ! -d $$(INSTDIR)/bin ]; then \'     >> $(MAKENAME)
	@echo '		mkdir $$(INSTDIR)/bin; \'            >> $(MAKENAME)
	@echo '	fi'                                     >> $(MAKENAME)
	@echo '	cp $(BUILDDIR)/*  $$(INSTDIR)/bin' >> $(MAKENAME)
	@echo ''       >> $(MAKENAME)
	@echo 'clean:' >> $(MAKENAME)
	@echo '	rm -f *.o ;'       >> $(MAKENAME)
	@echo '	rm -f *.mod ;'     >> $(MAKENAME)
	@echo '	rm -f *~ ;'        >> $(MAKENAME)
	@echo '	rm -f *.par ;'     >> $(MAKENAME)
	@echo '	rm -f *.diag ;'    >> $(MAKENAME)
	@echo '	rm -f *.a *.L'     >> $(MAKENAME)
	@echo '' >> $(MAKENAME)

install:
	if [ ! -d $$(INSTDIR) ]; then \
		mkdir $$(INSTDIR); \
	fi
	if [ ! -d $$(INSTDIR)/lib ]; then \
		mkdir $$(INSTDIR)/lib ; \
	fi
	if [ ! -d $$(INSTDIR)/bin ]; then \
		mkdir $$(INSTDIR)/bin; \
	fi
	if [ ! -d $$(INSTDIR)/include ]; then \
		mkdir $$(INSTDIR)/include; \
	fi
	@for dir in geoutil solver MHD surface marge_fluid; do \
	 (cd $${dir}; \
		$(GMAKE) \
		GMAKE='$(GMAKE)' \
		INSTDIR=$$(INSTDIR) RANLIB=$(RANLIB) \
		install; \
		) \
	done

target_list:
	@echo 'TARGET_PICK_RAYLEIGH_SPH = $$(BUILDDIR)/$(TARGET_PICK_RAYLEIGH_SPH)' >> $(MAKENAME)
	@echo >> $(MAKENAME)

target_task:
	@echo rayleigh_link: \
	'$$(TARGET_PICK_RAYLEIGH_SPH)' \
	>> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_PICK_RAYLEIGH_SPH)': '$$(MOD_PICK_RAYLEIGH_SPH)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90OPTFLAGS)' -o '$$(TARGET_PICK_RAYLEIGH_SPH)' '$$(MOD_PICK_RAYLEIGH_SPH)' \
	 >> $(MAKENAME)
	@echo '' >> $(MAKENAME)


mod_list:
	@echo  MOD_PICK_RAYLEIGH_SPH=  \\ >> $(MAKENAME)
	@echo  $(MOD_PICK_RAYLEIGH_SPH)   >> $(MAKENAME)
	@echo  '#' >> $(MAKENAME)


module:
	@cat Makefile.depends >> $(MAKENAME)

clean:
	rm -f *.o *~ *.mod *.diag *.L
