#
# NAVIER-STOKES FINITE-DIFFERENCES SOLVER MAKE FILE
#
# ######################################################################
#
#brief      MAKEFILE
#
#history    Sergio Castiblanco
#+          13/01/2021
########################################################################
#
# SET COMPILER OPTIONS
#
FC = gfortran
#
# COMPILER FLAGS
#
FCFLAGS = -m64 -O2
#
# MOD FLAG FOR .mod FILES LOCATION
#
MODFLAG = -J
#
# MAIN DIRECTORY
#
NSDIR = $(CURDIR)
#
# DIRECTORIES FOR SOURCE, OBJECTS, MODULES AND BINARIES
#
SRCDIR = $(NSDIR)/src
OBJDIR = $(NSDIR)/obj
MODDIR = $(NSDIR)/mod
BINDIR = $(NSDIR)/bin
#
# COMPILER AND LINKER FLAGS - GNU FORTRAN
#
CFLAGS = -c -g -Warray-bounds -fbacktrace -fbounds-check -Wall $(FCFLAGS)
LFLAGS = -g -Warray-bounds -fbacktrace -fbounds-check -Wall $(FCFLAGS)
#
# FLAGS FOR MODULE DIRECTORY
#
CFLAGS += $(MODFLAG) $(MODDIR)
LFLAGS += $(MODFLAG) $(MODDIR)
#########################################################################
#
# SOURCE FILES
#
SRCMAIN = $(SRCDIR)/home_diff.f
SRCNUME = $(SRCDIR)/declarations_numerical.f
SRCPHYS = $(SRCDIR)/declarations_physic.f
SRCCSCS = $(SRCDIR)/csc_storage.f
SRCDCSC = $(SRCDIR)/declarations_csc.f
SRCWRHE = $(SRCDIR)/write_headers.f
SRCPONS = $(SRCDIR)/point_diff.f
SRCALLC = $(SRCDIR)/all_csc.f
SRCDIAG = $(SRCDIR)/csc_diag.f
SRCKRON = $(SRCDIR)/csc_kron.f
SRCSUMC = $(SRCDIR)/csc_sum.f
SRCDIFF = $(SRCDIR)/diffusion_matrix.f
#
# OBJECT FILES
#
OBJMAIN = $(OBJDIR)/home_diff.o
OBJNUME = $(OBJDIR)/declarations_numerical.o
OBJPHYS = $(OBJDIR)/declarations_physic.o
OBJCSCS = $(OBJDIR)/csc_storage.o
OBJDCSC = $(OBJDIR)/declarations_csc.o
OBJWRHE = $(OBJDIR)/write_headers.o
OBJPONS = $(OBJDIR)/point_diff.o
OBJALLC = $(OBJDIR)/all_csc.o
OBJDIAG = $(OBJDIR)/csc_diag.o
OBJKRON = $(OBJDIR)/csc_kron.o
OBJSUMC = $(OBJDIR)/csc_sum.o
OBJDIFF = $(OBJDIR)/diffusion_matrix.o
#
OBJECTS = $(OBJDCSC) $(OBJCSCS) $(OBJPHYS) $(OBJNUME) $(OBJWRHE) \
			   	$(OBJPONS) $(OBJALLC) $(OBJDIAG) $(OBJKRON) $(OBJSUMC) \
				 	$(OBJDIFF) $(OBJMAIN)
#
#
# MODULE FILES
#
MODNUME = $(MODDIR)/declarations_numerical.mod
MODPHYS = $(MODDIR)/declarations_physic.mod
MODDCSC = $(MODDIR)/declarations_csc.mod
MODCSCS = $(MODDIR)/csc_storage.mod
#
MODULES = $(MODDCSC) $(MODCSCS) $(MODPHYS) $(MODNUME)
#
#
# BINARY FILES
#
BINNSFD = $(BINDIR)/diff_df.out
#########################################################################
#
# !!!COMPILING!!!
#
$(BINNSFD): $(OBJECTS)
	    $(FC) $(LFLAGS) $(OBJECTS) -o $(BINNSFD)

$(OBJDCSC): $(SRCDCSC)
	    $(FC) $(CFLAGS) $(SRCDCSC) -o $(OBJDCSC)

$(OBJCSCS): $(SRCCSCS)
	    $(FC) $(CFLAGS) $(SRCCSCS) -o $(OBJCSCS)

$(OBJPHYS): $(SRCPHYS)
	    $(FC) $(CFLAGS) $(SRCPHYS) -o $(OBJPHYS)

$(OBJNUME): $(SRCNUME)
	    $(FC) $(CFLAGS) $(SRCNUME) -o $(OBJNUME)

$(OBJWRHE): $(SRCWRHE)
	    $(FC) $(CFLAGS) $(SRCWRHE) -o $(OBJWRHE)

$(OBJPONS): $(SRCPONS)
	    $(FC) $(CFLAGS) $(SRCPONS) -o $(OBJPONS)

$(OBJALLC): $(SRCALLC)
	    $(FC) $(CFLAGS) $(SRCALLC) -o $(OBJALLC)

$(OBJDIAG): $(SRCDIAG)
	    $(FC) $(CFLAGS) $(SRCDIAG) -o $(OBJDIAG)

$(OBJKRON): $(SRCKRON)
	    $(FC) $(CFLAGS) $(SRCKRON) -o $(OBJKRON)

$(OBJSUMC): $(SRCSUMC)
	    $(FC) $(CFLAGS) $(SRCSUMC) -o $(OBJSUMC)

$(OBJDIFF): $(SRCDIFF)
	    $(FC) $(CFLAGS) $(SRCDIFF) -o $(OBJDIFF)

$(OBJMAIN): $(MODULES) $(SRCMAIN)
	    $(FC) $(CFLAGS) $(SRCMAIN) -o $(OBJMAIN)

clean:
	    rm -f $(OBJECTS)
			rm -f $(MODULES)
			rm -f $(BINNSFD)
#END



