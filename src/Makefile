#
# Makefile for GRANADA 2020
#
# The following variables will be set when the Makefile is invoked with
# a specific machine type.  The variables correspond to:
# CPP     => The C preprocessor on the machine
# DEBUG   => The flag used to invoke compiling in debug mode
# EXE     => The name of the executable that will be generated
# FC      => The name of the Fortran compiler
# FFLAGS  => A list of compile options
# LDFLAGS => A list of link or load options
# EXTRAS  => 'Special' routines used on this machine
#
#CPP	= cpp -P -traditional
DEBUG	= 
EXE     = granada4.5.x
FC	= gfortran
#FFLAGS  = -O3 -march=pentium -fbounds-check -fno-automatic
FFLAGS  = -O3 -march=native -fbounds-check -fno-automatic -fno-align-commons
#LDFLAGS	= second.o
LDFLAGS = -static 
#LDFLAGS =
MAKE	= make
EXTRAS  =
#
# This variable is set to the list of "include" files which is just "include.inc".
#
INCLUDES = include.inc
#
# OBJS is the list of all object files that need to be linked into the
# executable.  If a new source file is added, then the corresponding
# object file name must be added to the list below.
#
OBJS	= \
granadam.o \
input.o \
dresdenr.o \
dresdenc.o \
geome.o \
getgeo.o \
gettxt.o \
reada.o \
readint.o \
xyz2int.o \
upcase.o\
valstr.o
#
# For Sun arcstation 2 - version 4.1.3 of OS
#all:
#	@echo "Compiling for SPARC..."
#	@$(MAKE) FC='$(SPARCFC)' FFLAGS='$(SPARCFFLAGS)' CPP='$(SPARCCPP)' \
#	        EXE='$(SPARCEXE)' EXTRAS='$(SPARCEXTRAS)' linkit
#
# For Windows
#all:
#	@echo "Compiling for Windows..."
#	@$(MAKE) FC='$(FC)' FFLAGS='$(FFLAGS)' \
#		EXE='$(EXE)'  linkit
# For Linux
all:
	@echo "Compiling for Linux..."
	@$(MAKE) FC='$(FC)' FFLAGS='$(FFLAGS)' \
		EXE='$(EXE)'  linkit
#
# The default method of how to make object files from Fortran files.
#
.f.o:
#	$(CPP) $< TEMP.f
#	$(FC) -c $(FFLAGS) TEMP.f
	$(FC) -c $(FFLAGS) $<
#	/bin/mv TEMP.o $*.o
#	/bin/rm TEMP.[of]
#
#
linkit: $(OBJS)
	@echo " "
	@echo "Linking $(EXE)..."
	$(FC) $(LDFLAGS) -o $(EXE) $(OBJS)
	@echo " "
	@dir $(EXE)
#
# This removes all object files when "make clean" is entered by the user.
#
clean:
	rm -f *.o *.il
#
# Next comes all of the source file names and their dependencies (if any)
# and specific compile options (if any).  If no specific compile options
# are listed, then the default specified above is used.  If a new source file
# is added to GRANADA, then an entry should be included below
granadam.o:	granadam.f	$(INCLUDES)
input.o:	input.f		$(INCLUDES)
dresdenr.o:	dresdenr.f	$(INCLUDES)
dresdenc.o:	dresdenc.f	$(INCLUDES)
geome.o:	geome.f		$(INCLUDES)      
getgeo.o:	getgeo.f	$(INCLUDES)
gettxt.o:	gettxt.f
reada.o:	reada.f
readint.o:	readint.f	$(INCLUDES)   
xyz2int.o:	xyz2int.f
upcase.o:	upcase.f
valstr.o:	valstr.f
