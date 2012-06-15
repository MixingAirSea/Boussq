# Makefile for fluid solver
# uses the following libraries
# fftw3
# netcdf
# m
#
#
# Fortran Compiler
FC          = mpif90
#
# Fortran Flags
FFLAGS      = -ffree-form -fbounds-check -fdefault-real-8 -O5
#
# Local libraries
#MYLOC = /usr
MYLOC = /home/mhoecker/local
#
# Fortran include directories
FINC        = -I$(MYLOC)/include
#
# Fortran libraries
FLIB        = -L$(MYLOC)/lib -lfftw3 -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -lm
#
#

#-----------------------------------------------------------------

all : mpiffttest threedfluid

#------------------------------------------------------------------

mpiffttest : colrow.o mpiffttest.o
	$(FC) colrow.o mpiffttest.o $(FFLAGS) $(FINC) $(FLIB) -o mpiffttest

#------------------------------------------------------------------

threedfluid : threedstuff.o threedfluid.f threedNetCDF.o
	$(FC) $(FFLAGS) $(FINC) $(FLIB) -c threedfluid.f
	$(FC) threedstuff.o threedfluid.f threedNetCDF.o $(FFLAGS) $(FINC) $(FLIB) -o threedfluid

#------------------------------------------------------------------

threedstuff.o : threedstuff.f threedNetCDF.o
	$(FC) $(FFLAGS) $(FINC) $(FLIB) -c threedstuff.f

#------------------------------------------------------------------

threedNetCDF.o : threedNetCDF.f
	$(FC) $(FFLAGS) $(FINC) $(FLIB) -c threedNetCDF.f

#------------------------------------------------------------------

colrow.o : colrow.f
	$(FC) $(FFLAGS) $(FINC) $(FLIB) -c colrow.f

#------------------------------------------------------------------

doc: Doxyfile threedstuff.f threedfluid.f
	doxygen

#------------------------------------------------------------------
clean:
	rm -f core *.o *.mod
