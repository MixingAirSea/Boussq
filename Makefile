# Makefile for fluid solver
# uses the following libraries
#
# fftw   >=3
# netcdf >=4
# libz   >=1
# hdf    >=5
# m
#
#
# Fortran Compiler
FC          = mpif90
#
# Fortran Flags
FFLAGS      = -ffree-form -fbounds-check -fdefault-real-8 -O5
FPROFFLAGS  = -ffree-form -fbounds-check -fdefault-real-8 -pg -O 
#
# Local libraries
#MYLOC = /usr
MYLOC       = /home/mhoecker/local
#
# Fortran include directories
FINC        = -I$(MYLOC)/include
#
# Fortran libraries
FLIB        = -L$(MYLOC)/lib -lfftw3 -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -lm
#
#

#-----------------------------------------------------------------

all : threedfluid

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

#-----------------------------------------------------------------

profile : mainprofile

#------------------------------------------------------------------

mainprofile : stuffprofile.o threedfluid.f NetCDFprofile.o
	$(FC) $(FPROFFLAGS) $(FINC) $(FLIB) -c threedfluid.f
	$(FC) stuffprofile.o  NetCDFprofile.o threedfluid.f $(FPROFFLAGS) $(FINC) $(FLIB) -o mainprofile

#------------------------------------------------------------------

stuffprofile.o : threedstuff.f NetCDFprofile.o
	$(FC) $(FPROFFLAGS) $(FINC) $(FLIB) -c threedstuff.f -o stuffprofile.o

#------------------------------------------------------------------

NetCDFprofile.o : threedNetCDF.f
	$(FC) $(FPROFFLAGS) $(FINC) $(FLIB) -c threedNetCDF.f -o NetCDFprofile.o

#------------------------------------------------------------------

NetCDFtest : pres_temp_4D_wr.f90
	$(FC) $(FPROFFLAGS) $(FINC) $(FLIB) -c pres_temp_4D_wr.f90
	$(FC) pres_temp_4D_wr.o $(FPROFFLAGS) $(FINC) $(FLIB)  -o pres_temp_4D_wr

#------------------------------------------------------------------

colrow.o : colrow.f
	$(FC) $(FFLAGS) $(FINC) $(FLIB) -c colrow.f

#------------------------------------------------------------------

doc: Doxyfile threedstuff.f threedfluid.f
	doxygen

#------------------------------------------------------------------
clean:
	rm -f core *.o *.mod

