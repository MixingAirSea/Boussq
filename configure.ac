AC_INIT([hello],[1.0],[mhoecker@coas.oregonstate.edu])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([ foreign])
AC_LANG([Fortran])
AC_PROG_FC
AC_FC_FREEFORM
AC_FC_IMPLICIT_NONE
AC_FC_CHECK_BOUNDS
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
##
## Using local libraries
##
AC_ARG_WITH([local],
	[AS_HELP_STRING([--with-local=DIR],
	[Location of local libraries [default=/usr]])],
	,
	[withval=/usr]
)
case "$withval" in
	*,*)
		local_inc="`echo $withval |cut -f1 -d,`"
		local_lib="`echo $withval |cut -f2 -d, -s`"
	;;
	*)
		if test -n "$withval"; then
			local_inc="$withval/include"
			local_lib="$withval/lib"
		fi
	;;
esac
## -I/usr/include and -L/usr/lib are redundant
## can mess some compilers up.
if test "X$local_inc" = "X/usr/include"; then
	local_inc=""
fi
if test "X$local_lib" = "X/usr/lib"; then
	local_lib=""
fi
## Add local location to the compiler flags
if test -n "$local_inc"; then
	CPPFLAGS="$CPPFLAGS -I$local_inc"
	FCFLAGS="$FCFLAGS -I$local_inc"
	AM_CPPFLAGS="$AM_CPPFLAGS -I$local_inc"
fi
if test -n "$local_lib"; then
	LDFLAGS="$LDFLAGS -L$local_lib"
	AM_LDFLAGS="$AM_LDFLAGS -L$local_lib"
fi
##
## Test for libraries used by FFTW3
##
AC_LANG_PUSH([C])
AC_CHECK_LIB([m],[sin])
AC_LANG_POP([C])
##
## Using FFTW3
##
AC_ARG_WITH([fftw3],
	[AS_HELP_STRING([--with-fftw3=DIR],
	[Location of FFTW3 library [default=/usr]])],
	,
	[withval=/usr]
)
case "$withval" in
	*,*)
		fftw_inc="`echo $withval |cut -f1 -d,`"
		fftw_lib="`echo $withval |cut -f2 -d, -s`"
	;;
	*)
		if test -n "$withval"; then
			fftw_inc="$withval/include"
			fftw_lib="$withval/lib"
		fi
	;;
esac
## -I/usr/include and -L/usr/lib are redundant
## can mess some compilers up.
if test "X$fftw_inc" = "X/usr/include"; then
	fftw_inc=""
fi
if test "X$fftw_lib" = "X/usr/lib"; then
	fftw_lib=""
fi
## Add FFTW location to the compiler flags
if test -n "$fftw_inc"; then
	CPPFLAGS="$CPPFLAGS -I$fftw_inc"
	FCFLAGS="$FCFLAGS -I$fftw_inc"
	AM_CPPFLAGS="$AM_CPPFLAGS -I$fftw_inc"
fi
if test -n "$fftw_lib"; then
	LDFLAGS="$LDFLAGS -L$fftw_lib"
	AM_LDFLAGS="$AM_LDFLAGS -L$fftw_lib"
fi
##
## Test for FFTW3
##
AC_LANG_PUSH([C])
AC_CHECK_HEADERS([fftw3.h])
AC_LANG_POP([C])
AC_MSG_CHECKING([for fftw3.f])
AC_COMPILE_IFELSE(
[include "fftw3.f"],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])]
)
AC_MSG_CHECKING([for fftw3.f03])
AC_COMPILE_IFELSE(
[include "fftw3.f03"],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])]
)
AC_CHECK_LIB([fftw3],[dfftw_execute])

##
## Using zlib
##
AC_ARG_WITH([zlib],
	[AS_HELP_STRING([--with-zlib=DIR],
	[Location of zlib library [default=/usr]])],
	,
	[withval=/usr]
)
case "$withval" in
	*,*)
		zlib_inc="`echo $withval |cut -f1 -d,`"
		zlib_lib="`echo $withval |cut -f2 -d, -s`"
	;;
	*)
		if test -n "$withval"; then
			zlib_inc="$withval/include"
			zlib_lib="$withval/lib"
		fi
	;;
esac
## -I/usr/include and -L/usr/lib are redundant
## can mess some compilers up.
if test "X$zlib_inc" = "X/usr/include"; then
	zlib_inc=""
fi
if test "X$zlib_lib" = "X/usr/lib"; then
	zlib_lib=""
fi
## Add zlib location to the compiler flags
if test -n "$zlib_inc"; then
	CPPFLAGS="$CPPFLAGS -I$zlib_inc"
	AM_CPPFLAGS="$AM_CPPFLAGS -I$zlib_inc"
fi
if test -n "$zlib_lib"; then
	LDFLAGS="$LDFLAGS -L$zlib_lib"
	AM_LDFLAGS="$AM_LDFLAGS -L$zlib_lib"
fi
##
## Test for zlib
##
AC_LANG_PUSH([C])
AC_CHECK_HEADERS([zlib.h zconf.h])
AC_CHECK_LIB([z],[compress2])
AC_CHECK_FUNC([compress2])
AC_LANG_POP([C])
##
## Using HDF5
##
AC_ARG_WITH([hdf5],
	[AS_HELP_STRING([--with-hdf5=DIR],
	[Location of HDF5 library [default=/usr]])],
	,
	[withval=/usr]
)
case "$withval" in
	*,*)
		hdf_inc="`echo $withval |cut -f1 -d,`"
		hdf_lib="`echo $withval |cut -f2 -d, -s`"
	;;
	*)
		if test -n "$withval"; then
			hdf_inc="$withval/include"
			hdf_lib="$withval/lib"
		fi
	;;
esac
## -I/usr/include and -L/usr/lib are redundant
## can mess some compilers up.
if test "X$hdf_inc" = "X/usr/include"; then
	hdf_inc=""
fi
if test "X$hdf_lib" = "X/usr/lib"; then
	hdf_lib=""
fi
## Add hdf location to the compiler flags
if test -n "$hdf_inc"; then
	CPPFLAGS="$CPPFLAGS -I$hdf_inc"
	FCFLAGS="$FCFLAGS -I$hdf_inc"
	AM_CPPFLAGS="$AM_CPPFLAGS -I$hdf_inc"
fi
if test -n "$hdf_lib"; then
	LDFLAGS="$LDFLAGS -L$hdf_lib"
	AM_LDFLAGS="$AM_LDFLAGS -L$hdf_lib"
fi
##
## Test for HDF5
##
AC_LANG_PUSH([C])
AC_SEARCH_LIBS([H5Fflush], [hdf5dll hdf5], [], [AC_MSG_ERROR([Can not find or link to the hdf5 library])])
AC_SEARCH_LIBS([H5DSis_scale], [hdf5_hldll hdf5_hl], [], [AC_MSG_ERROR([Can not find or link to the hdf5 high-level library])])
AC_LANG_POP([C])
##
## Using NETCDF
##
AC_ARG_WITH([netcdf],
	[AS_HELP_STRING([--with-netcdf=DIR],
	[Location of NetCDF4 library [default=/usr]])],
	,
	[withval=/usr]
)
case "$withval" in
	*,*)
		netcdf_inc="`echo $withval |cut -f1 -d,`"
		netcdf_lib="`echo $withval |cut -f2 -d, -s`"
	;;
	*)
		if test -n "$withval"; then
			netcdf_inc="$withval/include"
			netcdf_lib="$withval/lib"
		fi
	;;
esac
## -I/usr/include and -L/usr/lib are redundant
## can mess some compilers up.
if test "X$netcdf_inc" = "X/usr/include"; then
	netcdf_inc=""
fi
if test "X$netcdf_lib" = "X/usr/lib"; then
	netcdf_lib=""
fi
## Add NetCDF location to the compiler flags
if test -n "$netcdf_inc"; then
	CPPFLAGS="$CPPFLAGS -I$netcdf_inc"
	FCFLAGS="$FCFLAGS -I$netcdf_inc"
	AM_CPPFLAGS="$AM_CPPFLAGS -I$netcdf_inc"
fi
if test -n "$netcdf_lib"; then
	LDFLAGS="$LDFLAGS -L$netcdf_lib"
	AM_LDFLAGS="$AM_LDFLAGS -L$netcdf_lib"
fi

##
## Check for NetCDF
##
AC_LANG_PUSH([C])
AC_CHECK_HEADER([netcdf.h])
AC_LANG_POP([C])
AC_MSG_CHECKING([for netcdf.inc])
AC_COMPILE_IFELSE(
[include "netcdf.inc"],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])]
)
AC_MSG_CHECKING([for netcdf.mod])
AC_COMPILE_IFELSE(
[include "netcdf.mod"],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])]
)
LIBS="-lnetcdf -lnetcdff $LIBS"
##
##
##
AC_ARG_WITH([mpi],
	[AS_HELP_STRING([--with-mpi=DIR],
	[Location of MPI implementation [default=/usr]])],
	,
	[withval=/usr]
)
case "$withval" in
	*,*)
		mpi_inc="`echo $withval |cut -f1 -d,`"
		mpi_lib="`echo $withval |cut -f2 -d, -s`"
	;;
	*)
		if test -n "$withval"; then
			mpi_inc="$withval/include"
			mpi_lib="$withval/lib"
		fi
	;;
esac
## -I/usr/include and -L/usr/lib are redundant
## can mess some compilers up.
if test "X$mpi_inc" = "X/usr/include"; then
	mpi_inc=""
fi
if test "X$mpi_lib" = "X/usr/lib"; then
	mpi_lib=""
fi
## Add mpi location to the compiler flags
if test -n "$mpi_inc"; then
	CPPFLAGS="$CPPFLAGS -I$mpi_inc"
	FCFLAGS="$FCFLAGS -I$mpi_inc"
	AM_CPPFLAGS="$AM_CPPFLAGS -I$mpi_inc"
fi
if test -n "$mpi_lib"; then
	LDFLAGS="$LDFLAGS -L$mpi_lib"
	AM_LDFLAGS="$AM_LDFLAGS -L$mpi_lib"
fi

##
## Check for MPI (this is lifted from the macro
## retrived from
## http://www.gnu.org/software/autoconf-archive/ax_mpi.html
## Feb 2013
##
AC_LANG_PUSH([C])
AC_CHECK_HEADER([mpi.h])
AC_LANG_POP([C])

AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
AC_CHECK_PROGS(MPIFC,[mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c],[$FC])
FC=$MPIFC
AC_SUBST(MPIFC)
AC_OUTPUT
