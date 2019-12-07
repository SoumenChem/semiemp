#
# $Id: GNUmakefile 25742 2014-06-08 07:38:13Z d3y133 $
#
   OBJ = \

   OBJ_OPTIMIZE = prop_chebyshev.o prop_os_superoperator.o prop_semiemp_euler.o prop_superoperator.o rt_semiemp_cs_cheb.o rt_semiemp_cs_prop.o rt_semiemp_os_cheb.o rt_semiemp_os_prop.o semiemp_cs_fock_build.o semiemp_os_fock_build.o\

   LIBRARY = libsemiemp.a

# If you add -DJEFF_OUTPUT to LIB_DEFINES, the dipole polarizability will be printed in the same format as the TCE.
# This is very useful (1) if you wrote scripts to parse that particular formatting and (2) if you want the anisotrophy.
#   THE LAPACK CALL ENABLED BY THIS IS BROKEN ON SOME ARCHITECTURES
#   LIB_DEFINES = -DJEFF_OUTPUT

   LIB_INCLUDES = -I../../../nwdft/include -I../../../ddscf

USES_BLAS = 


HEADERS = ../../../property/prop.fh ../../../property/atomdata.fh

include ../../../config/makefile.h
include ../../../config/makelib.h

#FDEBUG += -O0
