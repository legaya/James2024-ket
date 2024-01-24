FC = gfortran
CC = gcc
FFLAGS = -c -fPIC
CFLAGS = -c -fPIC
f2pyFFLAGS = "-O4 -fdefault-real-8 -fdefault-double-8"

.PHONY: doc clean clean_all

GENERATED_FILES = scm_oce_DG.F90

all : scm_gls_DG.o scm_kpp.o scm_tke.o scm_oce_DG.so
	@echo
	@echo ===================================================
	@echo       Fortran SCM tools OK
	@echo ===================================================
	@echo

scm_oce_DG.so : scm_oce_DG.F90
	f2py -c scm_oce_DG.F90 scm_gls_DG.o scm_kpp.o scm_tke.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scm_oce_DG

scm_gls_DG.o : scm_gls_DG.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -o scm_gls_DG.o -c scm_gls_DG.F90

scm_kpp.o : scm_kpp.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -o scm_kpp.o -c scm_kpp.F90

scm_tke.o : scm_tke.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -o scm_tke.o -c scm_tke.F90

clean:
	/bin/rm -f *.o *.mod *.so *.pyc *.msg* *~
clean_all:
	/bin/rm -f *.o *.mod *.so *.pyc *.msg* *~
