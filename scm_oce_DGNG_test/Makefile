FC = gfortran
CC = gcc
FFLAGS = -c -fPIC
CFLAGS = -c -fPIC
f2pyFFLAGS = "-O4 -fdefault-real-8 -fdefault-double-8 -fPIC"

.PHONY: doc clean clean_all

GENERATED_FILES = scm_oce_DGNG_test.F90

all : scm_gls_DGNG_test.o scm_kpp.o scm_tke.o scm_oce_DGNG_test.so
	@echo
	@echo ===================================================
	@echo       Fortran SCM tools OK
	@echo ===================================================
	@echo

scm_oce_DGNG_test.so : scm_oce_DGNG_test.F90
	f2py -c scm_oce_DGNG_test.F90 scm_gls_DGNG_test.o scm_kpp.o scm_tke.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scm_oce_DGNG_test

scm_gls_DGNG_test.o : scm_gls_DGNG_test.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -fPIC -o scm_gls_DGNG_test.o -c scm_gls_DGNG_test.F90

scm_kpp.o : scm_kpp.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -fPIC -o scm_kpp.o -c scm_kpp.F90

scm_tke.o : scm_tke.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -fPIC -o scm_tke.o -c scm_tke.F90

clean:
	/bin/rm -f *.o *.mod *.so *.pyc *.msg* *~
clean_all:
	/bin/rm -f *.o *.mod *.so *.pyc *.msg* *~
