PROG = maximize_penning
FF = gfortran
F90 = gfortran
#FF  = /opt/intel/composer_xe_2011_sp1.6.233/bin/ifort #/opt/intel/fce/9.0/bin/ifort
#F90  = /opt/intel/composer_xe_2011_sp1.6.233/bin/ifort #/opt/intel/fce/9.0/bin/ifort
#FFLAGS = -O3 
#FF=/opt/intel/Compiler/11.1/069/bin/intel64/ifort
#F90=/opt/intel/Compiler/11.1/069/bin/intel64/ifort
#FFLAGS=-mcmodel=medium -shared-intel -O3 -Wl,-rpath,/opt/intel/Compiler/11.1/069/lib/intel64/

#FFLAGS=-mcmodel=medium-mcmodel=medium -O0 -g3
#gfortran->
FFLAGS= -ffree-line-length-none -O3 -fopenmp

#FFLAGS =  -g -assume minus0 -debug extended -debug-parameters all   -inline-debug-info -std90 -fpe0 -traceback -vec-report0 	 -mp1 -prec_div -O0 -C -warn all -ftrapuv -check noarg_temp_created
#FFLAGS = -g -w95  
#LIB = -llapack -lblas \
#	 /usr/software/lapack-3.1.1/ -ldfftpack -lsfftpack -lg2c
#LIB= /home/rosario/Software/Lapack/lapack-3.1.1/blas_LINUX.a  /home/omiste/software/ARPACK/libarpack_PROT.a /home/rosario/Software/Lapack/lapack-3.1.1/lapack_LINUX.a #/home/rosario/Software/nagfl15df/libnag_old.a

#LIB= -lblas -llapack -larpack #~/Software/ARPACK/libarpack_PROT.a #/home/rosario/Software/nagfl15df/libnag_old.a

#ifneq ($(wildcard ~/Programas/software/lapack-3.1.1/.),)
#  LIB=~/Programas/software/lapack-3.1.1/lapack_LINUX.a ~/Programas/software/lapack-3.1.1/blas_LINUX.a   ~/Programas/software/ARPACK/libarpack_SUN4.a
#else
#  LIB=-lblas -llapack -larpack
#endif


.SUFFIXES : .o .f90 .f
.f90.o: ; $(F90) $(FFLAGS) -c $<
.f.o: ; $(FF) $(FFLAGS) -c $<

OBJ = global.o nrtype.o nrutil.o nr.o fff.o f1dim_mod.o brent.o mnbrak.o optimization.o #sillib.o 

#FFT = fft4g.o

#OBJ1 = $(OBJ) $(FFT) propagation_asymmetric.o

OBJ1 = $(OBJ) main.o

all :: $(PROG)

maximize_penning : $(OBJ1) 
	$(F90) $(FFLAGS) -o $@ $(OBJ1) $(LIB)

#newmap:	 global.o parameter.o inout.o spline.o mapping.o pot.o numrec.o \
	mapint.o mapping.o newmap.o
#	$(F90) $(FFLAGS) -o $@ global.o parameter.o inout.o spline.o \
	mapping.o  pot.o  numrec.o mapint.o mapping.o newmap.o $(LIB)

#overlap :: global.o parameter.o inout.o spline.o mapping.o pot.o \
	numrec.o mapint.o mapsin.o grid.o  overlap.o
#	$(F90) $(FFLAGS) -o $@ global.o parameter.o inout.o spline.o \
	mapping.o pot.o numrec.o mapint.o mapsin.o grid.o  overlap.o

#testpulse :: global.o parameter.o inout.o spline.o $(FFT) \
	pulse.o testpulse.o
#	$(F90) $(FFLAGS) -o $@ global.o parameter.o inout.o spline.o $(FFT) \
	pulse.o testpulse.o $(LIB)

#interpol:: global.o parameter.o inout.o  spline.o mapping.o pot.o \
	numrec.o mapint.o mapsin.o grid.o angular.o\
#	psi.o mapinout.o inoutpsi.o work.o diag.o $(FFT) interpol.o
#	$(F90) $(FFLAGS) -o $@ global.o parameter.o inout.o  spline.o mapping.o \
	pot.o numrec.o mapint.o  mapsin.o grid.o angular.o psi.o mapinout.o \
#	inoutpsi.o work.o diag.o \
	interpol.o $(FFT) $(LIB)
#continterpol:: global.o parameter.o inout.o  spline.o mapping.o pot.o \
	numrec.o mapint.o mapsin.o grid.o angular.o \
#	psi.o mapinout.o inoutpsi.o work.o diag.o $(FFT) continterpol.o
#	$(F90) $(FFLAGS) -o $@ global.o parameter.o inout.o  spline.o mapping.o \
	pot.o numrec.o mapint.o  mapsin.o grid.o angular.o psi.o mapinout.o \
	inoutpsi.o work.o diag.o \
	continterpol.o $(FFT) $(LIB)

#printwavepacket:: global.o parameter.o inout.o spline.o mapping.o pot.o \
	numrec.o mapint.o mapsin.o grid.o \
	psi.o mapinout.o inoutpsi.o work.o printwavepacket.o
#	$(F90) $(FFLAGS) -o $@ global.o parameter.o inout.o spline.o mapping.o \
	pot.o numrec.o mapint.o mapsin.o grid.o psi.o mapinout.o inoutpsi.o work.o \
	printwavepacket.o $(FFT) $(LIB)

#rotconst :: global.o parameter.o inout.o spline.o mapping.o pot.o \
#	numrec.o mapint.o mapsin.o \
#	grid.o  psi.o pulse.o  $(FFT) work.o hamil.o rotconst.o
#	$(F90) -o $@ global.o parameter.o inout.o spline.o mapping.o pot.o \
#	numrec.o mapint.o mapsin.o grid.o  psi.o pulse.o $(FFT) \
#	work.o hamil.o rotconst.o  $(LIB)

#wig  :  global.o parameter.o inout.o spline.o mapping.o \
	numrec.o mapint.o mapsin.o  pot.o grid.o psi.o \
#	mapinout.o inoutpsi.o rho.o wigner.o wig.o  $(FFT)
#	$(F90) $(FFLAGS) -o $@ global.o parameter.o inout.o spline.o mapping.o \
	numrec.o mapint.o mapsin.o  pot.o grid.o \
	psi.o mapinout.o inoutpsi.o rho.o wigner.o wig.o  $(FFT) $(LIB)

#fourierwin :	global.o fourierwin.o
#	$(F90) $(FFLAGS) -o $@ global.o fourierwin.o $(LIB)

clean::	
	rm -f $(PROG) *.mod *.out *% *~ core global.o optimization.o nrtype.o nrutil.o nr.o f1dim_mod.o brent.o mnbrak.o fff.o 
