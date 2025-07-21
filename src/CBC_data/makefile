#FC=gfortran
#FC=ifort -O3 -g -free -pc64 -check all -watch all -warn all -vec-report3 -debug all -debug-parameters all -fp-stack-check
FC=ifort -O3
gcrz: module.o splines.o gcrz.o
	$(FC) module.o splines.o gcrz.o -o gcrz
module.o: module.f90
	$(FC) -c module.f90
splines.o: splines.f90
	$(FC) -c splines.f90
gcrz.o: gcrz.f90 
	$(FC) -c gcrz.f90
clean:
	rm *.mod module.o splines.o gcrz.o gcrz
