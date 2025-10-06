SRCS =	gemx_com.F90 equil.F90 gemx.F90 outd.F90 fcnt.F90

OBJS =	gemx_com.o equil.o gemx.o ionPush.o outd.o fcnt.o

PLIB = pputil.o

F90 = ftn

#OPT = -f free -O0 -s real64 -eD -hlist=a -llapack # -lblas -cpp # -Minfo=accel -acc -mp
#OPT = -f free -O3 -s real64  -hlist=a -llapack # -lblas -cpp # -Minfo=accel -acc -mp
OPT = -acc -Mbounds -Mfree -r8 -Minfo=accel -llapack -lblas -I/global/cfs/cdirs/mp118/software/petsc/install/include -L/global/cfs/cdirs/mp118/software/petsc/install/lib -lpetsc

LDFLAGS = 

#all : gemx

gemx: equil.o gemx.o ionPush.o outd.o fcnt.o pputil.o gemx_com.o
	$(F90)  -o gemx $(OPT) $(OBJS) $(PLIB) $(LIBS) 

pputil.o: pputil.F90
	$(F90) -c $(OPT) pputil.F90

gemx_com.o: gemx_com.F90 pputil.o
	$(F90) -c $(OPT) gemx_com.F90

equil.o: equil.F90 pputil.o gemx_com.o
	$(F90) -c $(OPT)  equil.F90

gemx.o: gemx.F90 pputil.o gemx_com.o equil.o
	$(F90) -c $(OPT) gemx.F90

ionPush.o: ionPush.F90 gemx_com.o equil.o
	$(F90) -c $(OPT) ionPush.F90

outd.o: outd.F90 gemx_com.o equil.o
	$(F90) -c $(OPT) outd.F90

fcnt.o: fcnt.F90
	$(F90) -c $(OPT) fcnt.F90

clean:
	rm -f *.o *.lst *.mod gem out/* dump/* xpp matrix/*  testne testphi psi_test run.out flag_debug debug.dat test* mask*

