F90 := ftn
CPP := cpp

ROOT := $(shell git rev-parse --show-toplevel)
BINDIR := $(ROOT)/bin
JOBDIR := $(ROOT)/jobs

SRCS := $(wildcard *.F90)
INCS_C := $(wildcard *.hpp)
SRCS_C := $(wildcard *.cpp) #C++ sources
OBJS := $(SRCS:.F90=.o)
OBJS_C := $(SRCS_C:.cpp=.o)
LIBS := -L${PETSC_PATH}/lib -lstdc++ -lpetsc
INCS := -I${PETSC_PATH}/include
DIFF := $(BINDIR)/codeChanges.txt
EXE := gemx

#Setup flags for various situations.
DEBUG ?= 0 #Use 'make DEBUG=1' to build with debug capabilities.
CPU ?= 0   #Use 'make CPU=1' to build for cpu.
OPT_COMMON := -r8 -Mfree $(INCS)
OPT_COMMON_C := $(INCS_C)
OPT_COMMON_REL := -O4
OPT_COMMON_DBG := -O0 -g -C -traceback -Mchkptr
OPT_CPU := -acc=host #Default to serial for cpu. -acc=multicore for parallel but should test against OpenMP on system.
OPT_GPU := -acc -Minfo=accel
OPT_CPU_DBG := -Mbounds #Note, bounds checking disabled on gpu by OpenACC.

#If not on perlmutter use local run script and local compiler. Default to CPU compilation for now.
ifneq ($(NERSC_HOST),perlmutter)
	F90 := mpif90
	CPP := mpicxx
	JOB := job_local
	CPU := 1
	ifeq ($(DEBUG),1)
		OPT_C := -O0 -g -Wall $(OPT_COMMON_C)
	else
		OPT_C := -O4 $(OPT_COMMON_C)
	endif
endif

ifeq ($(DEBUG),1)
	ifeq ($(CPU),1)
		JOB := job_dbg_cpu
		OPT :=  $(OPT_COMMON_DBG) $(OPT_COMMON) $(OPT_CPU) $(OPT_CPU_DBG)
	else
		JOB := job_dbg_gpu
		OPT :=  $(OPT_COMMON_DBG) $(OPT_COMMON) $(OPT_GPU)
	endif
else
	ifeq ($(CPU),1)
		JOB := job_cpu
		OPT := $(OPT_COMMON_REL) $(OPT_COMMON) $(OPT_CPU)
	else
		JOB := job_gpu
		OPT := $(OPT_COMMON_REL) $(OPT_COMMON) $(OPT_GPU)
	endif
endif

#Do all make tasks here so that nothing happens if up to date and not cleaned.
$(EXE): $(OBJS) $(OBJS_C)
	$(F90) -o $(EXE) $(OPT) $(LIBS) $(OBJS) $(OBJS_C)
	$(MAKE) install
	$(MAKE) check_diff

#Fortran compilation.
pputil.o: pputil.F90
	$(F90) -c $(OPT) pputil.F90

gemx_com.o: gemx_com.F90 pputil.o
	$(F90) -c $(OPT) gemx_com.F90

equil.o: equil.F90 pputil.o gemx_com.o
	$(F90) -c $(OPT) equil.F90

gemx.o: gemx.F90 pputil.o gemx_com.o equil.o
	$(F90) -c $(OPT) gemx.F90

ionPush.o: ionPush.F90 gemx_com.o equil.o
	$(F90) -c $(OPT) ionPush.F90

outd.o: outd.F90 gemx_com.o equil.o
	$(F90) -c $(OPT) outd.F90

fcnt.o: fcnt.F90
	$(F90) -c $(OPT) fcnt.F90

#C++ compilation.
fcnt_c.o: fcnt_c.cpp
	$(CPP) -c $(OPT_C) fcnt_c.cpp

gemx_c.o: gemx_c.cpp 
	$(CPP) -c $(OPT_C) gemx_c.cpp

#Check for git diffs in the src directory
check_diff:
	@echo "Checking for changes in $(ROOT)/src"
	@git diff $(ROOT)/src >> $(DIFF)
	@git rev-parse HEAD >> $(DIFF)
	@git branch --show-current >> $(DIFF)
	@echo "$(JOB)" >> "$(DIFF)"
	@echo "Diff saved in runs folder."

clean:
	$(RM) *.o *.lst *.mod gemx
	$(RM) -r $(BINDIR)

install:
	mkdir -p $(BINDIR) $(BINDIR)/out $(BINDIR)/dump $(BINDIR)/matrix
	cp -f gemx $(BINDIR)
	cp -f gemx.in $(BINDIR)
	cp -f $(JOBDIR)/$(JOB) $(BINDIR)
