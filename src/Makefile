F90 := ftn
CPP := CC

ROOT := $(shell git rev-parse --show-toplevel)
BINDIR := $(ROOT)/bin
JOBDIR := $(ROOT)/jobs

SRCS := $(wildcard *.F90)
#INCS_C := $(wildcard *.hpp *.h) 
INCS_C := -I$(MPIPATH)/include
SRCS_C := $(wildcard *.cpp) #C++ sources
OBJS := $(SRCS:.F90=.o)
OBJS_C := $(SRCS_C:.cpp=.o)
LIBS := -L${PETSC_PATH}/lib -lstdc++ -lpetsc -lmpi_cxx
INCS := -I${PETSC_PATH}/include
DIFF := $(BINDIR)/codeChanges.txt
EXE := gemx
SUPPRESS := --diag_suppress set_but_not_used --diag_suppress declared_but_not_referenced

#Setup flags for various situations.
DEBUG ?= 0 #Use 'make DEBUG=1' to build with debug capabilities.
CPU ?= 0   #Use 'make CPU=1' to build for cpu.
OPT_COMMON := -r8 -Mfree $(INCS)
OPT_COMMON_C := $(INCS_C) $(INCS)
OPT_COMMON_REL := -O4
OPT_COMMON_DBG := -O0 -g -C -traceback -Mchkptr
OPT_CPU := -acc=host #Default to serial for cpu. -acc=multicore for parallel but should test against OpenMP on system.
OPT_GPU := -acc -Minfo=accel #-gpu=cuda12.2 #(mine)
OPT_CPU_DBG := -Mbounds #Note, bounds checking disabled on gpu by OpenACC.

#If not on perlmutter use local run script and local compiler. Default to CPU compilation for now.
ifneq ($(NERSC_HOST),perlmutter)
	F90 := mpif90
	CPP := mpicxx
	JOB := job_local
	CPU := 1
	ifeq ($(DEBUG),1)
		OPT_C := -O0 -g -Wall -fno-omit-frame-pointer $(OPT_COMMON_C) $(INCS) $(SUPPRESS)
	else
		OPT_C := -O4 $(SUPPRESS) $(OPT_COMMON_C) $(OPT_GPU)
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
# 		OPT_C := -O4 $(SUPPRESS) $(OPT_CPU)
	else
		JOB := job_gpu
		OPT := $(OPT_COMMON_REL) $(OPT_COMMON) $(OPT_GPU)
		OPT_C:= -O4 $(SUPPRESS) $(OPT_GPU) 
	endif
endif

#Do all make tasks here so that nothing happens if up to date and not cleaned.
$(EXE): $(OBJS_C)
	$(CPP) $(OPT_GPU) -o $(EXE) $(LIBS) $(OBJS_C)
	$(MAKE) install
	$(MAKE) check_diff

#C++ compilation.

pputil_c.o: pputil_c.cpp
	$(CPP) -c $(OPT_C) pputil_c.cpp

gemx_com_c.o: gemx_com_c.cpp
	$(CPP) -c $(OPT_C) gemx_com_c.cpp

equil_c.o: equil_c.cpp
	$(CPP) -c $(OPT_C) equil_c.cpp 

gemx_c.o: gemx_c.cpp
	$(CPP) -c $(OPT_C) gemx_c.cpp 

ionPush_c.o: ionPush_c.cpp 
	$(CPP) -c $(OPT_C) ionPush_c.cpp

outd_c.o: outd_c.cpp
	$(CPP) -c $(OPT_C) outd_c.cpp

fcnt_c.o: fcnt_c.cpp
	$(CPP) -c $(OPT_C) fcnt_c.cpp

readDatFiles.o: readDatFiles.cpp
	$(CPP) -c $(OPT_C) readDatFiles.cpp

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
	
