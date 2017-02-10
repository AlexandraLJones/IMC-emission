


## MUST SET VARIABLES 
## I3RCDIR and MPIDIR to run 
## make

### Instalation directory
#MPIDIR = /opt/apps/intel13/mvapich2/1.9
 I3RC_Home = /mnt/a/u/sciteam/aljones4/I3RC
 CodeDir    = ${I3RC_Home}/src
 IntegDir   = ${I3RC_Home}/Integrators
 
### Netcdf-specific entries

 NetcdfHome = ${NETCDF_DIR}
#NetcdfHome = /home1/00478/tg457444/netcdf4.1.3
 Netcdf_IncludeDir = ${NetcdfHome}/include 
 NetcdfLibs = -L${NetcdfHome}/lib -lnetcdf 
#NetcdfLibs = -L$(NetcdfHome)/lib -lnetcdf -lnetcdff

### General compiler macros 

 ModuleFlag = -I
#ModuleFlag = 
# Modules     =  ${ModuleFlag}${Netcdf_IncludeDir}
Modules     = 
# Libs        =  ${NetcdfLibs} 
Libs        = 
 Compile.F95 =  ${F95} ${F95Flags} -c
 Compile.F77 =  $(F77) ${FFlags}   -c
 Link.F95    =  $(F95) ${F95Flags} 


#### Compiler-specific entries (may override marcos above) 

# Macros are available for ifort, g95, xlf, absoft, ftn
 compiler=ftn
 debug=no

ifeq (${compiler},ftn)
#   fortran compiler on Blue Waters
   F95         = ftn
   F77         = ftn
  ifeq (${debug},no)
  #optimization flags
     F95Flags = -dynamic${Modules} 
     FFlags =
  else
     F95Flags = -g -O fp0 -Rb -Rc -Rp -rm${Modules}
     FFlags = -g -O fp0 -rm
  endif
endif

ifeq ($(compiler),ifort)
  #
  # Intel fortran 
  #
  F95         = ifort
  F77         = ifort
  Compile.F95 += -free -Tf
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O3 -diag-disable vec -ipo $(Modules) 
    FFlags   = -fast 
  else
    # Debugging flags
    F95Flags = -g -O2  -C -traceback -check bounds -diag-disable vec $(Modules) 
    FFLAGS   = -g -O2  -C -traceback -check bounds -diag-disable vec 
  endif
endif

ifeq ($(compiler),g95)
  #
  # GNU fortran 
  #
  F95         = g95
  F77         = g77
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O2 -std=f95 $(Modules) 
    FFLAGS   = -O2
  else
    # Debugging flags
    F95Flags = -g -std=f95 -fmodule-private -fimplicit-none -Wall $(Modules)
    FFLAGS   = -g                           -fimplicit-none -Wall 
  endif
endif

ifeq ($(compiler),xlf)
  #
  # IBM xlf 8 for Mac 
  #
  F95         = f95
  F77         = f77
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O3 -qlanglvl=95std $(Modules) 
    FFLAGS   = -O3
  else
    # Debugging flags
    F95Flags = -g -O2 -C -qlanglvl=95std -qmaxmem=-1  $(Modules)
    FFLAGS   = -g                      
  endif
endif

ifeq ($(compiler),absoft)
  #
  # absoft fortran 
  #
  F95        = f95
  F77        = f77
  ModuleFlag = -p
  Libs      += -lU77
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O3  -stack_size 0x10000000 -cpu:host $(Modules)  
    FFlags   = -O3  -stack_size 0x10000000 -cpu:host 
  else 
    ##Debugging flags
    F95Flags = -g $(Modules) 
    FFLAGS   = -g 
  endif
endif


 UseMPI = yes
ifeq (${UseMPI},yes)
#  MPI_Dir = $(MPIDIR)
#  MPI_IncludeDir=$(MPI_Dir)/include
#  MPI_LibDir=$(MPI_Dir)/lib
  #Libs += $(MPI_LibDir) -lmpi 
  # MPI and non-MPI versions live in $(CodeDir)
  multipleProcCode = multipleProcesses_mpi.o
#  F95Flags += $(ModuleFlag)$(MPI_IncludeDir)
   F95Flags += ${Modules}
#  F95=ifort -lmpi
#  F77=ifort -lmpi
#   F95=$(MPIDIR)/bin/mpif90
#   F77=$(MPIDIR)/bin/mpif77

else
  multipleProcCode = multipleProcesses_nompi.o
endif


### General rules - should not require editing

# Rules for bullding object files
%.o: %.f
	${Compile.F77} $<

%.o: %.f95
	${Compile.F95} $<

# Rule to build executables	- some compilers are happier if everything 
#   is an object file 
%: %.o
	${Link.F95} -o $@ $^ ${Libs}

