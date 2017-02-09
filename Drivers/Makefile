# $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
# $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Example-Drivers/Makefile $
include ../Makefile
Modules += $(ModuleFlag)$(NetcdfHome)/include $(ModuleFlag)$(CodeDir) $(ModuleFlag)$(IntegDir)

rtObjs =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
          $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
          $(CodeDir)/$(multipleProcCode) $(CodeDir)/scatteringPhaseFunctions.o \
          $(CodeDir)/inversePhaseFunctions.o $(CodeDir)/opticalProperties.o \
          $(CodeDir)/RandomNumbersForMC.o $(CodeDir)/monteCarloIllumination.o \
          $(CodeDir)/surfaceProperties.o $(IntegDir)/monteCarloRadiativeTransfer.o

Targets = monteCarloDriver planeParallel

all: $(Targets)

monteCarloDriver: monteCarloDriver.o $(rtObjs)
monteCarloDriver.o: monteCarloDriver.f95 $(rtObjs)

planeParallel: planeParallel.o $(rtObjs)
planeParallel.o: planeParallel.f95 $(rtObjs)

clean::
	/bin/rm *.lst *.mod *.o $(Targets)
