# $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
# $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Integrators/Makefile $

include ../Makefile
Modules += $(ModuleFlag)$(CodeDir)

uses =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
        $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
        $(CodeDir)/RandomNumbersForMC.o \
	$(CodeDir)/scatteringPhaseFunctions.o $(CodeDir)/inversePhaseFunctions.o \
	$(CodeDir)/opticalProperties.o $(CodeDir)/monteCarloIllumination.o \
	$(CodeDir)/surfaceProperties.o

monteCarloRadiativeTransfer.o: monteCarloRadiativeTransfer.f95 $(uses)

clean::
	-/bin/rm *.lst *.mod *.o
