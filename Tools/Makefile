# $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
# $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Tools/Makefile $
include ../Makefile
Modules += $(ModuleFlag)$(NetcdfHome)/include $(ModuleFlag)$(CodeDir) $(ModuleFlag)$(IntegDir)

mieObjs =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
           $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
           $(CodeDir)/scatteringPhaseFunctions.o $(CodeDir)/$(multipleProcCode)
           
optObjs = $(mieObjs) $(CodeDir)/inversePhaseFunctions.o $(CodeDir)/opticalProperties.o

all: MakeMieTable PhysicalPropertiesToDomain PhysicalPropertiesToDomainNew \
     netCDFPhysPropToDomainNew 
#OpticalPropertiesToDomain 

# Object files
#
MakeMieTable.o: MakeMieTable.f95  $(mieObjs)

PhysicalPropertiesToDomain.o: PhysicalPropertiesToDomain.f95 $(optObjs)

PhysicalPropertiesToDomainNew.o: PhysicalPropertiesToDomainNew.f95 $(optObjs)

#OpticalPropertiesToDomain.o: OpticalPropertiesToDomain.f95 $(optObjs)

netCDFPhysPropToDomainNew.o: netCDFPhysPropToDomainNew.f95 $(optObjs)

# Executables
#
MakeMieTable: MakeMieTable.o mieindsub.o RefractiveIndex-IceAndWater.o $(mieObjs)

PhysicalPropertiesToDomain: PhysicalPropertiesToDomain.o $(optObjs)

PhysicalPropertiesToDomainNew: PhysicalPropertiesToDomainNew.o $(optObjs)

#OpticalPropertiesToDomain: OpticalPropertiesToDomain.o $(optObjs)

netCDFPhysPropToDomainNew: netCDFPhysPropToDomainNew.o $(optObjs)

clean::
	-/bin/rm *.lst *.mod *.o MakeMieTable PhysicalPropertiesToDomain \
        PhysicalPropertiesToDomainNew netCDFPhysPropToDomainNew #OpticalPropertiesToDomain 
