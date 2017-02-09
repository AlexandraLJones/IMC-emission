! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

program ParticleFileToDomain
 ! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
 ! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Tools/PhysicalPropertiesToDomain.f95 $
 !  Modified by Alexandra Jones UIUC Fall 2011 to crreate a 3D temperature field from the vertical temperature distribution
!
 ! Reads a netCDF file containing mass contents and number concentrations (if provided)
 !   must also contain a temperature field (1D or 3D) and the vertical and horizontal 
 ! (or at least dx, dy) dimensions and writes the description to a
 !   "domain" object from the I3RC community model. Profiles of 
 !   molecular absoprtion can be added and Rayleigh scattering can 
 !   be included. 
 ! input mass should be in g/(m^3)  
! 
 ! Input parameters are specified with a namelist.
 !
 ! The scattering properties for each component are specified in
 ! files containing I3RC Monte Carlo phase function table objects
 ! (e.g. written by MakeMieTable).
 ! 
 ! An input file of molecular absorption extinction profile may be input.
 ! The Zlevels must be the same as the profile made by combining the
 ! levels in the particle file with the other levels specified in the
 ! namelist file.  Format (three lines): 
 !        nZ               [number of Z grid cells]
 !        Zlevels(1:nZ+1)  [increasing heights of cell boundaries in km]
 !        GasExt(1:nZ)     [molecular extinction in km^-1 for each cell]
 !
 ! Molecular Rayleigh scattering may be included (if RayleighWavelength > 0).
 ! The specified temperature profile is used to calculate the pressure
 ! profile with the hypsometric equation.  The Rayleigh extinction
 ! profile, which is proportional to air density, is calculated from
 ! the temperature and pressure profiles and the wavelength.  The average 
 ! extinction in each layer is calculated, assuming an exponential decay
 ! in air density.
 !

 !    Frank Evans    University of Colorado     December 2005
 !      Modified (for organization) by Robert Pincus, Climate Diagnostics Center


   ! Modules from the I3RC community Monte Carlo model framework
  use ErrorMessages
  use CharacterUtils
  use numericUtilities
  use scatteringPhaseFunctions
  use opticalProperties
  use UserInterface
  use netcdf
  
  implicit none

   ! Input parameters
  character(len=256)   :: NamelistFileName = ""
  character(len=256)   :: ParticleFileName = ""
  integer              :: numScatTables=0
  integer, parameter   :: maxNumComps=5
  character(len=256)   :: ScatTableFiles(maxNumComps) = ""
  character(len=256)   :: MolecAbsFileName = ""
  
  integer, parameter   :: maxOtherLevels = 20
  integer              :: numOtherLevels = 0
  real                 :: OtherHeights(maxOtherLevels) = 0, &
                          OtherTemps(maxOtherLevels) = 0
                          
  real                 :: RayleighWavelength = 0
  real                 :: DropNumConc(maxNumComps) = 0.
  character(len=256)   :: outputFileName = ""
  
  character(len=512)   :: tableDescription
  
  namelist /fileNames/ ParticleFileName, ScatTableFiles, MolecAbsFileName, &
                       outputFileName
  namelist /profile/  OtherHeights, OtherTemps
  namelist /physicalProperties/  DropNumConc, RayleighWavelength
  
   ! Local variables
  integer              :: nX, nY, nZp, nZt, Pfile_type, Tdim
  integer              :: i, j, k, n, iscat, ix, iy, iz, il
  integer              :: maxNretab, izLevelBase
  real(8)                 :: deltaX, deltaY, x, y, z
  real                  :: f
  real(8), allocatable    :: Zpar(:), TempPar(:), xEdges(:), yEdges(:), zEdges(:)
  integer, allocatable :: nComp(:,:,:), ptype(:,:,:,:)
  real, allocatable    :: MassCont(:,:,:,:), Reff(:,:,:,:), Pres(:,:,:)
  integer, allocatable :: Nretab(:)
  real, allocatable    :: ReffTable(:,:)
  real(8), allocatable    ::  ExtinctTable(:,:), ssaTable(:,:), temps(:,:,:), tempsOut(:,:,:)
  real(8), allocatable    :: Zlevels(:), Temp(:)
  real(8), allocatable    :: GasExt(:), RaylExt(:,:,:)
  real(8), allocatable    :: ssaProf(:,:,:), ssaProfGas(:)
  integer, allocatable :: phaseIndex(:,:,:), phaseIndexGas(:)
  real                 :: LegendreCoefs(2)
  real(8),    allocatable :: extinct(:,:,:,:), ssa(:,:,:,:)
  integer, allocatable :: phaseFuncIndex(:,:,:,:)

   ! I3RC Monte Carlo code derived type variables
  type(phaseFunction)                   :: phaseFuncObject
  type(phaseFunction)                   :: PhaseFuncs(2)
  type(phaseFunctionTable)              :: OnePhaseFuncTable
  type(phaseFunctionTable), allocatable :: phaseFuncTables(:)
  type(domain)                          :: thisDomain
  type(ErrorMessage)                    :: status


  ! -----------------------------------------
  ! Get the input variables from the namelist file
  !
  
  namelistFileName = getOneArgument()
  open (unit=1, file=trim(namelistFileName), status='OLD')
  read (1, nml = fileNames)
PRINT *, "read filenames"
  read (1, nml = profile)
PRINT *, "read profile"
  read (1, nml = physicalProperties)
PRINT *, "read properties"
  close (1)
  
  !
  ! Check input arguments as much as possible 
  !   Some checks can't be done until you know what kind of file you're reading
  !
  if(len_trim(ParticleFileName) == 0)   stop "Must specify particle file name." 
  if(len_trim(outputFileName) == 0) stop "Must specify output file name." 
 ! if(all(len_trim(ScatTableFiles) == 0))    stop "Must specify as many scattering tables as there are components"

  if(any(DropNumConc < 0)) stop "DropNumConc must be positive" 
  if(RayleighWavelength < 0) stop "RayleighWavelength must be non-negative." 
  
  
  numOtherLevels = count(otherTemps(:) > 0) 
  if(numOtherLevels == maxOtherLevels) &
    print *, "Read only the first ", maxOtherLevels, " other levels." 
  numScatTables = count(len_trim(ScatTableFiles) > 0)
  if(numScatTables == maxNumComps) &
    print *, "Read only the first ", maxNumComps, " scattering tables." 
   

  ! Read in the particle properties file
 
  call read_particle_file_size (ParticleFileName, nX, nY, nZt, Pfile_type, Tdim)
  allocate (zEdges(nZt+1), xEdges(nX+1), yEdges(nY+1), temps(nx, nY, nzt), Pres(nx, nY, nzt))
  allocate (nComp(nx,nY,nzt), ptype(numScatTables,nx,nY,nzt))
  allocate (MassCont(numScatTables,nx,nY,nzt), Reff(numScatTables,nx,nY,nzt))

  if (Tdim .eq. 1) then
   allocate(Temp(nZt))
   call read_particle_file1DT (ParticleFileName, Pfile_type, nX, nY, nZt, numScatTables, &
                           DropNumConc(1:numScatTables),  xEdges, yEdges, zEdges, Temp, Pres, &
                           nComp, ptype, MassCont, Reff)
   call create_temp_field (nZt, nX, nY, Temp, temps)
  else if (Tdim .eq. 3)then
   call read_particle_file3DT (ParticleFileName, Pfile_type, nX, nY, nZt, numScatTables, &
                           DropNumConc(1:numScatTables),  xEdges, yEdges, zEdges, temps, Pres, &
                           nComp, ptype, MassCont, Reff)
  else
   PRINT *, "Check dimensions of temperature field. must be either 1 or 3 but is: ", Tdim
  end if

!Reff(2,:,:,:) = 100.0  

if(any(temps(:,:,:) .lt. 0.0_8)) PRINT *, "temperatures negative"

  zEdges(nZt+1) = zEdges(nZt)+zEdges(2)
  yEdges(nY+1) = yEdges(nY)+yEdges(2)
  xEdges(nX+1) = xEdges(nX)+xEdges(2)

 ! convert from m to km
 PRINT *, "converting from meters to km"
 zEdges = zEdges/1000.0_8
 yEdges = yEdges/1000.0_8
 xEdges = xEdges/1000.0_8
! mass should stay in grams/(meter cubed) !
!MassCont = MassCont/(1000.0**3)

   ! Read the molecular absorption extinction file if there is one
  allocate (GasExt(nZt))
  GasExt(:) = 0
  if(len_trim(MolecAbsFileName) > 0) &
    call read_molec_abs_file (MolecAbsFileName, nZt, Zlevels, GasExt)

   ! Calculate the molecular Rayleigh scattering extinction profile
  allocate (RaylExt(nx,ny,nZt))
  RaylExt = 0
  if(RayleighWavelength > 0.) &
    call rayleigh_extinct (nx, ny, nzt, temps, Pres, RayleighWavelength, RaylExt)


  ! -----------------------------------------
  !  Read in the scattering tables
  !
 if(any(len_trim(ScatTableFiles) > 0))then
  allocate (phaseFuncTables(numScatTables), Nretab(numScatTables))
   ! Read the scattering tables and get the number of entries in each table
  do i = 1, numScatTables
PRINT *, ScatTableFiles(i)
    call read_PhaseFunctionTable(fileName = ScatTableFiles(i), &
                                 table = phaseFuncTables(i),   &
                                 status = status) 
    call printStatus(status)
    call getInfo_PhaseFunctionTable (phaseFuncTables(i), nEntries=Nretab(i), &
                                     tabledescription = tableDescription,    &
                                     status=status)
    call printStatus(status)
  enddo

   ! Get the effective radius, extinction, and single scattering albedo for
   !   all the entries in each scattering table
  maxNretab = maxval(Nretab(:))
  allocate (ReffTable(maxNretab,numScatTables))
  allocate (ExtinctTable(maxNretab,numScatTables), ssaTable(maxNretab,numScatTables))
  do i = 1, numScatTables
    call getInfo_PhaseFunctionTable (phaseFuncTables(i), &
                                     key                    = ReffTable(1:Nretab(i),i), &
                                     extinction             = ExtinctTable(1:Nretab(i),i), &
                                     singleScatteringAlbedo = ssaTable(1:Nretab(i),i), &
                                     status=status)
    call printStatus(status)
  enddo



  ! -----------------------------------------
  !  Use the effective radius for each particle type in each grid cell 
  !  to index into the scattering tables to calculate the extinction, 
  !  single scattering albedo, and phase function index fields.
  !  
  allocate (extinct(nX,nY,nZt,numScatTables))
  allocate (ssa(nX,nY,nZt,numScatTables))
  allocate (phaseFuncIndex(nX,nY,nZt,numScatTables))

  extinct(:,:,:,:) = 0.0
  ssa(:,:,:,:) = 0.0
  phaseFuncIndex(:,:,:,:) = 1
PRINT *, "max/min of ptype: ", maxval(ptype), minval(ptype)
PRINT *, "max/min of ncomp: ", maxval(ncomp), minval(ncomp)
PRINT *, "extinct dimension sizes, 1-4: ", size(extinct,1), size(extinct,2), size(extinct,3), size(extinct,4)
PRINT *, "numscatTables= ", numScatTables

  do iz = 1, nZt
   do iy = 1, nY
    do ix = 1, nX
      do k = 1, nComp(ix,iy,iz)
        ! Get the scattering table number (iscat) for this particle type
        iscat = ptype(k,ix,iy,iz)
        if(Reff(k,ix,iy,iz) <= maxval(ReffTable(1:NReTab(iscat),iscat)) .and. &
           Reff(k,ix,iy,iz) >  minval(ReffTable(1:NReTab(iscat),iscat))) then 
           
          ! Binary search to find effective radius entry in table
          il = findIndex(Reff(k,ix,iy,iz), ReffTable(1:NReTab(iscat),iscat))

           ! Interpolate optical properties linearly in Reff for this component
          f = (Reff(k,ix,iy,iz)-ReffTable(il,iscat)) / (ReffTable(il+1,iscat)-ReffTable(il,iscat))
          extinct(ix,iy,iz,iscat) = MassCont(k,ix,iy,iz) * &  
                                    ((1-f)*ExtinctTable(il,iscat) + f*ExtinctTable(il+1,iscat))
          ssa(ix,iy,iz,iscat) = (1-f)*ssaTable(il,iscat) + f*ssaTable(il+1,iscat)
	  if(ssa(ix,iy,iz,iscat) .gt. 1.0_8)then
	    ssa(ix,iy,iz,iscat) = 1.0_8
!PRINT *, "ix, iy, iz, iscat, ssa", ix, iy, iz, iscat, ssa(ix,iy,iz,iscat)
!STOP
	  end if
          ! Chose the closest phase function
          if (f < 0.5) then
            phaseFuncIndex(ix,iy,iz,iscat) = il
          else
            phaseFuncIndex(ix,iy,iz,iscat) = il+1
          endif
        
        else if(MassCont(k,ix,iy,iz) > 0.0) then 
          ! This effective radius isn't in the table
          print *, 'Warning: effective radius outside of table (ix,iy,iz,type,Reff):'
          print '(4(1x,i3),1x,f6.2)', ix, iy, iz, ptype(k,ix,iy,iz), Reff(k,ix,iy,iz)
          extinct(ix,iy,iz,iscat) = 0.
          ssa(ix,iy,iz,iscat) = 0.
        end if
      enddo
    enddo
   enddo
  enddo
  

  deallocate (MassCont, Reff, ptype, Ncomp)
  deallocate (Nretab, ReffTable, ExtinctTable, ssaTable)
 end if

  ! -----------------------------------------
  ! Package the optical properties in a domain object
  !
  !
  ! Create the domain
  !
  thisDomain = new_Domain (xEdges, yEdges, zEdges,temps, status)
  call printStatus(status)


   ! Add the optical properties for each particle component to the domain
  do i = 1, numScatTables
    print *, "Adding component ", i
    call addOpticalComponent (thisDomain, "Particle type " // trim(IntToChar(i)), &
                              extinct(:,:,:,i), &
                              ssa(:,:,:,i), phaseFuncIndex(:,:,:,i), &
                              phaseFuncTables(i),  &
                              status = status)
PRINT *, "returned from adding optical component"
    call printStatus(status)
    call finalize_PhaseFunctionTable (phaseFuncTables(i))
  end do
  if (ALLOCATED (extinct)) deallocate (extinct)
  if (ALLOCATED (ssa)) deallocate ( ssa) 
  if (ALLOCATED (phaseFuncIndex)) deallocate (phaseFuncIndex)
  if (ALLOCATED (phaseFuncTables)) deallocate ( phaseFuncTables)

   ! Add the Rayleigh scattering optical properties
  
  if (any(RaylExt(:,:,:) > 0.0)) then
    allocate (ssaProf(nx,ny,nZt), phaseIndex(nx,ny,nZt))
    print *, "Adding Rayleigh scattering"
    ssaProf = 1.0
    phaseIndex = 1
    LegendreCoefs(1:2) = (/ 0.0, 0.5 /) / (/ 2.*1. + 1., 2.*2. + 1. /)
    PhaseFuncs(1) = new_PhaseFunction (LegendreCoefs(1:2), status=status)
    call printStatus(status)
    OnePhaseFuncTable = new_PhaseFunctionTable (PhaseFuncs(1:1), key=(/ 0.0 /),&
                                                tableDescription = "Rayleigh scattering", &
                                                status=status)
PRINT *, "returned from adding Rayleigh"
    call printStatus(status)
    call finalize_phaseFunction (phaseFuncObject)
    call addOpticalComponent (thisDomain, 'Rayleigh scattering', &
                              RaylExt(:,:,:), ssaProf(:,:,:), phaseIndex(:,:,:), &
                              OnePhaseFuncTable, status = status)
    call printStatus(status)
    call finalize_PhaseFunctionTable (OnePhaseFuncTable)
  endif

   ! Add the molecular absorption optical properties
  if (any(GasExt(:) > 0.0)) then
    allocate (ssaProfGas(nZt), phaseIndexGas(nZt))
    print *, "Adding molecular absorption" 
    ssaProfGas(:) = 0.0
    phaseIndexGas(:) = 1
    LegendreCoefs(1:1) = 0.0
    PhaseFuncs(2) = new_PhaseFunction (LegendreCoefs(1:1), status=status)
    call printStatus(status)
    OnePhaseFuncTable = new_PhaseFunctionTable (PhaseFuncs(2:2), key=(/ 0.0 /),&
                                                tableDescription = "Molecular absorption", &
                                                status=status)
    call printStatus(status)
    call finalize_phaseFunction (phaseFuncObject)
    call addOpticalComponent (thisDomain, 'Molecular absorption', &
                              GasExt(:), ssaProfGas(:), phaseIndexGas(:), &
                              OnePhaseFuncTable, status = status)
PRINT *, "returned from adding gasExt"
    call printStatus(status)
    call finalize_PhaseFunctionTable (OnePhaseFuncTable)
  endif
  if(allocated(ssaProf)) deallocate (ssaProf)
  if(allocated(phaseIndex)) deallocate (phaseIndex)
  if(allocated(RaylExt)) deallocate(RaylExt)
  if(allocated(GasExt)) deallocate(GasExt)
  if(allocated(Zpar)) deallocate(Zpar)
  if(allocated(TempPar)) deallocate(TempPar)
  
  !
  ! Write the domain to a file
  !
  call write_Domain(thisDomain, outputFileName, status)
  call printStatus(status)
  call finalize_Domain(thisDomain)
end program ParticleFileToDomain
! --------------------------------------------------------------------------

subroutine read_particle_file_size (parfile, nx, ny, nzp, ftype, dims)

  use netcdf

!  implicit none
  character(len=*), intent(in) :: parfile
  integer, intent(out) :: nx, ny, nzp, ftype, dims

  integer, dimension(10)        :: ncStatus
  integer                       :: ncFileID, ncDimID, ncVarID
 
  ncStatus(:) = nf90_NoErr 
  if(nf90_open(trim(parfile), nf90_NoWrite, ncFileID) /= nf90_NoErr) then
      PRINT *, "read_particle_file_size: Can't open file " 
      STOP
  end if

!! Check to see if these are supposed to be number of edges or cells and adjust accordingly
!! make sure the names of the fields are correct
!! ALso, I may have to use varID commands instead of dimID commands, depending on how they're stored in the netCDF file
  ncStatus( 1) = nf90_inq_dimid(ncFileId, "NX", ncDimId)
  ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nx)
  ncStatus( 3) = nf90_inq_dimid(ncFileId, "NY", ncDimId)
  ncStatus( 4) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = ny)
  ncStatus( 5) = nf90_inq_dimid(ncFileId, "NZ", ncDimId)
  ncStatus( 6) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nzp)
  if(any(ncStatus(:) /= nf90_NoErr))then
    PRINT *, "read_particle_file_size: problem reading dimensions from netCDF file"
    STOP
  end if
  ncStatus( 7) = nf90_inq_varid(ncFileId, "NCw", ncVarId)  
  if( ncStatus( 7) /= nf90_NoErr) then
    PRINT *, "read_particle_file_size: problem reading number concentration. Assuming it is not provided and continuing with assumed number concentration provided in namelist"
    ftype = 1
  else
    PRINT *, "read_particle_file_size: detected presence of number conentration."
    ftype = 2
  end if
  ncStatus( 8) = nf90_inq_varid(ncFileId, "TMc", ncVarId)
  ncStatus( 9) = nf90_inquire_variable(ncFIleID, ncVarID,ndims=dims) 
  if( ncStatus( 8) /= nf90_NoErr .or. ncStatus( 9) /= nf90_NoErr) PRINT *, "read_particle_file_size: problem reading Temperature dimensions"

end subroutine read_particle_file_size

! --------------------------------------------------------------------------

subroutine read_particle_file1DT (parfile, filekind, nx, ny, nzp, nscattab, &
                               DropNumConc,  X, Y, Z, temppar, Pres, &
                               ncomp, ptype, masscont, reff)
  use netcdf

 ! currently set up to only deal with 1 scattering component

 ! The droplet number concentration (cm^-3) is used to derive the 
 ! effective radius from LWC for the 1 parameter LWC file (assumes a
 ! gamma distribution with alpha=7).
!  implicit none
  character(len=*), intent(in) :: parfile
  integer, intent(in) :: nx, ny, nzp, nscattab, filekind
  real,    intent(in) :: DropNumConc(nscattab)
  real(8),    intent(out) :: X(nx+1), Y(ny+1), Z(nzp+1), temppar(nzp)
  integer, intent(out) :: ncomp(nx,ny,nzp), ptype(nscattab,nx,ny,nzp)
  real,    intent(out) :: masscont(nscattab,nx,ny,nzp), reff(nscattab,nx,ny,nzp)
  real, dimension(nx,ny,nzp), intent(out)   :: Pres
  
  ! Local variables
  integer :: i, ic, ix, iy, iz, nc, pt(nscattab), ncFileID, ncVarID
  real    :: mass(nscattab), re(nscattab)
  real, parameter               :: Rd = 287.0
  integer, dimension(20)        :: ncStatus
  real, allocatable, dimension(:,:,:,:)  :: numConc

   ! Initialize the output arrays (in case there are missing grid points)
  ncomp(:,:,:) = 0
  ptype(:,:,:,:) = 0
  masscont(:,:,:,:) = 0.0
  reff(:,:,:,:) = 0.0
  
  ncStatus(:) = nf90_NoErr

  if(nf90_open(trim(parfile), nf90_NoWrite, ncFileID) /= nf90_NoErr) then
      PRINT *, "read_particle_file: Can't open file "
      STOP
  end if  

  ncStatus( 1) = nf90_inq_varid(ncFileId, "X", ncVarId)   ! units of meters
  ncStatus( 2) = nf90_get_var(ncFileId, ncVarId, X(1:nx))
  ncStatus( 3) = nf90_inq_varid(ncFileId, "Y", ncVarId)   ! units of meters
  ncStatus( 4) = nf90_get_var(ncFileId, ncVarId, Y(1:ny))
  ncStatus( 5) = nf90_inq_varid(ncFileId, "Z", ncVarId)   ! units of meters
  ncStatus( 6) = nf90_get_var(ncFileId, ncVarId, Z(1:nzp))
  ncStatus( 7) = nf90_inq_varid(ncFileId, "QCw", ncVarId)  ! units of g/kg
  ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, masscont(1,:,:,:))
  ncStatus( 9) = nf90_inq_varid(ncFileId, "TMc", ncVarId)  ! units of celcius
  ncStatus(10) = nf90_get_var(ncFileId, ncVarId, temppar)
  ncStatus(11) = nf90_inq_varid(ncFileId, "Prs", ncVarId) ! units of hPa
  ncStatus(12) = nf90_get_var(ncFileId, ncVarId, Pres)
  ncStatus(13) = nf90_inq_varid(ncFileId, "QCi", ncVarId)  ! units of g/kg
  ncStatus(14) = nf90_get_var(ncid=ncFileId, varid=ncVarId, values=masscont(2,:,:,:))

  if(any(ncStatus(:) /= nf90_NoErr)) then
     PRINT *, "read_particle_file1DT: problem reading  properties file."
     STOP
  end if

 ! convert pressure, temperature, mass content to proper units: Pa, K, g m^-3
 ! to convert mass you multiply by dry air density which can be determined from T, P, and Rd
  Pres = Pres *100.0
  temppar = temppar + 273.15
  if(any(temppar(:) .lt. 0.0_8)) PRINT *, "read_particle_file1DT: temperatures negative"
  forall (ix=1:nx, iy=1:ny, ic=1:nscattab)
    masscont(ic,:,iy,ix) = masscont(ic,:,iy,ix) * Pres(:,iy,ix) / (Rd * temppar)
  end forall

  select case(filekind)
    case(1) ! number concentration not provided in file so derive Reff from constant value
      forall (ic=1:nscattab)
        reff(ic,:,:,:) = 100* ( masscont(ic,:,:,:) *0.75*1.3889/(3.14159*DropNumConc(ic)) )**(1.0/3)
        ptype(ic,:,:,:) = ic
      end forall
      ncomp = nscattab

    case(2) ! number concentration provided in file, so derive Reff from that
      allocate(numConc(nscattab,nzp,ny,nx))
      ncStatus(1) = nf90_inq_varid(ncFileId, "NCw", ncVarId) ! units of hPa
      ncStatus(2) = nf90_get_var(ncFileId, ncVarId, numConc(1,:,:,:))
      ncStatus(3) = nf90_inq_varid(ncFileId, "NCi", ncVarId) ! units of hPa
      ncStatus(4) = nf90_get_var(ncFileId, ncVarId, numConc(2,:,:,:))

      if(any(ncStatus(:) /= nf90_NoErr)) then
        PRINT *, "read_particle_file1DT: problem reading number concentration."
        STOP
      end if
      reff(:,:,:,:) = 100* ( masscont(:,:,:,:) *0.75*1.3889/(3.14159*numConc(:,:,:,:)) )**(1.0/3)
      ncomp = 2 
      ptype(1,:,:,:) = 1
      ptype(2,:,:,:) = 2
 
  case default
    stop 'read_particle_file1DT: must be recognized type particle properties file.'
  end select 
  
end subroutine read_particle_file1DT
! -------------------------------------------------------------------------
subroutine read_particle_file3DT (parfile, filekind, nx, ny, nzp, nscattab, &
                               DropNumConc,  X, Y, Z, temppar, Pres, &
                               ncomp, ptype, masscont, reff)
 ! currently set up to only deal with 1 scattering component

 ! The droplet number concentration (cm^-3) is used to derive the
 ! effective radius from LWC for the 1 parameter LWC file (assumes a
 ! gamma distribution with alpha=7).
  use netcdf
!  implicit none
  character(len=*), intent(in) :: parfile
  integer, intent(in) :: nx, ny, nzp, nscattab, filekind
  real,    intent(in) :: DropNumConc(nscattab)
  real(8),    intent(out) :: X(nx+1), Y(ny+1), Z(nzp+1), temppar(nx,ny,nzp)
  integer, intent(out) :: ncomp(nx,ny,nzp), ptype(nscattab,nx,ny,nzp)
  real,    intent(out) :: masscont(nscattab,nx,ny,nzp), reff(nscattab,nx,ny,nzp)
  real, dimension(nx,ny,nzp), intent(out)   :: Pres

  ! Local variables
  integer :: i, ix, iy, iz, nc, pt(nscattab), ncFileID, ncVarID
  real    :: mass(nscattab), re(nscattab)
  real, parameter               :: Rd = 287.0
  integer, dimension(20)        :: ncStatus
  real, allocatable, dimension(:,:,:,:)  :: numConc

   ! Initialize the output arrays (in case there are missing grid points)
  ncomp(:,:,:) = 0
  ptype(:,:,:,:) = 0
  masscont(:,:,:,:) = 0.0
  reff(:,:,:,:) = 0.0

  ncStatus(:) = nf90_NoErr

  if(nf90_open(trim(parfile), nf90_NoWrite, ncFileID) /= nf90_NoErr) then
      PRINT *, "read_particle_file3DT: Can't open file "
      STOP
  end if

  ncStatus( 1) = nf90_inq_varid(ncFileId, "X", ncVarId)   ! units of meters
  ncStatus( 2) = nf90_get_var(ncFileId, ncVarId, X(1:nx))
  ncStatus( 3) = nf90_inq_varid(ncFileId, "Y", ncVarId)   ! units of meters
  ncStatus( 4) = nf90_get_var(ncFileId, ncVarId, Y(1:ny))
  ncStatus( 5) = nf90_inq_varid(ncFileId, "Z", ncVarId)   ! units of meters
  ncStatus( 6) = nf90_get_var(ncFileId, ncVarId, Z(1:nzp))
  ncStatus( 7) = nf90_inq_varid(ncFileId, "QCw", ncVarId)  ! units of g/kg
  ncStatus( 8) = nf90_get_var(ncid=ncFileId, varid=ncVarId, values=masscont(1,:,:,:))
  ncStatus( 9) = nf90_inq_varid(ncFileId, "TMc", ncVarId)  ! units of celcius
  ncStatus(10) = nf90_get_var(ncFileId, ncVarId, temppar(:,:,:))
  ncStatus(11) = nf90_inq_varid(ncFileId, "Prs", ncVarId) ! units of hPa
  ncStatus(12) = nf90_get_var(ncFileId, ncVarId, Pres(:,:,:))
  ncStatus(13) = nf90_inq_varid(ncFileId, "QCi", ncVarId)  ! units of g/kg
  ncStatus(14) = nf90_get_var(ncid=ncFileId, varid=ncVarId, values=masscont(2,:,:,:))

  if(any(ncStatus(:) /= nf90_NoErr)) then
     PRINT *, "read_particle_file3DT: problem reading  properties file. ", ncStatus(1:14)
!     PRINT *, trim(nf90_strerror(ncStatus(8)))
     STOP
  end if

 ! convert pressure, temperature, mass content to proper units: Pa, K, g m^-3
 ! to convert mass you multiply by dry air density which can be determined from T, P, and Rd
  Pres = Pres *100.0 ! convert from hPa to Pa
  temppar = temppar + 273.15 ! convert from Celcius to Kelvin
  if(any(temppar(:,:,:) .lt. 0.0_8))PRINT *, "read_particle_file3DT: temperatures negative"
  forall (i=1:nscattab)
    masscont(i,:,:,:) = masscont(i,:,:,:) * Pres / (Rd * temppar) ! convert from g/kg to g/(m^3)
  end forall

  select case(filekind)
    case(1) ! number concentration not provided in file so derive Reff from constant value
      forall (i=1:nscattab)
        reff(i,:,:,:) = 100* ( masscont(i,:,:,:) *0.75*1.3889/(3.14159*DropNumConc(i)) )**(1.0/3)
        ptype(i,:,:,:) = i
      end forall
      ncomp = 2

    case(2) ! number concentration provided in file, so derive Reff from that
      allocate(numConc(nscattab,nx,ny,nzp))
      ncStatus(1) = nf90_inq_varid(ncFileId, "NCi", ncVarId) ! units of hPa
      ncStatus(2) = nf90_get_var(ncFileId, ncVarId, numConc(2,:,:,:))
      ncStatus(3) = nf90_inq_varid(ncFileId, "NCw", ncVarId) ! units of hPa
      ncStatus(4) = nf90_get_var(ncFileId, ncVarId, numConc(1,:,:,:))
      if(any(ncStatus(:) /= nf90_NoErr)) then
        PRINT *, "read_particle_file3DT: problem reading number concentration."
        STOP
      end if
      reff(:,:,:,:) = 100* ( masscont(:,:,:,:) *0.75*1.3889/(3.14159*numConc(:,:,:,:)) )**(1.0/3)
      ncomp = 2
      ptype(1,:,:,:) = 1
      ptype(2,:,:,:) = 2

  case default
    stop 'read_particle_file3DT: must be recognized type particle properties file.'
  end select
end subroutine read_particle_file3DT
! --------------------------------------------------------------------------

subroutine organize_levels (nZp, Zpar, TempPar, &
                            numOtherLevels, OtherHeights, OtherTemps, &
                            nZt, Zlevels, Temp, izLevelBase)
 ! Adds the "other" layers to the layers from the particle file.
 ! The other heights must be outside the range of heights in the particle file.
 ! The base level index of the particle levels is returned in izLevelBase.
  implicit none
  integer, intent(in) :: nZp, numOtherLevels, nZt
  real(8),    intent(in) :: Zpar(nZp+1), TempPar(nZp+1)
  real,    intent(in) :: OtherHeights(numOtherLevels), OtherTemps(numOtherLevels)
  real(8),    intent(out) :: Zlevels(nZt+1), Temp(nZt+1)
  integer, intent(out) :: izLevelBase
  integer :: i, j, k

  
  if(any(Zpar(2:nZp+1) - zpar(1:nzp) <= 0)) &
    stop 'organize_levels: Zpar must increase'
  if(any(OtherHeights(:) >= Zpar(1) .and. OtherHeights(:) <= Zpar(nZp+1))) &
    stop 'organize_levels: OtherHeights must be outside particle file height range'

  do j = 1, numOtherLevels-1
    if (OtherHeights(j) >= OtherHeights(j+1)) then
      stop 'organize_levels: OtherHeights must increase'
    endif
  enddo
  
  k = 1  ;  j = 1
  do while (j<=numOtherLevels)
    if(OtherHeights(j) > Zpar(1)) exit
    Zlevels(k) = OtherHeights(j)
    Temp(k) = OtherTemps(j)
    k = k + 1  ;  j = j + 1
  enddo
  izLevelBase = k
  do i = 1, nZp+1
    Zlevels(k) = Zpar(i)
    Temp(k) = TempPar(i)
    k = k + 1
  enddo
  do while (j<=numOtherLevels)
    Zlevels(k) = OtherHeights(j)
    Temp(k) = OtherTemps(j)
    k = k + 1  ;  j = j + 1
  enddo
end subroutine organize_levels


! --------------------------------------------------------------------------

subroutine read_molec_abs_file (MolecAbsFileName, nZt, Zlevels, GasExt)
 ! Reads the three line ascii molecular absorption file which contains
 ! the extinction profile in the format:
 !        nZ               [number of Z grid cells]
 !        Zlevels(1:nZ+1)  [increasing heights of cell boundaries in km]
 !        GasExt(1:nZ)     [molecular extinction in km^-1 for each cell]
 ! If MolecAbsFileName is 'NONE' or '' then no file is read and the
 ! GasExt profile is set to zeros.
  implicit none
  character(len=*), intent(in) :: MolecAbsFileName
  integer, intent(in) :: nZt
  real(8),    intent(in) :: Zlevels(1:nZt+1)
  real(8),    intent(out) :: GasExt(1:nZt+1)
  integer :: nZ
  real(8), allocatable :: Zlevin(:)

  GasExt(:) = 0.0
  if (trim(MolecAbsFileName) /= 'NONE' .and. len_trim(MolecAbsFIlename) > 0) then
    open (unit=2, file=trim(MolecAbsFileName), status='old')
    read (2,*) nZ
    allocate (Zlevin(nZ+1))
    read (2,*) Zlevin(1:nZ+1)
    if (nZ /= nZt .or. any(abs(Zlevin(:) - Zlevels(:)) > spacing(Zlevels))) then
      print *, 'read_molec_abs_file: input Z levels do not match'
      print *, 'Zlevin=', Zlevin(:), 'Zlevels=', Zlevels(:)
      stop
    endif
    deallocate (Zlevin)
    read (2,*) GasExt(1:nZt)
  endif
end subroutine read_molec_abs_file


! --------------------------------------------------------------------------

subroutine rayleigh_extinct (nx, ny, nzt, Temps, Pres, wavelen, RaylExt)
 ! Computes the molecular Rayleigh extinction profile RaylExt [/km]
 ! from the temperature profile Temp [k] at Zlevels [km].  Assumes
 ! a linear lapse rate between levels to compute the pressure at
 ! each level.  The Rayleigh extinction is proportional to air
 ! density and depends on the wavelength [um].  The Rayleigh
 ! extinction is calculated at grid cell boundaries and then
 ! interpolated to get the correct average extinction assuming 
 ! an exponential behavior.
 ! The extinction profile is returned with zeros if wavelen<=0.
  implicit none
  integer, intent(in) :: nx, ny, nzt
  real ,    intent(in) :: wavelen
  real, dimension(nx, ny, nzt), intent(in)     :: Pres
  real(8), dimension(nx, ny, nzt), intent(in)     :: Temps
  real(8),    intent(out) :: RaylExt(nx,ny,nzt)
  
  real    :: raylcoef, Pres_hPa(nx,ny,nzt)
  integer  :: ix, iy, iz

! pressure needs to be in hPa 
  Pres_hPa = Pres/100.0
  raylcoef = 2.97E-4*wavelen**(-4.15+0.2*wavelen)
    RaylExt = raylcoef*Pres_hPa/Temps
  if(any(RaylExt(:,:,:)<0.0_8)) PRINT *, "WARNING: Rayleigh extinction negative"
  if(any(Pres_hpa(:,:,:)<0.0_8)) PRINT *, "WARNING: pressure negative"
  if(any(Temps(:,:,:)<0.0_8)) PRINT *, "WARNING: temperature negative"
end subroutine rayleigh_extinct

! --------------------------------------------------------------------------

subroutine create_temp_field(nZt, nx, ny, temp_in, temp_out)
 !fills in a 3D matrix of temperature field the same dimensions as the domain
 ! by averaging the vertical boundary temperature vector to grid cell centers.
 ! there is no horizontal variability

 implicit none
 integer, intent(in)   :: nZt, nx, ny
 real(8), intent(in)      :: temp_in(1:nZt)
 real(8), intent(out)     :: temp_out(1:nx,1:ny,1:nZt)

 integer               :: ix, iy

 forall (ix=1:nx, iy=1:ny)
    temp_out(ix,iy,:)=temp_in
 end forall

end subroutine create_temp_field
