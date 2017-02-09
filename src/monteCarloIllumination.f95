! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Code/monteCarloIllumination.f95 $
module monteCarloIllumination
  ! Provides an object representing a series of photons. The object specifies the 
  !   initial x, y position and direction. 
  ! In principle any arbitrary illumination condition can be specified. 
  ! Initial x, y positions are between 0 and 1, and must be scaled by the Monte Carlo
  !   integrator. 
  ! On input, azimuth is in degrees and zenith angle is specified as the cosine. 
  ! On output, azimuth is in radians (0, 2 pi) and solar mu is negative (down-going). 
  
  use ErrorMessages
  use RandomNumbers
  implicit none
  private

  !------------------------------------------------------------------------------------------
  ! Constants
  !------------------------------------------------------------------------------------------
  logical, parameter :: useFiniteSolarWidth = .false. 
  real,    parameter :: halfAngleDiamaterOfSun = 0.25 ! degrees
  
  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  type photonStream
    integer                     :: currentPhoton = 0
    real(8), dimension(:), pointer :: xPosition    => null()
    real(8), dimension(:), pointer :: yPosition    => null()
    real(8), dimension(:), pointer :: zPosition    => null()
    real, dimension(:), pointer :: solarMu      => null()
    real, dimension(:), pointer :: solarAzimuth => null()
  end type photonStream
  
  !------------------------------------------------------------------------------------------
  ! Overloading
  !------------------------------------------------------------------------------------------
  interface new_PhotonStream
    module procedure newPhotonStream_Directional, newPhotonStream_RandomAzimuth, &
                     newPhotonStream_Flux, newPhotonStream_Spotlight, newPhotonStream_LWemission
  end interface new_PhotonStream
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------
  public :: photonStream 
  public :: new_PhotonStream, finalize_PhotonStream, morePhotonsExist, getNextPhoton, emission_weighting
contains
  !------------------------------------------------------------------------------------------
  ! Code
  !------------------------------------------------------------------------------------------
  ! Initialization: Routines to create streams of incoming photons
  !------------------------------------------------------------------------------------------
  function newPhotonStream_Directional(solarMu, solarAzimuth, &
                                       numberOfPhotons, randomNumbers, status) result(photons)
    ! Create a set of incoming photons with specified initial zenith angle cosine and  
    !   azimuth. 
    real,                       intent(in   ) :: solarMu, solarAzimuth
    integer                                   :: numberOfPhotons
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(ErrorMessage),         intent(inout) :: status
    type(photonStream)                        :: photons
    
        ! Local variables
    integer :: i
    
    ! Checks: are input parameters specified correctly? 
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
    if(solarAzimuth < 0. .or. solarAzimuth > 360.) &
      call setStateToFailure(status, "setIllumination: solarAzimuth out of bounds")
    if(abs(solarMu) > 1. .or. abs(solarMu) <= tiny(solarMu)) &
      call setStateToFailure(status, "setIllumination: solarMu out of bounds")
    
    if(.not. stateIsFailure(status)) then
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons), photons%solarAzimuth(numberOfPhotons))
                        
      do i = 1, numberOfPhotons
        ! Random initial positions 
        photons%xPosition(   i) = getRandomReal(randomNumbers)
        photons%yPosition(   i) = getRandomReal(randomNumbers)
      end do
      photons%zPosition(:) = 1. - spacing(1.) 
      ! Specified inital directions
      photons%solarMu( :) = -abs(solarMu)
      photons%solarAzimuth(:) = solarAzimuth * acos(-1.) / 180. 
      photons%currentPhoton = 1
      
      call setStateToSuccess(status)
   end if   
  end function newPhotonStream_Directional
  ! ------------------------------------------------------
  function newPhotonStream_RandomAzimuth(solarMu, numberOfPhotons, randomNumbers, status) &
           result(photons)
    ! Create a set of incoming photons with specified initial zenith angle cosine but
    !  random initial azimuth. 
    real,                       intent(in   ) :: solarMu
    integer                                   :: numberOfPhotons
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(ErrorMessage),         intent(inout) :: status
    type(photonStream)                        :: photons
    
    ! Local variables
    integer :: i
    
    ! Checks: are input parameters specified correctly? 
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
    if(abs(solarMu) > 1. .or. abs(solarMu) <= tiny(solarMu)) &
      call setStateToFailure(status, "setIllumination: solarMu out of bounds")
    
    if(.not. stateIsFailure(status)) then
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons), photons%solarAzimuth(numberOfPhotons))
      do i = 1, numberOfPhotons
        ! Random initial positions 
        photons%xPosition(   i) = getRandomReal(randomNumbers)
        photons%yPosition(   i) = getRandomReal(randomNumbers)
        ! Random initial azimuth
        photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.) 
      end do
      photons%zPosition(:) = 1. - spacing(1.) 
      ! but specified inital mu
      photons%solarMu( :) = -abs(solarMu)
      photons%currentPhoton = 1
      
      call setStateToSuccess(status)
    end if   
 end function newPhotonStream_RandomAzimuth
 ! ------------------------------------------------------
 function newPhotonStream_Flux(numberOfPhotons, randomNumbers, status) result(photons)
    ! Create a set of incoming photons with random initial azimuth and initial
    !  mus constructed so the solar flux on the horizontal is equally weighted
    !  in mu (daytime average is 1/2 solar constant; this is "global" average weighting)
    integer                                   :: numberOfPhotons
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(ErrorMessage),         intent(inout) :: status
    type(photonStream)                        :: photons
    
    ! Local variables
    integer :: i
    
    ! Checks
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
     
    if(.not. stateIsFailure(status)) then
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons),  photons%solarAzimuth(numberOfPhotons))
               
      do i = 1, numberOfPhotons
        ! Random initial positions
        photons%xPosition(   i) = getRandomReal(randomNumbers)
        photons%yPosition(   i) = getRandomReal(randomNumbers)
        ! Random initial directions
        photons%solarMu(     i) = -sqrt(getRandomReal(randomNumbers)) 
        photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.)
      end do
      photons%zPosition(:) = 1. - spacing(1.) 
      
      photons%currentPhoton = 1
      call setStateToSuccess(status)  
    end if     
  end function newPhotonStream_Flux
  !------------------------------------------------------------------------------------------
  function newPhotonStream_Spotlight(solarMu, solarAzimuth, solarX, solarY, &
                                     numberOfPhotons, randomNumbers, status) result(photons)
    ! Create a set of incoming photons with specified initial zenith angle cosine and  
    !   azimuth. 
    real,                       intent(in   ) :: solarMu, solarAzimuth, solarX, solarY
    integer                                   :: numberOfPhotons
    type(randomNumberSequence), optional, &
                                intent(inout) :: randomNumbers
    type(ErrorMessage),         intent(inout) :: status
    type(photonStream)                        :: photons
        
    ! Checks: are input parameters specified correctly? 
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
    if(solarAzimuth < 0. .or. solarAzimuth > 360.) &
      call setStateToFailure(status, "setIllumination: solarAzimuth out of bounds")
    if(abs(solarMu) > 1. .or. abs(solarMu) <= tiny(solarMu)) &
      call setStateToFailure(status, "setIllumination: solarMu out of bounds")
    if(abs(solarX) > 1. .or. abs(solarX) <= 0. .or. &
       abs(solarY) > 1. .or. abs(solarY) <= 0. )    &
      call setStateToFailure(status, "setIllumination: x and y positions must be between 0 and 1")
    
    if(.not. stateIsFailure(status)) then
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons), photons%solarAzimuth(numberOfPhotons))
                        
      ! Specified inital directions and position
      photons%solarMu( :) = -abs(solarMu)
      photons%solarAzimuth(:) = solarAzimuth * acos(-1.) / 180. 
      photons%xPosition(:) = solarX
      photons%yPosition(:) = solarY
      photons%zPosition(:) = 1. - spacing(1.) 
      photons%currentPhoton = 1
      
      call setStateToSuccess(status)
   end if   
  end function newPhotonStream_Spotlight

!-------------------------------------------------------------------------------------------
  function newPhotonStream_LWemission(numberOfPhotons, atms_photons, voxel_weights, col_weights, level_weights, nx, ny, nz, randomNumbers, status, option1) result(photons)
    ! Create a set of emitted photons with random initial azimuth, random 
    ! mus, and random x,y,z location within the domain. This is the LW source from the atmosphere
    ! and surface. The x,y,z locations are weighted based on the power emitted 
    ! from each  voxel, pixel, in both the atmosphere and surface.
    ! Written by Alexandra Jones at the University of Illinois, Urbana-Champaign. Fall 2011
    ! Updated Fall 2012 to remove predetermination of number of photons emitted per column
    implicit none
    
    integer, intent(in)                             :: numberOfPhotons, nx, ny, nz
    type(randomNumberSequence), intent(inout)       :: randomNumbers
    type(ErrorMessage),         intent(inout)       :: status
    type(photonStream)                              :: photons
    integer,                       intent(in)       :: atms_photons
    real(8), dimension(nx,ny,nz), intent(in)         :: voxel_weights
    real(8), dimension(ny,nz), intent(in)         :: col_weights
    real(8), dimension(nz), intent(in)         :: level_weights
    integer, dimension(nx,ny,nz), optional, intent(out)   :: option1

    ! Local variables
    integer :: i, numberOfAtmsPhotons, startPhoton,  ii, ij, ik
    real    :: RN!, test
!    real, dimension(0:nx*ny*nz)                      :: temp1, temp2
    
    if(present(option1))option1=0

    ! Checks
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
       
    if(.not. stateIsFailure(status)) then
!PRINT *, "number of photons > 0"
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons),  photons%solarAzimuth(numberOfPhotons))

              
     ! divide the photons into sfc vs. atms sources 
      numberOfAtmsPhotons=MIN(numberOfPhotons,atms_photons)
      if ( numberOfAtmsPhotons .gt. 0)then
!PRINT *, "numberOfAtmsPhotons .gt. 0"
          startPhoton = 1
     
          do i = startPhoton,  numberOfAtmsPhotons ! Loop over photons from atms source 
!           do while (i .le. numberOfPhotoins) 
            ! Random initial positions
            RN = getRandomReal(randomNumbers)
!if(i .eq. 1)PRINT *, "i=", i, " RN=", RNi          
            DO ik=1,nz
              if(RN .le. level_weights(ik))then
                DO ij=1,ny
                  if(RN .le. col_weights(ij,ik))then
                    DO ii=1,nx
                      if(RN .le. voxel_weights(ii,ij,ik))then
                        exit
                      endif     
                    ENDDO
                    exit
                  endif
                ENDDO
                exit
              endif
            ENDDO

!PRINT *, RN, ii, ij, ik
!write(16,"(I5 ,2X, E30.20)") i, voxel_weights(ii,ij,ik)-RN

 !if(present(option1))then
 !   option1(ii,ij,ik)=option1(ii,ij,ik)+1
 !end if 

            photons%zPosition(i) = ((ik-1)*(1.0_8)/nz) + dble(getRandomReal(randomNumbers)/nz) ! The first term represents the fractional position of the layer bottom in the column, such that ik=1 corresponds to a position of 0. The second term respresents the position within the layer.
            if(ik .eq. 1 .and. photons%zPosition(i) .eq. 0.0_8) photons%zPosition(i)=0.0_8+spacing(1.0_8)
            if(ik .eq. nz .and. photons%zPosition(i) .gt. 1.0_8-2.0_8*spacing(1.0_8)) photons%zPosition(i)=photons%zPosition(i) - (2.0_8*spacing(1.0_8))
            photons%xPosition(i) = ((ii -1)*1.0_8/nx) + dble(getRandomReal(randomNumbers)*(1.0/nx)) 
            photons%yPosition(i) = ((ij -1)*1.0_8/ny) + dble(getRandomReal(randomNumbers)*(1.0/ny)) 
!if(i .eq. 1) PRINT *, 'ind= ', ind, 'i= ', ik, 'j= ', ij, 'k= ', ik, 'xPos= ', photons%xPosition(i)
            ! Random initial directions
!            test=getRandomReal(randomNumbers)
!PRINT *, 'i=', i, ' atms_rand=', test 
!            photons%solarMu(i) = 1-(2.*test)
          DO
            photons%solarMu(i) = 1-(2.*getRandomReal(randomNumbers))    ! This formula is from section 10.5 of '3D radiative transfer in cloudy atmospheres'. These angles will stay the same...truly random numbers, since LW emission is isotropic. But the Mu no longer has to be negative. The name "solarMu" should really just be "sourceMu" 
            if(abs(photons%solarMu(i)) > 2 * tiny (photons%solarMu(i)))exit  ! This ensures that there is some vertical component to the photon trajectory, so the photon doesn't get permanently stuck in the layer it's initialized in
          END DO
            photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.)
            if (numberOfPhotons < numberOfAtmsPhotons)exit
!           end do
          end do  !loop over i

     end if

     if (numberOfPhotons > numberOfAtmsPhotons)then
      do i = numberOfAtmsPhotons+1, numberOfPhotons ! Loop over sfc source
        ! Random initial positions
        photons%xPosition(   i) = getRandomReal(randomNumbers) ! Assuming the surface temperature and emissivities are the same everywhere the x,y values don't need to be weighted
        photons%yPosition(   i) = getRandomReal(randomNumbers)
        ! Random initial directions
!        test = getRandomReal(randomNumbers)
!PRINT *, 'i=', i, ' sfc_rand=', test
!         photons%solarMu(     i) = sqrt(test)
       DO
        photons%solarMu(     i) = sqrt(getRandomReal(randomNumbers))    ! This formula is from section 10.5 of '3D radiative trasnfer in cloudy atmospheres'. These angles will stay the same...truly random numbers, since LW emission is isotropic. But the Mu has to be positive for a sfc source. The name "solarMu" should really just be "sourceMu"
        if(abs(photons%solarMu(i)) > 2 * tiny (photons%solarMu(i)))exit  ! This ensures that there is some vertical component to the photon trajectory, so the photon doesn't get permanently stuck in the lowest layer. this is especially bad when there is no atmospheric extinction because then the photon will never leave the domain or hit the surface or become extinct. Thus we enter an infinite loop
       END DO
        photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.)
      end do  !loop over i
        photons%zPosition(numberOfAtmsPhotons+1:numberOfPhotons) = 0.0_8 ! must be from the sfc. Make sure this value makes sense for how the position is interpretted
     end if
!PRINT *,  ' photons%solarMu=', photons%solarMu(1)
      photons%currentPhoton = 1
      call setStateToSuccess(status) 
    end if    
  end function newPhotonStream_LWemission

   
  !------------------------------------------------------------------------------------------
  ! Are there more photons? Get the next photon in the sequence
  !------------------------------------------------------------------------------------------
  function morePhotonsExist(photons)
    type(photonStream), intent(inout) :: photons
    logical                           :: morePhotonsExist
    
    morePhotonsExist = photons%currentPhoton > 0 .and. &
                       photons%currentPhoton <= size(photons%xPosition)
  end function morePhotonsExist
  !------------------------------------------------------------------------------------------
  subroutine getNextPhoton(photons, xPosition, yPosition, zPosition, solarMu, solarAzimuth, status)
    type(photonStream), intent(inout) :: photons
    real(8),                 intent(  out) :: xPosition, yPosition, zPosition
    real,                 intent(  out) ::  solarMu, solarAzimuth
    type(ErrorMessage),   intent(inout) :: status
    
    ! Checks 
    ! Are there more photons?
    if(photons%currentPhoton < 1) &
      call setStateToFailure(status, "getNextPhoton: photons have not been initialized.")  
    if(.not. stateIsFailure(status)) then
      if(photons%currentPhoton > size(photons%xPosition)) &
        call setStateToFailure(status, "getNextPhoton: Ran out of photons")
    end if
      
    if(.not. stateIsFailure(status)) then
      xPosition    = photons%xPosition(photons%currentPhoton) 
      yPosition    = photons%yPosition(photons%currentPhoton)
      zPosition    = photons%zPosition(photons%currentPhoton)
      solarMu      = photons%solarMu(photons%currentPhoton)
!PRINT *, 'solarMu=', solarMu
      solarAzimuth = photons%solarAzimuth(photons%currentPhoton)
      photons%currentPhoton = photons%currentPhoton + 1
    end if
     
  end subroutine getNextPhoton
  !------------------------------------------------------------------------------------------
  ! Finalization
  !------------------------------------------------------------------------------------------
  subroutine finalize_PhotonStream(photons)
    type(photonStream), intent(inout) :: photons
    ! Free memory and nullify pointers. This leaves the variable in 
    !   a pristine state
    if(associated(photons%xPosition))    deallocate(photons%xPosition)
    if(associated(photons%yPosition))    deallocate(photons%yPosition)
    if(associated(photons%zPosition))    deallocate(photons%zPosition)
    if(associated(photons%solarMu))      deallocate(photons%solarMu)
    if(associated(photons%solarAzimuth)) deallocate(photons%solarAzimuth)
    
    photons%currentPhoton = 0
  end subroutine finalize_PhotonStream
!---------------------------------------------------------------------------------------------
! Compute number of photons emmitted from domain and surface when LW illumination is used
!---------------------------------------------------------------------------------------------
  subroutine emission_weighting(nx, ny, nz, numComponents, xPosition, yPosition, zPosition, &
				lambda_u, totalPhotons, atmsPhotons, voxel_weights, col_weights, &
				level_weights, atmsTemp, ssas, cumExt, ext, sfcTemp, &
				emiss, totalFlux)
!Computes Planck Radiance for each surface and atmosphere pixel to determine the wieghting for the distribution of photons.
!Written by ALexandra Jones, University of Illinois, Urbana-Champaign, Fall 2011
! Updated Fall 2012 to remove predetermination of number of photons emitted per column
     implicit none

     integer, intent(in)                               :: nx, ny, nz, numComponents, totalPhotons
     real(8), dimension(1:nx, 1:ny, 1:nz), intent(in)     :: atmsTemp, cumExt
     real(8), dimension(1:nx,1:ny,1:nz,1:numComponents), intent(in) :: ssas, ext
     real(8), dimension(1:nz+1), intent(in)               :: zPosition
     real(8), dimension(1:ny+1), intent(in)               :: yPosition
     real(8), dimension(1:nx+1), intent(in)               :: xPosition
     real(8), intent(in)                                  :: sfcTemp, emiss, lambda_u
     integer,  intent(out)                             :: atmsPhotons
     real(8), dimension(nx,ny,nz), intent(out)          :: voxel_weights
     real(8), dimension(ny,nz), intent(out)          :: col_weights
     real(8), dimension(nz), intent(out)          :: level_weights
     real(8),                                intent(out)  :: totalFlux
     !Local variables
     integer                                           :: ix, iy, iz!, last
     real(8)                                              ::  sfcPlanckRad, sfcPower,  atmsPower, totalPower, totalAbsCoef, b, lambda
     real(8)                                          :: previous, corr_contrib,corr,temp_sum, prev_exact
     real(8), dimension(1:nz)                          :: dz
     real(8), dimension(1:ny)                             :: dy
     real(8), dimension(1:nx)                             :: dx
     real(8)                                               :: atmsPlanckRad
!     real, dimension(1:nx, 1:ny)                       :: atmsColumn_power

     real(8), parameter                                   :: h=6.62606957e-34 !planck's constant [Js]
     real(8), parameter                                   :: c=2.99792458e+8 !speed of light [ms^-1]
     real(8), parameter                                   :: k=1.3806488e-23 !boltzman constant [J/K molecule]
     real(8), parameter                                   :: a=2.0_8*h*c**2.0_8  
     real(8), parameter                                   :: Pi=4*DATAN(1.0_8)

     lambda=lambda_u/(10.0_8**6.0_8) ! convert lambda from micrometers to meters
     b=h*c/(k*lambda)

!PRINT *, h, c, k, lambda, Pi, a, b
!calculate arrays of depths from the position arrays in km
     dz(1:nz)=zPosition(2:nz+1)-zPosition(1:nz)
     dy(1:ny)=yPosition(2:ny+1)-yPosition(1:ny)
     dx(1:nx)=xPosition(2:nx+1)-xPosition(1:nx)


!dz(1:nz)= 0.04					! be sure to remove this line after debugging FOR DIAGNOSTIC PURPOSES ONLY!

!     last=nx*ny*nz ! the value of the index of the last element of the voxel_weights array

!first compute atms planck radiances then combine algorithms from mcarWld_fMC_srcDist and mcarWld_fMC_srcProf to determine the  weights of each voxel taking into consideration the ones that would be emitted from the surface instead.
     if (emiss .eq. 0.0_8 .or. sfcTemp .eq. 0.0_8)then
        sfcPower=0.0_8
     else
        sfcPlanckRad=(a/((lambda**5.0_8)*(exp(b/sfcTemp)-1.0_8)))/(10.0_8**6.0_8)
!PRINT *, "sfcPlanckRad=", sfcPlanckRad, "sfcTemp", sfcTemp
        sfcPower= Pi*emiss*sfcPlanckRad*(xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8)     ! [W] factor of 1000^2 needed to convert area from km to m
     end if

     atmsPower=0.0_8
     voxel_weights=0.0_8  
     level_weights=0.0_8
     col_weights=0.0_8
     previous=0.0_8
     corr_contrib=0.0_8
     temp_sum=0.0_8
     corr=0.0_8
     prev_exact=0.0_8
    if(COUNT(atmsTemp .le. 0.0_8) .eq. 0)then
     do iz = 1, nz
       do iy = 1, ny
         do ix = 1, nx
           atmsPlanckRad= (a/((lambda**5.0_8)*(exp(b/atmsTemp(ix,iy,iz))-1.0_8)))/(10.0_8**6.0_8) ! the 10^-6 factor converts it from Wsr^-1m^-3 to Wm^-2sr^-1micron^-1
           totalAbsCoef=cumExt(ix,iy,iz)-sum(ssas(ix,iy,iz,:) * ext(ix,iy,iz,:))
!PRINT *, "atmsPlanckRad=", atmsPlanckRad, "totalAbsCoef=", totalAbsCoef, "lambda=", lambda, "atmsTemp=", atmsTemp(ix,iy,iz), "Pi=", Pi
!if (ix .eq. 1 .and. iy .eq. 1 .and. iz .eq. 1) PRINT *, cumExt(ix,iy,iz), ssas(ix,iy,iz,:), ext(ix,iy,iz,:), sum(ssas(ix,iy,iz,:) * ext(ix,iy,iz,:)), totalAbsCoef
           corr_contrib = (4.0_8*Pi* atmsPlanckRad * totalAbsCoef*dz(iz))-corr     ! [Wm^-2] 
           temp_sum = previous + corr_contrib
           corr = (temp_sum - previous)-corr_contrib
           previous = temp_sum
           voxel_weights(ix,iy,iz) = previous
           prev_exact=prev_exact + dble(1.0_8/(nx*ny*nz))
!           write(11, "(6E30.20)") atmsTemp(ix,iy,iz), atmsPlanckRad, totalAbsCoef, 4.0*Pi* atmsPlanckRad * totalAbsCoef*dz(iz), dz(iz), voxel_weights(ix,iy,iz) 
!            write(11, "(9E30.20)") atmsTemp(ix,iy,iz), atmsPlanckRad, totalAbsCoef, 4.0_8*Pi* atmsPlanckRad * totalAbsCoef*dz(iz), dz(iz), voxel_weights(ix,iy,iz), dble( ((iz-1)*nx*ny)+((iy-1)*nx)+ix  )/dble(nx*ny*nz), prev_exact,corr
         end do ! i loop
         col_weights(iy,iz)= previous
!          write(10, "(3I5, A, E30.20, A, E30.20)" ) ix, iy, iz, 'voxel_weights= ', voxel_weights(ix-1,iy,iz), 'col_weights= ', col_weights(iy,iz)
       end do   ! j loop
       level_weights(iz)= previous
!       write(10, "(3I5, A, E30.20, A, E30.20, A, E30.20)" ) ix, iy, iz, 'voxel_weights= ', voxel_weights(ix-1,iy-1,iz), 'col_weights= ', col_weights(iy-1,iz), 'level_weights= ', level_weights(iz)
     end do     ! k loop
    end if
          if (voxel_weights(nx,ny,nz) .gt. 0.0_8) then
               atmsPower = voxel_weights(nx,ny,nz)*(xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8)/dble(nx*ny)  ! [W] total power emitted by atmosphere. Factor of 1000^2 is to convert dx and dy from km to m
               voxel_weights(:,:,:)=voxel_weights(:,:,:)/voxel_weights(nx,ny,nz)     ! normalized
!               do iz = 1, nz
!                  do iy = 1, ny
!                     write(17, "(100E35.25)") voxel_weights(:,iy,iz)
!                  end do
!               end do    
               col_weights(:,:)=col_weights(:,:)/col_weights(ny,nz)
               level_weights(:)=level_weights(:)/level_weights(nz)

               voxel_weights(nx,ny,nz)=1.0_8     ! need this to be 1 for algorithm used to select emitting voxel
               col_weights(ny,nz)=1.0_8
               level_weights(nz)=1.0_8      
          end if
      
!PRINT *, 'level_weights= ', level_weights, 'col_weights= ', col_weights

     totalPower=sfcPower + atmsPower
     if (totalPower .eq. 0.0_8)then
        PRINT *, 'Neither surface nor atmosphere will emitt photons since total power is 0. Not a valid solution'
     end if
     totalFlux=totalPower/((xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8))  ! We want the units to be [Wm^-2] but the x and y positions are in km
!PRINT *, 'atmsPower= ',atmsPower, 'sfcPower= ', sfcPower, ' totalFlux=', totalFlux, ' totalArea=', (xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1)), &
!         ' average column area=', (SUM(dx)/dble(nx))*(SUM(dy)/dble(ny)), (xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))/dble(nx*ny), ' expected radiance=', atmsPlanckRad*(1.0_8-exp(-1.0_8*totalAbsCoef*(zPosition(nz+1)-zPosition(1))))
     if (atmsPower .eq. 0.0_8)then
         atmsPhotons=0
     else   
        atmsPhotons=ceiling(totalPhotons * atmsPower / totalPower)
     end if
  end subroutine emission_weighting

end module monteCarloIllumination
