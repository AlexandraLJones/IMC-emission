! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Code/multipleProcesses_mpi.f95 $

module multipleProcesses
  !
  ! Module encapsulating all MPI calls needed in the I3RC model. 
  !   These can be replaced with stubs and the code compiled without MPI. 
  !
  implicit none
  include "mpif.h"
!  use MPI
  
  logical :: MasterProc = .true. 
  
  interface sumAcrossProcesses
    module procedure sumAcrossProcesses_Real_Scalar,  &
                     sumAcrossProcesses_Real_1D, sumAcrossProcesses_Double_1D, sumAcrossProcesses_Real_2D, &
                     sumAcrossProcesses_Real_3D, sumAcrossProcesses_Int_3D, sumAcrossProcesses_Real_4D, &
                     sumAcrossProcesses_Real_5D, sumAcrossProcesses_Double_2D, sumAcrossProcesses_Double_3D,&
                     sumAcrossProcesses_Double_4D
  end interface sumAcrossProcesses
contains
  ! -----------------------------------------------------------
  subroutine initializeProcesses(numProcs, thisProcNum)
    implicit none
    include 'mpif.h'
!    use MPI

    integer, intent(out) :: numProcs, thisProcNum
    !
    !Initial MPI calls; how many processors are being used 
    !   and which number is this one? 
    !
    integer :: ierr
!    external :: MPI_Comm_f2c ! added 07/24/12 by Alexandra Jones to fix a bug
!PRINT*, 'entered subroutine initializeProcesses'
    call MPI_INIT(ierr)
!PRINT*, 'error after MPI_INIT', ierr
!    call MPI_COMM_RANK( MPI_COMM_WORLD, MPI_Comm_f2c(thisProcNum), ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, thisProcNum, ierr )
!PRINT*, 'error after MPI_COMM_RANK', ierr
!    call MPI_COMM_SIZE( MPI_Comm_f2c(MPI_COMM_WORLD), numProcs, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numProcs, ierr )
!PRINT*, 'error after MPI_COMM_SIZE', ierr
    
    MasterProc = (thisProcNum == 0)
  end subroutine initializeProcesses
  ! -----------------------------------------------------------
  subroutine synchronizeProcesses
    !
    ! Wait for all processors to get to this point
    ! 
    integer :: ierr

    CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
    
  end subroutine synchronizeProcesses
  ! -----------------------------------------------------------
  subroutine finalizeProcesses
    integer :: ierr
    
    call MPI_FINALIZE(ierr)
  end subroutine finalizeProcesses
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_Scalar(x) 
    !
    ! Add values across all processors
    !
    real, intent(in) :: x
    real             :: sumAcrossProcesses_Real_Scalar
    
    real    :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_Scalar = temp
    
  end function sumAcrossProcesses_Real_Scalar
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_1D(x) 
    !
    ! Add values across all processors
    !
    real, dimension(:), intent(in) :: x
    real, dimension(size(x))       :: sumAcrossProcesses_Real_1D
    
    real, dimension(size(x)) :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_1D(:) = temp(:)
    
  end function sumAcrossProcesses_Real_1D

  ! -----------------------------------------------------------
  function sumAcrossProcesses_Double_1D(x)
    !
    ! Add values across all processors
    !
    real(8), dimension(:), intent(in) :: x
    real(8), dimension(size(x))       :: sumAcrossProcesses_Double_1D

    real(8), dimension(size(x)) :: temp
    integer :: ierr

    call MPI_REDUCE(x, temp, size(x), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Double_1D(:) = temp(:)

  end function sumAcrossProcesses_Double_1D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_2D(x) 
    !
    ! Add values across all processors
    !
    real, dimension(:, :),       intent(in) :: x
    real, dimension(size(x, 1), size(x, 2)) :: sumAcrossProcesses_Real_2D
    
    real, dimension(size(x, 1), size(x, 2)) :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_2D(:, :) = temp(:, :)
    
  end function sumAcrossProcesses_Real_2D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Double_2D(x)
    !
    ! Add values across all processors
    !
    real(8), dimension(:,:), intent(in) :: x
    real(8), dimension(size(x, 1), size(x, 2))       :: sumAcrossProcesses_Double_2D

    real(8), dimension(size(x, 1), size(x, 2)) :: temp
    integer :: ierr

    call MPI_REDUCE(x, temp, size(x), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Double_2D(:,:) = temp(:,:)

  end function sumAcrossProcesses_Double_2D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Double_3D(x)
    !
    ! Add values across all processors
    !
    real(8), dimension(:,:,:), intent(in) :: x
    real(8), dimension(size(x, 1), size(x, 2), size(x, 3))       :: sumAcrossProcesses_Double_3D

    real(8), dimension(size(x, 1), size(x, 2), size(x, 3)) :: temp
    integer :: ierr

    call MPI_REDUCE(x, temp, size(x), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Double_3D(:,:,:) = temp(:,:,:)

  end function sumAcrossProcesses_Double_3D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_3D(x) 
    !
    ! Add values across all processors
    !
    real, dimension(:, :, :), intent(in) :: x
    real, dimension(size(x, 1), size(x, 2), &
                    size(x, 3))          :: sumAcrossProcesses_Real_3D
    
    real, dimension(size(x, 1), size(x, 2), &
                    size(x, 3))             :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_3D(:, :, :) = temp(:, :, :)
    
  end function sumAcrossProcesses_Real_3D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Int_3D(x)
    !
    ! Add values across all processors
    !
    integer, dimension(:, :, :), intent(in) :: x
    integer, dimension(size(x, 1), size(x, 2), &
                    size(x, 3))          :: sumAcrossProcesses_Int_3D

    integer, dimension(size(x, 1), size(x, 2), &
                    size(x, 3))             :: temp
    integer :: ierr

    call MPI_REDUCE(x, temp, size(x), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Int_3D(:, :, :) = temp(:, :, :)

  end function sumAcrossProcesses_Int_3D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_4D(x) 
    real, dimension(:, :, :, :), intent(in) :: x

    real, dimension(size(x, 1), size(x, 2), &
                    size(x, 3), size(x, 4)) :: sumAcrossProcesses_Real_4D
    
    real, dimension(size(x, 1), size(x, 2), &
                    size(x, 3), size(x, 4)) :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_4D(:, :, :, :) = temp(:, :, :, :)
    
  end function sumAcrossProcesses_Real_4D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Double_4D(x)
    !
    ! Add values across all processors
    !
    real(8), dimension(:,:,:,:), intent(in) :: x
    real(8), dimension(size(x, 1), size(x, 2), size(x, 3), size(x, 4))       :: sumAcrossProcesses_Double_4D

    real(8), dimension(size(x, 1), size(x, 2), size(x, 3), size(x, 4)) :: temp
    integer :: ierr

    call MPI_REDUCE(x, temp, size(x), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Double_4D(:,:,:,:) = temp(:,:,:,:)

  end function sumAcrossProcesses_Double_4D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_5D(x) 
    real, dimension(:,:,:,:,:), intent(in) :: x
    real, dimension(size(x,1), size(x,2), size(x,3), &
                    size(x,4), size(x,5)) ::  sumAcrossProcesses_Real_5D

    real, dimension(size(x,1), size(x,2), size(x,3), &
                    size(x,4), size(x,5)) ::  temp

    integer :: ierr

    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_5D(:,:,:,:,:) = temp(:,:,:,:,:)
  end function sumAcrossProcesses_Real_5D
end module multipleProcesses
