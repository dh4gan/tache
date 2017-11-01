module tachedata

  !-----------------------------------------------------------------------
  ! Data module for input parameters to TACHE and tensors
  ! PJC 20/05/2008
  ! DHF 01/09/2009
  ! DHF 26/02/2015 to incorporate modern sphNG output
  !-----------------------------------------------------------------------

  implicit none
  save


  !-------------------Single values---------------------------------------
  ! Integers
  integer :: nfiles

  ! Reals
  real(kind=8) :: threshold

  ! Characters
  character(100) :: listfile, filetype,fileformat,tensorchoice
  character(1) :: splitdump

  character(100),allocatable,dimension(:) :: filename,gravfile,potfile,eigenfile,vectorfile

! Logicals


		      
!-------------------Array values----------------------------------------
! Integers
  ! Reals
  
  real,allocatable,dimension(:,:) :: eigenvalues
  real,allocatable,dimension(:,:,:) :: eigenvectors
    
end module tachedata
