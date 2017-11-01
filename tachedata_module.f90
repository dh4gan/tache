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

  integer :: ielement, nelement, neigen

  ! Reals
  real(kind=8) :: threshold

  ! Characters
  character(8),parameter :: paramfile='tache.params'
  character(100) :: listfile, filetype,fileformat,tensorchoice
  character(1) :: splitdump, tensorchar

  ! Parameter decides how neighbours are calculated:
  ! 'g' - particle neighbours found using regular grid
  ! 'o' - neighbours found by building an octree (O(N log N))
  ! 'b' - neighbours found by brute force (O(N^2))

  character(1),parameter :: use_octree_grid = 'g' 

  ! Parameter decides how gravitational force is calculated
  ! 'p' - if potential exists, take its derivative
  ! 'o' - use an octree
  ! 'b' - use brute force (on neighbour list)

  character(1),parameter :: grav_calc_choice = 'b'

  character(100),allocatable,dimension(:) :: filename,gravfile,potfile,eigenfile,vectorfile

! Logicals


		      
!-------------------Array values----------------------------------------
! Integers
  ! Reals
  
  real,allocatable,dimension(:,:) :: eigenvalues
  real,allocatable,dimension(:,:,:) :: eigenvectors
  real, allocatable,dimension(:,:,:) :: tensor
    
end module tachedata
