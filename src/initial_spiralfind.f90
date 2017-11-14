subroutine initial
!
! Subroutine reads in parameter file
! and sets up analysis
!

use tachedata
implicit none

character(100),parameter :: paramfile = 'spiralfind.params'

integer :: ifile
character(100) :: suffix,eigensuffix

! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"
  print*, " SPIRAL FINDING on tensor classified data"
  print*, "     dh4gan (current version: 14th Nov 2017)    "
  print*, "	                                       	  "
  print*, "-----------------------------------------------"
  print*, " "
  
  print*, " "
  print*, " "
  print*, " input parameters in ./",trim(paramfile)
  print*, "-----------------------------------------------"
  print*, " "



  open(10,file=paramfile, status='old')
  read(10,*) listfile ! File containing list of eigenvalue files
  read(10,*) fileformat ! sphNG_wkmr, sphNG_iab
  read(10,*) tensorchoice ! tidal, velocity
  read(10,*) mindist ! Minimum distance from the origin
  read(10,*) D ! Linking length between spiral segments
  read(10,*) xpercentile ! Limit search to this density percentile

  close(10)

  filetype = 'SPH'
  fileformat = 'sphNG_wkmr'

  call checksetup

  ! Setup kernel interpolation tables if this is an SPH run
  if(filetype=='SPH') call ktable

  ! Read listfile and generate array of filenames

  open(20, file=listfile, form='formatted')

  read(20,*) nfiles

  allocate(eigenfile(nfiles))
  do ifile=1,nfiles
     read(20,*) eigenfile(ifile)
  enddo

  ! 

end subroutine initial
