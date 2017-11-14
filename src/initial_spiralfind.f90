subroutine initial
!
! Subroutine reads in parameter file for spiralfind program
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
  read(10,*) spiralclass ! Which class of element is to be analysed
  read(10,*) threshold ! Eigenvalue threshold for classification
  read(10,*) mindist ! Minimum distance from the origin
  read(10,*) D ! Linking length between spiral segments
  read(10,*) xpercentile ! Limit search to this density percentile (%)
  read(10,*) angcrit ! Maximum angle between segments before spiral broken (deg)

  close(10)

  angcrit = angcrit*pi/180.0

  call checksetup

  ! Read listfile and generate array of filenames

  open(20, file=listfile, form='formatted')
  read(20,*) nfiles

  allocate(eigenfile(nfiles))
  do ifile=1,nfiles
     read(20,*) eigenfile(ifile)
  enddo

end subroutine initial
