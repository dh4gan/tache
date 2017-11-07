subroutine initial
!
! Subroutine reads in parameter file
! and sets up analysis
!

use tachedata
implicit none


integer :: ifile
character(3) :: num
character(100) :: suffix

! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"
  print*, " TACHE: TensoriAl Classification of Hydrodynamic Elements"
  print*, "     dh4gan (current version: 1st Nov 2017)    "
  print*, "	                                       	  "
  print*, "-----------------------------------------------"
  print*, " "
  
  print*, " "
  print*, " "
  print*, " input parameters in ./",trim(paramfile)
  print*, "-----------------------------------------------"
  print*, " "



  open(10,file=paramfile, status='old')
  read(10,*) listfile ! File containing list of dumps to analyse
  read(10,*) filetype ! SPH, grid or mesh
  read(10,*) fileformat ! sphNG_wkmr, sphNG_iab
  read(10,*) tensorchoice ! tidal, velocity
  read(10,*) splitdumpchoice ! (y/n)
  read(10,*) threshold

  close(10)

  call checksetup

  ! Setup kernel interpolation tables if this is an SPH run
  if(filetype=='SPH') call ktable

  ! Read listfile and generate array of filenames

  open(20, file=listfile, form='formatted')

  read(20,*) nfiles

  allocate(filename(nfiles))
  do ifile=1,nfiles
     read(20,*) filename(ifile)
  enddo

  ! Generate relevant filenames:

  if(tensorchoice=='tidal') then
     tensorchar = 'T'
  else if(tensorchoice=='velocity') then
     tensorchar='V'
  endif
  
  allocate(gravfile(nfiles))
  allocate(potfile(nfiles))
  allocate(eigenfile(nfiles))
  allocate(vectorfile(nfiles))
  allocate(memberfile(nfiles))

  do ifile=1,nfiles

     write(suffix,'(A1,"_",A)'),tensorchar,trim(filename(ifile))

     write(num, '(I3.3)') ifile

     write(gravfile(ifile),'("grav",A3)') num
     write(potfile(ifile), '("pot",A3)') num
     
     write(eigenfile(ifile), '("eig",A)') trim(suffix)
     write(vectorfile(ifile),'("evc",A1,"_",A)') tensorchar,trim(filename(ifile)    )
     write(memberfile(ifile), '("class",A1,"_",A)') tensorchar,trim(filename(ifile))

  enddo

end subroutine initial
