subroutine initial
! Subroutine reads in parameter file
! and sets up simulation
!

use tachedata
implicit none


integer :: ifile
character(3) :: num




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
  read(10,*) threshold ! 
  read(10,*) splitdump ! (y/n)

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

  do ifile=1,nfiles
     write(num, '(I4.3)')ifile

     write(gravfile(ifile),'("grav",A3)') num
     write(potfile(ifile), '("pot",A3)') num
     
     write(eigenfile(ifile), '("eig",A1,A3)') tensorchar,num
     write(vectorfile(ifile),'("evc",A1,A3)') tensorchar,num    

  enddo

end subroutine initial
