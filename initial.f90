subroutine initial
! Subroutine reads in parameter file
! and sets up simulation
!

character(8),parameter :: paramfile='tache.params'


! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"
  print*, " TACHE: Tensor clAssification of Hydrodynamic Elements"
  print*, "     dh4gan (current version: 1st Nov 2017)    "
  print*, "	                                       	  "
  print*, "-----------------------------------------------"
  print*, " "
  print*, " If SPH read fails, check the following:"
  print*, "    - File endianness"
  print*, "    - Default real size"
  print*, " "
  print*, " For the gfortran compiler, inserting or "
  print*, " removing the following flags should correct"
  print*, " the problem:"
  print*, "    -fconvert=swap"
  print*, "       (swaps endianness during read-in)"
  print*, "    -fdefault-real-8"
  print*, "       (sets default real to double precision)"
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

  do ifile=1,nfiles
     read(20,*) filename(ifile)
  enddo

  ! Generate relevant filenames:

  if(tensorchoice=='tidal') then
     tensorchar = 'T'
  else if(tensorchoice=='velocity') then
     tensorchar='V'
  endif
  
  do ifile=1,nfiles
     write(num, '(I4.3)')i

     write(gravfile(ifile),'("grav",A3)') num
     write(potfile(ifile), '("pot",A3)') num
     
     write(eigenfile(ifile), '("eig",A1,A3)') tensorchar,num
     write(vectorfile(ifile),'("evc",A1,A3)') tensorchar,num    

  enddo

end subroutine initial
