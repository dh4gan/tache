PROGRAM tache
  ! This program reads in hydrodynamic simulation files and 
  ! performs tensor classification
  ! Eigenvalues of the tensor determine what environment 
  ! an hydrodynamic fluid element is in
  
  use sphgravdata
  use treedata
  
  implicit none
  
  real,parameter :: pi = 3.141592653
  real, parameter :: twopi = 2.0*pi
  
  integer, parameter:: it_max = 10
  
  integer,allocatable,dimension(:) :: eigenpart
  
  real, allocatable, dimension (:,:,:) :: velocityshear, eigenvectors,eigenvecbin
  real, allocatable, dimension(:,:) :: eigenvalues, eigenbin
  real, allocatable, dimension(:) :: xbin, ybin, zbin
  
  real, dimension (3,3) :: tensor, eigenvec
  real, dimension(3) :: eigen
  
  integer :: ipart, i,j,k,n, start, finish, check, skip
  integer :: it_num, rot_num, ngas, counter
  
  real :: vmean, percentcount
  
  character(1) :: use_octree, scalevelocity
  character(4) :: fileroot
  character(3),dimension(1000) :: num
  character(7),dimension(1000) :: filename,gravfile,potfile
  character(7), dimension(1000) :: eigenfile,vectorfile
  character(100) :: inputfile, filedir
  character(18) :: neighbourfile
  logical :: existneigh
  
  ! Initialise parameters
  
  call initial
      
  ! Now begin tensor classification, file by file
  DO ifile=1,nfiles

     ! *********************
     ! 1. Read in data file
     ! *********************

     call read_dump(filename(ifile),skipdump)

     if(skipdump.eqv..true.) cycle
    

     !********************************************
     ! 2. Compute derivatives and construct tensor
     !********************************************

     call compute_tensor

     !*******************************************
     ! 3. Compute eigenvalues of the tensor
     !*******************************************

     call calc_eigenvalues

     !******************************************
     ! 4. Write data to files
     !******************************************

     call write_eigendata

     ! If we want to split dump into tensor classified components, do it here
     if(splitdump=='y') call splitdump

     ! Deallocate memory ready for the next run
     call deallocate_memory


     ! End of loop over files
  enddo

     !
     ! 4. Write data to file 
     !


     ! Find maximum and minimum values for all coordinates
     xmax = 0.0
     ymax = 0.0
     zmax = 0.0

     allocate(isort(npart))
     allocate(iorig(npart))

     ngas = 0
     do ipart=1,npart
       ! rho(ipart) = xyzmh(4,ipart)/(xyzmh(5,ipart)**3) ! DEBUG LINE
        isort(ipart) = ipart
        iorig(ipart) = ipart

        if(iphase(ipart)==0) ngas = ngas +1

        IF(abs(xyzmh(1,ipart))> xmax) xmax = abs(xyzmh(1,ipart))
        IF(abs(xyzmh(2,ipart))> ymax) ymax = abs(xyzmh(2,ipart))
        IF(abs(xyzmh(3,ipart))> zmax) zmax = abs(xyzmh(3,ipart))

     enddo

     print*,'-----------------------------------------'
     print*, "Maximum values for spatial co-ordinates: "
     print*, "x: ", xmax
     print*, "y: ", ymax
     print*, "z: ", zmax



     
     ! ********************************************
     ! 2. Calculate Velocity Shear Tensor for all particles
     ! ********************************************

     
     ! ****************************************
     ! 3. Calculate Eigenvalues of the Tidal Tensor
     ! ****************************************

   

   


    
     print*, '----------------------'
     print*, 'Writing to file ', TRIM(eigenfile(n))

     ! Write data to file - for now, simple formatted file

     ! Write data to binary file

     ! Prepare special arrays for binary write

     allocate(xbin(ngas))
     allocate(ybin(ngas))
     allocate(zbin(ngas))
     allocate(eigenbin(3,ngas))
     allocate(eigenvecbin(3,3,ngas))
     allocate(eigenpart(ngas))

     counter =1
     do ipart=1,npart
        if(iphase(ipart)/=0) cycle
        if(allocated(iunique)) then
           eigenpart(counter) = iunique(ipart)
        else
           eigenpart(counter) = ipart
        endif
        xbin(counter) = xyzmh(1,ipart)
        ybin(counter) = xyzmh(2,ipart)
        zbin(counter) = xyzmh(3,ipart)
        do k=1,3
           eigenbin(k,counter) = eigenvalues(k,ipart)
           do j=1,3
              eigenvecbin(j,k,counter) = eigenvectors(j,k,ipart)
           enddo
        enddo
        
        counter = counter +1
     enddo

     open(27,file=eigenfile(n), status='unknown',form='unformatted')
     write(27) ngas
     write(27) (eigenpart(i),i=1,ngas)
     write(27) (xbin(i), i=1,ngas)
     write(27) (ybin(i), i=1,ngas)
     write(27) (zbin(i), i=1,ngas)
     write(27) (eigenbin(1,i), i=1,ngas)
     write(27) (eigenbin(2,i), i=1,ngas)     
     write(27) (eigenbin(3,i), i=1,ngas)
     close(27)

     ! Now write the eigenvectors to file
     open(27,file=vectorfile(n),status='unknown', form='unformatted')
     write(27) ngas
     write(27) (eigenpart(i),i=1,ngas)
     write(27) (eigenvecbin(1,1:3,i),i=1,ngas)
     write(27) (eigenvecbin(2,1:3,i),i=1,ngas)
     write(27) (eigenvecbin(3,1:3,i),i=1,ngas)


     deallocate(xbin,ybin,zbin,eigenbin,eigenpart, eigenvecbin)

     call deallocate_memory

     deallocate(eigenvalues,velocityshear, eigenvectors)
     print*, 'Tensor memory deallocated'

  enddo

END PROGRAM tache
