PROGRAM tache
  ! This program reads in hydrodynamic simulation files and 
  ! performs tensor classification
  ! Eigenvalues of the tensor determine what environment 
  ! an hydrodynamic fluid element is in
    
  use tachedata

  implicit none
      

  integer :: ifile
  logical :: skipdump

  ! Initialise parameters
  
  call initial
      
  ! Now begin tensor classification, file by file
  DO ifile=1,nfiles

     ! *********************
     ! 1. Read in data file
     ! *********************

     call read_dump(ifile,skipdump)

     if(skipdump.eqv..true.) cycle
    
     !********************************************
     ! 2. Compute derivatives and construct tensor
     !********************************************

     call compute_tensor(ifile)

     !*******************************************
     ! 3. Compute eigenvalues of the tensor
     !*******************************************

     call calc_eigenvalues

     !******************************************
     ! 4. Write data to files
     !******************************************

     call write_eigendata(ifile)

     ! If we want to split dump into tensor classified components, do it here
     if(splitdump) call write_splitdump(ifile)

     ! Deallocate memory ready for the next run
     call deallocate_memory

     ! End of loop over files
  enddo
    

print*, 'TACHE run complete'
END PROGRAM tache
