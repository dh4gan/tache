SUBROUTINE deallocate_memory
!**************************************************************************
! Subroutine deallocates all memory used for analysing a specific dumpfile
!**************************************************************************
  use sphdata
  use sphneighbourdata
  use tachedata, only: tensorchoice

  implicit none

  print*, 'Deallocating Memory'
  
  deallocate(nneigh,neighb)
  if(allocated(occ)) deallocate(occ)
  if(allocated(n_occ)) deallocate(n_occ)
  
  print*, 'Tree memory deallocated'
  
  deallocate(iphase,isteps,isort)
  if(allocated(iorig)) deallocate(iorig)
  deallocate(xyzmh,vxyzu,rho)	
  deallocate(listpm,spinx,spiny,spinz)
  deallocate(angaddx,angaddy,angaddz)
  deallocate(spinadx,spinady,spinadz)
  
  if(allocated(dgrav)) deallocate(dgrav)
  if(allocated(alphaMM)) deallocate(alphaMM)
  if(allocated(gradh)) deallocate(gradh)
  if(allocated(gradhsoft)) deallocate(gradhsoft)
  
  if(tensorchoice=='tidal') deallocate(gravxyz,poten)
  
  print*, 'SPH memory deallocated'
  
  ! Some MPI arrays might need to be deallocated too
  if(allocated(nelementblocks)) deallocate(nelementblocks)
  if(allocated(iunique)) deallocate(iunique)
  
END SUBROUTINE deallocate_memory
