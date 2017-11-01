SUBROUTINE deallocate_memory
! Subroutine simply deallocates all memory used in loop

use sphgravdata
use treedata

implicit none

   print*, 'Deallocating Memory'

     deallocate(nneigh,neighb,n_occ)
     if(allocated(occ)) deallocate(occ)

     print*, 'Tree memory deallocated'

     deallocate(iphase,isteps,iorig,isort)
     deallocate(xyzmh,vxyzu,rho)	
     deallocate(listpm,spinx,spiny,spinz)
     deallocate(angaddx,angaddy,angaddz)
     deallocate(spinadx,spinady,spinadz)

     if(allocated(dgrav)) deallocate(dgrav)
     if(allocated(alphaMM)) deallocate(alphaMM)
     if(allocated(gradh)) deallocate(gradh)
     if(allocated(gradhsoft)) deallocate(gradhsoft)
     deallocate(gravxyz,poten)

     print*, 'SPH memory deallocated'

     ! Some MPI arrays might need to be deallocated too
     if(allocated(nelementblocks)) deallocate(nelementblocks)
     if(allocated(iunique)) deallocate(iunique)

END SUBROUTINE deallocate_memory
