SUBROUTINE write_neighbours(neighbourfile)
  !**************************************************
  ! Subroutine writes SPH neighbour data to file 
  ! (avoids re-computation every time tache is run)
  !*************************************************

  use sphdata
  use sphneighbourdata
  use tachedata, only: nelement
  
  implicit none
  
  integer :: i,j
  character(18)::neighbourfile
  
  neighbourfile = TRIM(neighbourfile)
  
  print*, 'Writing to file ', neighbourfile
  
  OPEN (2, file=neighbourfile, form='unformatted')
  
  WRITE(2)  neighmax, tolerance, meanneigh,sdneigh,neighcrit
  WRITE(2) (nneigh(i), i=1,nelement)
  do i=1,nelement
     WRITE(2) (neighb(i,j), j=1,nneigh(i))
  enddo
  
  close(2)
  
END SUBROUTINE write_neighbours
