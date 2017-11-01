SUBROUTINE read_neighbours(neighbourfile)
! Subroutine reads in a neighbours file
  use sphgravdata
  use treedata
  implicit none

  integer :: i,j,neighcheck, tolcheck
  character(18)::neighbourfile

 

  open(2, file= neighbourfile,  form = 'UNFORMATTED')

  READ(2)  neighcheck, tolcheck, meanneigh,sdneigh,neighcrit

  if(neighcheck/=neighmax) print*, 'WARNING: mismatch in neighmax: ', neighmax, neighcheck
  if(tolerance/=tolcheck) print*, 'WARNING: mismatch in tolerance: ', tolerance, tolcheck
  READ(2) (nneigh(i), i=1,nelement)
  do i=1,nelement
     READ(2) (neighb(i,j), j=1,nneigh(i))
  enddo
  close(2)

 ! Calculate mean and standard deviation of neighbour counts
 meanneigh = sum(nneigh)/REAL(nelement)
 sdneigh = 0.0

!$OMP PARALLEL &
!$OMP shared(nneigh,meanneigh,nelement)&
!$OMP private(i) &
!$OMP reduction(+:sdneigh)
!$OMP DO SCHEDULE(runtime)
 do i=1,nelement
     sdneigh = sdneigh+(nneigh(i)-meanneigh)**2
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 sdneigh = sqrt(sdneigh/REAL(nelement))

 print*, 'Mean neighbour number is ', meanneigh
 print*, 'Standard Deviation: ', sdneigh
 neighcrit = meanneigh-5.0*sdneigh
     
 print*, 'Clumps created only if neighbour number greater than ', neighcrit


END SUBROUTINE read_neighbours
