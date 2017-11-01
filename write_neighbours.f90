
SUBROUTINE write_neighbours(neighbourfile)

use sphgravdata
use treedata
implicit none

 integer :: i,j
 character(18)::neighbourfile

neighbourfile = TRIM(neighbourfile)

print*, 'Writing to file ', neighbourfile

 OPEN (2, file=neighbourfile, form='unformatted')

 WRITE(2)  neighmax, tolerance, meanneigh,sdneigh,neighcrit
 WRITE(2) (nneigh(i), i=1,npart)
 do i=1,npart
     WRITE(2) (neighb(i,j), j=1,nneigh(i))
 enddo

 close(2)

END SUBROUTINE write_neighbours
