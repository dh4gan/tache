subroutine checksetup
!
! Subroutine checks the inputs to spiralfind are correct
!

use spiraldata
implicit none

logical :: existlist

! Does the filename list exist?

inquire(file=listfile,exist= existlist)

if(existlist.eqv..false.) then
   print*, 'ERROR in spiralfind: List of filenames ',listfile,' not found'
   STOP
endif

if(spiralclass <0 .or.spiralclass>4) then
   print*, 'ERROR in spiralfind: element classification of choice (spiralclass) out of range (1,4)'
   STOP
endif

if(xpercentile<=0.0 .or. xpercentile>=100.0) then
   print*, 'ERROR in spiralfind: percentile choice out of range (0,100)'
   STOP
endif


end subroutine checksetup
