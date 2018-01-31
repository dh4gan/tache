SUBROUTINE mark_children(ipar,used)
!**************************************************************
! Written 01/03/2016 by dh4gan
! Subroutine marks all the descendants of node ipar as 'used'
! Subroutine begins by constructing a list of all descendants
! Then marks all nodes in the list
!**************************************************************

use sphneighbourdata

implicit none
integer, dimension(n_node) :: used, mark
integer :: ipar,ichild,jchild, markcounter,marktotal

if(n_child(ipar)==0) return

mark(:) = 0
marktotal = 0

! Construct the list of descendants
! Begin the list with the children of ipar

do ichild = 1,n_child(ipar)
   marktotal = marktotal+1
   mark(marktotal) = child(ipar,ichild)
enddo

! Now loop over the descendants of the eight children

! marktotal holds the total number of nodes in the list
! markcounter points to where the selection algorithm is in the list

! Once they are equal, this means all the descendant nodes have been found

markcounter = 1

do while(markcounter < marktotal)
   ichild = mark(markcounter)
   
   if(n_child(ichild)/=0) then
      do jchild = 1,n_child(ichild)
         marktotal = marktotal +1
         mark(marktotal) =child(ichild,jchild)
      enddo
   endif

   markcounter = markcounter +1

enddo

! Now go through the list and set them all as used

!print*, 'Node ', ipar, 'has ',marktotal, ' children ', markcounter
do ichild = 1,marktotal
   used(mark(ichild))=1
enddo



END SUBROUTINE mark_children
