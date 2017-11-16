SUBROUTINE classify_by_eigenvalues(classification,classnum, values,threshold)
!--------------------------------------------------------------------------
! This subroutine compares three eigenvalues to a threshold, and determines
! how many eigenvalues exceed that threshold (N_above +1)
!
! Classification = N_above +1  
! (+1 because fortran can't deal with zero array indices)
!---------------------------------------------------------------------------
implicit none

integer,intent(inout) :: classification
real, intent(in) :: threshold
integer,dimension(4), intent(inout) :: classnum
real,dimension(3) :: values

integer :: k

classification = 1

! Find classification
do k=1,3
   if(values(k)<threshold) then
      classification = classification +1
   endif
enddo

! Add to totals for each class
classnum(classification) = classnum(classification)+1

return

END SUBROUTINE classify_by_eigenvalues
