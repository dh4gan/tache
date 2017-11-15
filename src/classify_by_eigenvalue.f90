SUBROUTINE classify_by_eigenvalues(classification, values,)
! This subroutine compares three eigenvalues to a threshold, and determines
! how many eigenvalues exceed that threshold

! Classification = N_above +1  (+1 because fortran can't deal with zero array indices)

!use tachedata
!TODO - make this usable for tache and spiralfind

integer :: k,classification
real,dimension(3) :: values


classification = 1

do k=1,3
   if(values(k)<threshold) then
      classification = classification +1
   endif
enddo

classnum(classification) = classnum(classification)+1

return

END SUBROUTINE classify_by_eigenvalues
