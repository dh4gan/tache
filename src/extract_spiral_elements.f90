subroutine extract_spiral_elements(filename,skipdump,spiralclass,threshold)
!---------------------------------------------
! Written 14/11/17 by dh4gan
! This subroutine loads an eigenvalue file,
! Classifies the elements by a given threshold
! and then extracts one class of elements as "spiral-like"
!---------------------------------------------

use spiraldata

implicit none

prefixes(1) = "cluster"
prefixes(2) = "filament"
prefixes(3) = "sheet"
prefixes(4) = "void"


real, allocatable,dimension(:,:) :: xyzfull,eigenvalfull
real,allocatable,dimension(:) :: rhofull,mfull

! Read in full eigenvalue file

call read_eigenvalue_file(filename,skipdump,nfull,xyzfull,rhofull,mfull,eigenvalfull)

 allocate(class(nfull))

   class(:) = -1
   classnum(:) = 0

   do i=1,nfull
      eigensingle(:) = eigenvalfull(:,ielement)
      CALL classify_by_eigenvalues(class(ielement), eigensingle,threshold)

      ! Only want particles of given class

      if(class(ielement)/=spiralclass) tested(ielement)=-1

   enddo

print*, 'Classification complete'

print*, 'Spiralfind will operate on ',classnum(spiralclass), ' ',trim(prefixes(spiralclass), ' elements' 

nelement = classnum(spiralclass)

allocate(xyz(3,nelement))
allocate(eigenvalues(3,nelement))
allocate(mass(nelement))
allocate(rho(nelement))


! Extract elements corresponding to the specified class into main arrays
counter =0
do ielement=1,nfull   

   if(class(ielement)==spiralclass) then
      counter = counter+1
      xyz(:,counter) = xyzfull(:,ielement)
      eigenvalues(:,counter) = eigenvalfull(:,ielement)
      rho(counter) = rhofull(ielement)
      mass(counter) = mfull(ielement)
   endif

enddo

! Deallocate full arrays
deallocate(xyzfull,rhofull,mfull,eigenvalfull)

end subroutine extract_spiral_elements
