subroutine extract_spiral_elements(filename,skipdump)
!---------------------------------------------
! Written 14/11/17 by dh4gan
! This subroutine loads an eigenvalue file,
! Classifies the elements by a given threshold
! and then extracts one class of elements as "spiral-like"
!---------------------------------------------

use spiraldata

implicit none

character(100), intent(inout) :: filename
logical, intent(inout) :: skipdump

real, allocatable,dimension(:,:) :: xyzfull,eigenvalfull
real,allocatable,dimension(:) :: rhofull,mfull
integer,allocatable,dimension(:) :: class,eigenelementfull
integer,dimension(4) ::classnum

integer :: i,ielement,counter,nfull

real,dimension(3) :: eigensingle

character(10), dimension(4) :: prefixes

prefixes(1) = "cluster"
prefixes(2) = "filament"
prefixes(3) = "sheet"
prefixes(4) = "void"


! Read in full eigenvalue file

call count_eigenvalue_entries(filename,skipdump,nfull)

allocate(eigenelementfull(nfull))
allocate(rhofull(nfull))
allocate(mfull(nfull))
allocate(xyzfull(3,nfull))
allocate(eigenvalfull(3,nfull))

call read_eigenvalue_file(filename,skipdump,nfull,eigenelementfull,xyzfull,rhofull,mfull,eigenvalfull)

 allocate(class(nfull))

   class(:) = -1
   classnum(:) = 0

   do i=1,nfull
      eigensingle(:) = eigenvalfull(:,i)
      CALL classify_by_eigenvalues(class(i), classnum,eigensingle,threshold)
   enddo

print*, 'Classification complete'

print*, 'Spiralfind will operate on ',classnum(spiralclass), ' ',&
     trim(prefixes(spiralclass)), ' elements' 

nelement = classnum(spiralclass)

allocate(xyz(3,nelement))
allocate(eigenelement(nelement))
allocate(eigenvalues(3,nelement))
allocate(mass(nelement))
allocate(rho(nelement))


! Extract elements corresponding to the specified class into main arrays
counter =0
do ielement=1,nfull   

   if(class(ielement)==spiralclass) then
      counter = counter+1
      eigenelement(counter) = eigenelementfull(ielement)
      xyz(:,counter) = xyzfull(:,ielement)
      eigenvalues(:,counter) = eigenvalfull(:,ielement)
      rho(counter) = rhofull(ielement)
      mass(counter) = mfull(ielement)
   endif

enddo

print*, counter, classnum(spiralclass)
print*, maxval(rho), minval(rho)
! Deallocate full arrays
deallocate(xyzfull,rhofull,mfull,eigenvalfull)

end subroutine extract_spiral_elements
