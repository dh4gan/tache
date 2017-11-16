subroutine count_eigenvalue_entries(filename,skipdump,n)
!
!
!
!

implicit none

character(100),intent(in) :: filename
logical,intent(inout) :: skipdump
integer,intent(inout) :: n

logical :: fileexist

skipdump=.false.

inquire(file=filename,exist=fileexist)

if(fileexist.eqv..false.) then
   print*, 'Skipping missing file ',trim(filename)
   skipdump=.true.
   return
endif

print*, 'Reading eigenvalue file ',trim(filename)
open(10, file=filename, form='unformatted')
read(10) n
close(10)

print*, 'There are ', n, 'entries'

return
end subroutine count_eigenvalue_entries


subroutine read_eigenvalue_file(filename,skipdump,n,eigenelement,xyz,rho,mass,eigenvalues)

implicit none

character(100), intent(in):: filename
logical,intent(inout) :: skipdump
integer, intent(inout) :: n
real,dimension(n) :: rho,mass
real,dimension(3,n) :: xyz,eigenvalues
integer,dimension(n) :: eigenelement

integer :: i,m
logical::fileexist

skipdump=.false.

inquire(file=filename,exist=fileexist)

if(fileexist.eqv..false.) then
   print*, 'Skipping missing file ',trim(filename)
   skipdump=.true.
   return
endif

print*, 'Reading eigenvalue file ',trim(filename)
open(10, file=filename, form='unformatted')
read(10) m

if(m/=n) then
   print*, 'ERROR in read_eigenvalue_file: mismatch in entry counts'
   print*, m,n
endif

read(10) (eigenelement(i),i=1,n)

read(10) (xyz(1,i), i=1,n)
read(10) (xyz(2,i), i=1,n)
read(10) (xyz(3,i), i=1,n)

read(10) (eigenvalues(1,i),i=1,n)
read(10) (eigenvalues(2,i),i=1,n)
read(10) (eigenvalues(3,i), i=1,n)
read(10) (rho(i), i=1,n)
read(10) (mass(i), i=1,n)
close(10)

print*, 'Max, min density: ',maxval(rho),minval(rho)

print*, 'Eigenvalue File Read complete '

end subroutine read_eigenvalue_file
