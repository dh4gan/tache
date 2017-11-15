subroutine read_eigenvalue_file(filename,skipdump,n,xyz,rho,mass,eigenvalues)

implicit none

character(100), intent(in):: filename
logical,intent(inout) :: skipdump
integer, intent(inout) :: n
real, allocatable,dimension(:),intent(inout) :: rho,mass
real,allocatable,dimension(:,:),intent(inout) :: xyz,eigenvalues
integer,allocatable,dimension(:) :: eigenelement

integer :: i
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
read(10) n

print*, 'There are ', n, 'entries'

allocate(eigenelement(n))

allocate(xyz(3,n))
allocate(rho(n))
allocate(mass(n))
allocate(eigenvalues(3,n))

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

print*, 'Eigenvalue File Read complete '

end subroutine read_eigenvalue_file
