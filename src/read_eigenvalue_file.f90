subroutine read_eigenvalue_file(filename,skipdump,n,xyz,rho,m,eigenvalues)

implicit none

character(100), intent(in):: filename
logical,intent(inout) :: skipdump
real, allocatable,dimension(:),intent(inout) :: rho,m
real,allocatable,dimension(:,:),intent(inout) :: xyz,eigenvalues
integer,allocatable,dimension(:) :: eigenelement

logical::filexist

skipdump=.false.

inquire(file=filename,exist=filexist)

if(fileexist.eqv..false) then
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
allocate(eigenvalues(3,n))

read(10) (eigenelement(i),i=1,n)

read(10) (xyzmh(1,i), i=1,n)
read(10) (xyzmh(2,i), i=1,n)
read(10) (xyzmh(3,i), i=1,n)

read(10) (eigenvalues(1,i),i=1,n)
read(10) (eigenvalues(2,i),i=1,n)
read(10) (eigenvalues(3,i), i=1,n)
read(10) (rho(i), i=1,n))
read(10) (xyzmh(4,i), i=1,n)
close(10)

print*, 'Eigenvalue File Read complete '

end subroutine read_eigenvalue_file
