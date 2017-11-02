SUBROUTINE find_number_entries(filename,n)
!F2PY INTENT(IN) :: filename
!F2PY INTENT(OUT) :: n

character(100) :: filename
integer :: n

print*, 'Finding number of entries in file'

open(10, file=filename, form='unformatted')
read(10) n
print*, 'The file has ', n, 'entries'
close(10)

END SUBROUTINE find_number_entries


SUBROUTINE read_file(filename,n,eigenpart, eigenvectors)
! This subroutine reads a simple binary file
! (To be wrapped by f2py)

!F2PY INTENT(IN) :: filename
!F2PY INTENT(IN) :: n
!F2PY INTENT(OUT) :: eigenpart
!F2PY INTENT(OUT) :: eigenvectors

integer :: m,n
integer,dimension(n) :: eigenpart
real*8,dimension(3,3,n) :: eigenvectors

character(100) :: filename

! Check that the specified number of entries (n) 
! matches that in the file (m)

open(10, file=filename, form='unformatted')
read(10) m
print*, 'read_file detects ', m, 'entries'

if(m/=n) then
    print*, "ERROR: mismatch between input entry number and read_file entry number"
    print*, m,n
endif

read(10) (eigenpart(i),i=1,n)

read(10) (eigenvectors(1,1:3,i),i=1,n)
read(10) (eigenvectors(2,1:3,i),i=1,n)
read(10) (eigenvectors(3,1:3,i),i=1,n)
close(10)

print*, 'Eigenvector File Read complete '

END SUBROUTINE read_file
