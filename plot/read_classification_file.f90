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


SUBROUTINE read_file(filename,n,time,class)
! This subroutine reads a simple binary file
! (To be wrapped by f2py)

!F2PY INTENT(IN) :: filename
!F2PY INTENT(IN) :: n
!F2PY INTENT(OUT) :: time
!F2PY INTENT(OUT) :: class

integer :: m,n
integer,dimension(n) :: class
real*8 :: time

character(100) :: filename

! Check that the specified number of entries (n) 
! matches that in the file (m)

open(10, file=filename, form='unformatted')
read(10) m
print*, 'read_file detects ', m, 'entries in file ',TRIM(filename)

if(m/=n) then
    print*, "ERROR: mismatch between input entry number and read_file entry number"
    print*, m,n
endif

read(10) time
read(10) (class(i),i=1,n)

close(10)

print*, 'Classification File Read complete '

END SUBROUTINE read_file
