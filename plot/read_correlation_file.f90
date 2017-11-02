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


SUBROUTINE read_file(filename,n,nentries,time,class,rclosest,rho,vmag,poten)
! This subroutine reads a simple binary file
! (To be wrapped by f2py)

!F2PY INTENT(IN) :: filename
!F2PY INTENT(IN) :: n
!F2PY INTENT(OUT) :: nentries
!F2PY INTENT(OUT) :: time
!F2PY INTENT(OUT) :: class
!F2PY INTENT(OUT) :: rclosest
!F2PY INTENT(OUT) :: rho
!F2PY INTENT(OUT) :: vmag
!F2PY INTENT(OUT) :: poten

integer :: m,n,nentries
integer,dimension(n) :: class
real,dimension(4,n) :: rclosest
real,dimension(n) :: rho
real,dimension(n) :: vmag
real,dimension(n) :: poten

real*8 :: time

character(100) :: filename

! Check that the specified number of entries (n) 
! matches that in the file (m)

open(10, file=filename, form='unformatted')
read(10) m, nentries
print*, 'read_file detects ', m, 'entries in file ',TRIM(filename)

if(m/=n) then
    print*, "ERROR: mismatch between input entry number and read_file entry number"
    print*, m,n
endif

read(10) time
read(10) (class(i),i=1,n)
read(10) (rclosest(1,i),i=1,n)
read(10) (rclosest(2,i),i=1,n)
read(10) (rclosest(3,i),i=1,n)
read(10) (rclosest(4,i),i=1,n)
read(10) (rho(i),i=1,n)
read(10) (vmag(i),i=1,n)

! If the file contains potential data, then read it as well
if(nentries==7) then
  read(10) (poten(i),i=1,n)
endif

close(10)

print*, 'Correlation File Read complete '

END SUBROUTINE read_file
