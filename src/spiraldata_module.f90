MODULE spiraldata
! Module stores data referring to clumps produced in spiralfind.f90
! These clumps are found for a series of annuli (segments), and then stitched into a spiral
! We derive types for both the clumps, and the spirals we produce

implicit none

save

real,parameter :: pi = 3.1415926585
integer, parameter :: nsegmax = 1000
integer, parameter :: nspiralmax =100

type sphspiral

integer :: nseg,num
real, dimension(3,nsegmax) :: r
!real, dimension(nsegmax) :: h
character(100) :: spiralfile

end type sphspiral

type(sphspiral),dimension(nspiralmax) :: spirals(nspiralmax)

integer :: nelement
integer, allocatable,dimension(:) :: eigenelement,isort,spiralmember
real,allocatable,dimension(:) :: rho, mass
real,allocatable,dimension(:,:) :: xyz, eigenvalues 


! TODO - add input parameters here

END MODULE spiraldata
