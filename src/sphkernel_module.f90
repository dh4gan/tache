module sphkerneldata

  !-----------------------------------------------------------------------
  ! Data module for SPH kernels
  ! 
  ! DHF 17/10/2009
  !-----------------------------------------------------------------------
  
  implicit none
  save
  
  !-------------------Single values---------------------------------------
  ! Integers
  integer, parameter :: itable=40001
  integer :: index,index1
  integer :: igphi,iptintree,iptsoft,isoft 
  ! Reals
  real :: dvtable,radkernel,part1kernel,part2kernel,v2max
  real :: selfnormkernel, cnormk, part1potenkernel,part2potenkernel
  real :: poteni, phi,dfmassdx,dphi,dphiti,dpotdh
  real :: dx,dy,dz,pmassj,fm,dhmean
  real :: dgrwdx, dwdx, dxx, grwtij, wtij, hmean,hmean21,hmean31,hmean41
  real :: rij, rijgrav,rij2grav,sep,xmasj
  real :: rij1, rij2,v2, hi,hj, dfptdx
  
  real,parameter:: psoft = 0.01
  
  !-------------------Array values----------------------------------------
  ! Integers
  
  ! Reals
  real,allocatable,dimension(:) :: wij,grwij,fmass,fpoten,dphidh
  
end module sphkerneldata
