module sphgravdata

  !-----------------------------------------------------------------------
  ! Data module for saving data from sph, grav and pot files
  ! PJC 20/05/2008
  ! DHF 01/09/2009
  ! DHF 26/02/2015 to incorporate modern sphNG output
  !-----------------------------------------------------------------------

  implicit none
  save

  !-------------------Single values---------------------------------------
  ! Integers
  integer :: npart,naccrete,n1,n2,nreassign,nkill,nblocks

  ! Can't allocate sink arrays at runtime for MPI as we don't store nsinktotal
  integer, parameter :: nptmax = 1000 
  
  integer :: numberarray, icount, icountsink
  integer :: nptmass
  integer, dimension(8) :: nums

  integer,parameter :: nsinkmax = 2000 ! Required for contiguous MPI reads



  ! Modern sphNG parameters
  integer, parameter :: igradh = 1 !0 = off, 1 = on
  integer, parameter :: imhd = 0 !0 = off, 1 = on
  integer, parameter :: iexf = 0 !external force
  integer, parameter :: imigrate = 0 !0 = off, 1 = on
  real :: pmrate, rorbitmax

  integer :: iblock !iblock is the block we're currently working on

  integer :: nparttot !total number of particles
  integer, allocatable :: npartblocks(:) !number of particles in each block

  integer :: iuniquemax

  !npart contains number of current data, whether for 1 block or a whole simulation



  ! Reals
  real :: gt, dtmaxdp,gamma,rhozero,RK2,escap,tkin,tgrav,tterm
  real :: anglostx,anglosty,anglostz,specang,ptmassin
  real :: udist,umass,utime,udisti, umassi, utimei

  ! modern sphNG variables
  real :: tmag, Bextx, Bexty, Bextz
  real :: hzero, uzero_n2, hmass, gapfac, sdprof, umagfd
  real :: rorbit_orig, min_rplan, max_rplan, planetesimalmass, coremass_orig, coremass

  ! Characters
  character(100) :: fileident		      

! Logicals
  logical,parameter :: contiguous=.true.
  logical :: tagged
		      
!-------------------Array values----------------------------------------
! Integers
      integer,allocatable,dimension(:) :: tphase,uphase,isteps,listpm,isort
      integer*1,allocatable,dimension(:) :: iphase,iorig

      integer*8, allocatable,dimension(:) :: iunique
! Reals
      real,allocatable,dimension(:) :: spinx,spiny,spinz
      real,allocatable,dimension(:) :: angaddx,angaddy,angaddz
      real,allocatable,dimension(:) :: spinadx,spinady,spinadz
      real,allocatable,dimension(:,:) :: xyzmh,vxyzu
      real*4,allocatable,dimension(:) :: rho,dgrav,poten
      real*4,allocatable,dimension(:,:) :: gravxyz	        

! Extra arrays for modern sphNG dumps

      real*4, allocatable,dimension(:) :: alphaMM,gradh,gradhsoft

!RT data

integer*2,allocatable :: radneigh(:)
real, allocatable :: e(:), rkappa(:), cv(:), rlambda(:), edd(:)
real, allocatable :: force(:,:)
real*4, allocatable :: dlnTdlnP(:), adiabaticgradient(:)
real*4, allocatable :: pressure(:,:), viscosity(:,:), gravity(:,:), radpres(:,:)

      

contains

!---------------------------------------------------------------
!extract subroutine taken directly from rdump and reformatted to fit in
!with the rest of rbin.
  SUBROUTINE extract(tag,rval,rarr,tags,ntags,ierr)    

  character(len=*), intent(in)  :: tag
  real, intent(out) :: rval
  real, intent(in)  :: rarr(:)
  character(len=16), intent(in)  :: tags(:)
  integer, intent(in)  :: ntags
  integer, intent(out) :: ierr
  logical :: matched
  integer :: i

  ierr = 1
  matched = .FALSE.
  rval = 0.0 ! default if not found
  over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
      if (size(rarr) >= i) then
        rval = rarr(i)
          matched = .true.
      endif
      exit over_tags  ! only match first occurrence
    endif
  enddo over_tags
  if (matched) ierr = 0
  if (ierr.NE.0) print "(a)", &
    ' ERROR: could not find '//trim(adjustl(tag))//' in header'
  END SUBROUTINE extract

      end module sphgravdata
