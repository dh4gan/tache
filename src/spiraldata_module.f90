MODULE spiraldata

!-----------------------------------------------------------------
! Module for spiralfind.f90 - contains data and subroutines used
! in spiral spine identification program, which acts on eigenvalue files
!-----------------------------------------------------------------

! Module stores data referring to spirals identified in spiralfind.f90
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

integer :: nelement,nfiles

real,dimension(3) :: xcom

integer, allocatable,dimension(:) :: eigenelement,isort,spiralmember
real,allocatable,dimension(:) :: rho, mass
real,allocatable,dimension(:,:) :: xyz, eigenvalues 

character(100),allocatable,dimension(:) :: eigenfile


! Input parameters
character(100) :: listfile
integer :: spiralclass
real :: mindist, D,xpercentile,threshold,angcrit

contains


subroutine read_parameters
!
! Subroutine reads in parameter file for spiralfind program
! and sets up analysis
!

implicit none

character(100),parameter :: paramfile = 'spiralfind.params'

integer :: ifile

! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"
  print*, " SPIRAL FINDING on tensor classified data"
  print*, "     dh4gan (current version: 14th Nov 2017)    "
  print*, "	                                       	  "
  print*, "-----------------------------------------------"
  print*, " "
  
  print*, " "
  print*, " "
  print*, " input parameters in ./",trim(paramfile)
  print*, "-----------------------------------------------"
  print*, " "

  open(10,file=paramfile, status='old')
  read(10,*) listfile ! File containing list of eigenvalue files
  read(10,*) spiralclass ! Which class of element is to be analysed
  read(10,*) threshold ! Eigenvalue threshold for classification
  read(10,*) mindist ! Minimum distance from the origin
  read(10,*) D ! Linking length between spiral segments
  read(10,*) xpercentile ! Limit search to this density percentile (%)
  read(10,*) angcrit ! Maximum angle between segments before spiral broken (deg)

  close(10)

  angcrit = angcrit*pi/180.0

  call check_parameters

  ! Read listfile and generate array of filenames

  open(20, file=listfile, form='formatted')
  read(20,*) nfiles

  allocate(eigenfile(nfiles))
  do ifile=1,nfiles
     read(20,*) eigenfile(ifile)
  enddo

end subroutine read_parameters

subroutine check_parameters
!
! Subroutine checks the inputs to spiralfind are correct
!

implicit none

logical :: existlist

! Does the filename list exist?

inquire(file=listfile,exist= existlist)

if(existlist.eqv..false.) then
   print*, 'ERROR in spiralfind: List of filenames ',listfile,' not found'
   STOP
endif

if(spiralclass <0 .or.spiralclass>4) then
   print*, 'ERROR in spiralfind: element classification of choice (spiralclass) out of range (1,4)'
   STOP
endif

if(xpercentile<=0.0 .or. xpercentile>=100.0) then
   print*, 'ERROR in spiralfind: percentile choice out of range (0,100)'
   STOP
endif


end subroutine check_parameters


!---------------------------------------------------
!+
! Calculates the position, velocity of centre of mass
!+
!----------------------------------------------------
subroutine calc_centre_of_mass

implicit none

integer :: k,ielement
real :: totmass

xcom(:) = 0.0

totmass = 0.0

do ielement=1,nelement

   do k=1,3
      xcom(k) = xcom(k) + xyz(k,ielement)*mass(ielement)
      totmass = totmass + mass(ielement)
   enddo

enddo

if(totmass>1.0e-30) then
   xcom(:) = xcom(:)/totmass
else
   xcom(:) = 0.0
endif

print '(A,3(es10.2,1X))', 'Centre of Mass: Position', xcom

end subroutine calc_centre_of_mass

subroutine sort_by_density
  !
  ! subroutine sorts elements in order of increasing density   
  

  implicit none

  integer :: ielement
  real,allocatable,dimension(:) :: rhohold

  print*, 'Sorting elements by density'

  allocate(isort(nelement))
     
   do ielement=1,nelement
      isort(ielement) = ielement
   enddo
     
   rhohold = rho

   print*,minval(rho),maxval(rho)
   print*, rho(1:10)
   print*, minval(isort),maxval(isort)
   CALL sort2(nelement,rhohold,isort,nelement)
   print*, minval(rho),maxval(rho)
   print*, rho(1:10)
   print*, minval(isort),maxval(isort)

   deallocate(rhohold)
   print*, 'Elements sorted by Density'
   print*, "-----------------------------------------------"
   
end subroutine sort_by_density


subroutine apply_percentile_cut
!
! Removes all but the top x% of sorted elements from the analysis
!

implicit none 

integer :: i,ielement,ipercentile

  ! Only consider the top x percent

   print '(a,f5.1,a)', 'Only considering the top ',xpercentile,' density percentile'

   ipercentile = int(xpercentile*real(nelement)/100.0)

   print'(a,I7,a)', 'This constitutes ', ipercentile, ' elements'
   
   do ielement= ipercentile,nelement
      i = isort(nelement-ielement+1)
      spiralmember(i)=-1
   enddo




end subroutine apply_percentile_cut


!----------------------------------------------------
!+
! Calculates separation of two elements
!+
!-----------------------------------------------------
subroutine calc_separation(i,j,r)

implicit none

integer,intent(in) :: i,j
real,intent(inout) :: r
integer :: k

 r = 0.0
 do k=1,3
    r = r + (xyz(k,i)-xyz(k,j))**2
 enddo

 r = sqrt(r)

end subroutine calc_separation


!----------------------------------------------------
!+
! Calculates distance to the origin
!+
!-----------------------------------------------------
subroutine calc_origin_distance(i, r_origin,xcom)

implicit none

integer, intent(in) :: i
real, intent(inout) :: r_origin
real,intent(in) :: xcom(3)

integer :: k
 r_origin = 0.0
 do k=1,3
    r_origin = r_origin + (xyz(k,i)-xcom(k))**2
 enddo

 r_origin = sqrt(r_origin)

end subroutine calc_origin_distance

!-----------------------------------------------------
!+
! Calculate angle between current spiral line segment and next potential segment
!+
!-----------------------------------------------------
subroutine calc_position_angle(i,j, ispiral,pos_ang)

implicit none
integer,intent(in) :: i,j,ispiral
integer :: k, lastseg,thisseg
real, intent(out) :: pos_ang
real, dimension(3) :: rcurrent,rnext
real :: rcurrmag, rnextmag

pos_ang = 0.0

! Previous line segment
thisseg = spirals(ispiral)%nseg
!print*, thisseg
! If this is the first line segment, skip the calculation
if(thisseg ==1) return

! Otherwise find the current line segment's beginning point
lastseg = thisseg-1

rcurrent(:) = 0.0
rcurrmag = 0.0
do k=1,3
   rcurrent(k) = spirals(ispiral)%r(k,thisseg)-spirals(ispiral)%r(k,lastseg)
   rcurrmag = rcurrmag + rcurrent(k)*rcurrent(k)
enddo
rcurrent(:) = rcurrent(:)/rcurrmag

! Next (potential) line segment

rnext(:) = 0.0
rnextmag = 0.0
do k=1,3
   rnext(k) = xyz(k,j) - spirals(ispiral)%r(k,thisseg)
   rnextmag = rnextmag + rnext(k)*rnext(k)
enddo

rnext(:) = rnext(:)/rnextmag

! Get angle between them

pos_ang = 0.0
do k=1,3
   pos_ang = pos_ang + rcurrent(k)*rnext(k)
enddo

pos_ang = acos(pos_ang)
!print*, rcurrent(:), rnext(:), pos_ang

end subroutine calc_position_angle


!----------------------------------------------------
!+
! Adds position of an element  to the spiral
!+
!-----------------------------------------------------
subroutine add_spiral_segment(ispiral, ielement)

implicit none

integer,intent(in) :: ispiral,ielement

 spirals(ispiral)%nseg = spirals(ispiral)%nseg+1
 spirals(ispiral)%r(1:3,spirals(ispiral)%nseg) = xyz(1:3,ielement)

end subroutine add_spiral_segment

!-----------------------------------------------------
!+
! Writes spiral segment to file
!+
!------------------------------------------------------

subroutine write_spiral_data(ispiral,dumpfile)

implicit none

integer,intent(in) :: ispiral
integer :: iseg
character(len=4) :: spiralnum
character(len=100) :: spiralfile
character(len=100) :: dumpfile


write(spiralnum,'(I4.4)') ispiral

spiralfile = trim(dumpfile)//'_spiral_'//trim(spiralnum)//'.dat'

print*, 'Writing to file ', trim(spiralfile),' : ',spirals(ispiral)%nseg, ' segments'

open(10, file=spiralfile, status='unknown')

do iseg = 1,spirals(ispiral)%nseg
   write(10,*) spirals(ispiral)%r(1:3,iseg)
enddo

close(10)

end subroutine write_spiral_data

subroutine write_spiralmember_data(dumpfile)

implicit none
integer :: i
character(len=100) :: dumpfile,memberfile

memberfile = trim(dumpfile)//'_spiralmembers.dat'

OPEN(10,file=memberfile)
write(10,*) (eigenelement(i),i=1,nelement)
write(10,*) (spiralmember(i),i=1,nelement)
close(10)

end subroutine write_spiralmember_data

!----------------------------------------
!+
! Deallocates spiraldata memory ready for next file
!+
!----------------------------------------
subroutine deallocate_spiraldata_memory

deallocate(eigenelement,eigenvalues)
deallocate(xyz,rho,mass,isort)
deallocate(spiralmember)


end subroutine deallocate_spiraldata_memory

END MODULE spiraldata
