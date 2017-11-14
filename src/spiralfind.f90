program spiralfind

use spiraldata

implicit none

!---------------------------
!+
! Written 30/1/17 by dh4gan
! Code attempts to detect spiral structure in SPH data
! Uses tensor classified data (Forgan et al 2016)
!
! Spiral spine elucidated by friends of friends approach
!
!--------------------------



integer,allocatable,dimension(:) :: tested
character(100),parameter :: paramfile = 'spiralfind.params'

integer :: skip,check,jmax,ispiral

integer :: i,ipart,jpart,n,k
integer :: start,finish

real :: rhomax,percent
real :: ri,rj,rjmax

real,dimension(3) :: xcom,vcom

logical :: skipdump

! Initialise parameters
call initial

do ifile=1,nfile

   ! ********************************************************
   ! 1. Read in eigenvalue file (x,y,z,rho,mass,eigenvalues)
   !    And extract spiral-like elements
   ! *********************************************************
   call extract_spiral_elements(eigenfile(ifile),skipdump,spiralclass,threshold)

   if(skipdump.eqv..true.) cycle

   !**************************************
   ! 2.	Order elements by their density
   !**************************************

   call sort_by_density(nelement,isort)

   call apply_percentile_cut(xpercentile)

   ispiral = 0
   percentcount = 0.0

   do ielement=1,nelement

      call element_percent_complete(ielement,nelement,percentcount,increment)

        !*****************************************************
        ! 2. Pick a dense particle (i) a sufficient distance from the origin
        ! to pick the start of a spiral
        !*****************************************************

        i = isort(npart-ipart+1)

        ! If we've already tested this particle, then cycle back
        if(spiralmember(i)/=0) cycle

        ! If we're too close to the origin, cycle
        call calc_origin_distance(i,ri,xcom)
        if(ri<mindist) then
           spiralmember=-1
           cycle
        endif

        ! Otherwise, start a new spiral
        
        ispiral = ispiral + 1
        spirals(ispiral)%nseg = 0

        spiralmember(i) = ispiral

        print*, 'Beginning spiral ',ispiral
        print '(a,3(es10.2,1X))', 'Location: ',xyz(1:3,i)
        call add_spiral_segment(ispiral, i)

        ! Begin loop to define spiral ispiral

        jmax = 1
        do while (jmax>0)


        !*****************************************************
        ! 3. Search neighbours of i (within some search radius D)
        !    Flag all searched particles
        !    Find j: rho_j is maximum, and |rj-origin| > |ri-origin|
        !    Also demand that change in velocity direction between i and j is small (angcrit)
        !*****************************************************

        rhomax = -1.0e30
        jmax = -1
        do jpart=1,npart

           ! Ignore all elements already tested
           if(spiralmember(jpart)/=0) cycle

           ! If element within linking length D, it's in spiral ispiral
           call calc_separation(i,jpart,rpart)
           if(rpart< D) then

              spiralmember(jpart)=ispiral
              call calc_origin_distance(jpart,rj,xcom)

              ! Calculate angle between velocity vectors of i and jpart
              call calc_position_angle(i,jpart,ispiral,pos_ang)
                 
              if(rj > ri .and. rho(jpart)>rhomax .and. pos_ang < angcrit) then
                 rhomax = rho(jpart)
                 jmax = jpart
                 rjmax = rj
              endif
           endif

        enddo


        ! If we have a successful match, add it to the spiral

        if(jmax>0) then         
              
!           print*, 'Next link: ', jmax
!           print '(a,4(es10.2,1X))', 'Location, density: ',xyzmh(1:3,i), rhomax
           call add_spiral_segment(ispiral, jmax)
           ! Now use this new position to find the next point
           i = jmax
           ri = rjmax

        endif

     enddo ! End loop to define spiral

     ! Write spiral data to file (only write spirals with more than two links)
     
     if(spirals(ispiral)%nseg>2) then
        call write_spiral_data(ispiral,filename(n))
     else
        print*, 'Spiral ', ispiral, ' has insufficient segments'
        print*, 'Trying again'

        ! Remove all particles that were once part of this spiral from analysis

        do jelement=1,nelement
           if(spiralmember(jelement) == ispiral) spiralmember(jelement)=-1
        enddo
        ispiral = ispiral-1

     endif

     !if(ispiral>0) stop
  enddo  ! End loop over all particles

! Write spiral membership data to file

OPEN(10,file='spiral_members')
write(10,*) (eigenelement(i),i=1,nelement)
write(10,*) (spiralmember(i),i=1,nelement)
close(10)

! Deallocate arrays ready for next loop
call deallocate_memory
deallocate(spiralmember)

enddo ! End loop over all files


end program spiralfind

!----------------------------------------------------
!+
! Calculates separation of elements
!+
!-----------------------------------------------------
subroutine calc_separation(i,j,r)

use spiraldata,only:xyz

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

use spiraldata,only:xyz

implicit none

integer, intent(in) :: i
real, intent(inout) :: r_origin
real,intent(in) :: xcom(3)

integer :: k
 r_origin = 0.0
 do k=1,3
    r_origin = r_origin + (xyzmh(k,i)-xcom(k))**2
 enddo

 r_origin = sqrt(r_origin)

end subroutine calc_origin_distance


!-----------------------------------------------------
!+
! Calculate angle between current spiral line segment and next potential segment
!+
!-----------------------------------------------------
subroutine calc_position_angle(i,j, ispiral,pos_ang)

use spiraldata

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
use spiraldata
implicit none

integer,intent(in) :: ispiral,ipart

 spirals(ispiral)%nseg = spirals(ispiral)%nseg+1
 spirals(ispiral)%r(1:3,spirals(ispiral)%nseg) = xyz(1:3,ipart)
! spirals(ispiral)%h(spirals(ispiral)%nseg) = xyzmh(5,ipart)


end subroutine add_spiral_segment


subroutine write_spiral_data(ispiral,dumpfile)
use spiraldata

implicit none

integer,intent(in) :: ispiral
integer :: iseg
character(len=4) :: spiralnum
character(len=100) :: spiralfile
character(len=7) :: dumpfile


write(spiralnum,'(I4.4)') ispiral

spiralfile = trim(dumpfile)//'_spiral_'//trim(spiralnum)//'.dat'

print*, dumpfile
print*, 'Writing to file ', spiralfile

open(10, file=spiralfile, status='unknown')

do iseg = 1,spirals(ispiral)%nseg
   write(10,*) spirals(ispiral)%r(1:3,iseg)
enddo

close(10)

end subroutine write_spiral_data
