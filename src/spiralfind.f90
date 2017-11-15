program spiralfind

!---------------------------
!+
! Written 30/1/17 by dh4gan
! Code attempts to detect spiral structure in SPH data
! Uses tensor classified data (Forgan et al 2016)
!
! Spiral spine elucidated by friends of friends approach
! All subroutines used contained in spiraldata_module.f90
!+
!--------------------------

use spiraldata

implicit none

integer :: i,k
integer :: ielement,jelement,ifile
integer :: jmax,ispiral

real :: rhomax,percent,increment,pos_ang
real :: ri,rj,rjmax,rsep

logical :: skipdump

! Initialise parameters
call initial

do ifile=1,nfiles

   ! ********************************************************
   ! 1. Read in eigenvalue file (x,y,z,rho,mass,eigenvalues)
   !    And extract spiral-like elements
   ! *********************************************************
   call extract_spiral_elements(eigenfile(ifile),skipdump,spiralclass,threshold)

   if(skipdump.eqv..true.) cycle

   ! Find centre of mass for appropriate system origin
   call calc_centre_of_mass

   !**************************************
   ! 2.	Order elements by their density
   !**************************************

   call sort_by_density(nelement,isort)

   call apply_percentile_cut(xpercentile)

   ispiral = 0
   percent = 0.0

   !***************************************
   ! 3. Begin friends of friends algorithm
   !***************************************
   do ielement=1,nelement

      call element_percent_complete(ielement,nelement,percent,increment)

        !*****************************************************
        ! 3a. Pick a dense particle (i) a sufficient distance from the origin
        ! to pick the start of a spiral
        !*****************************************************

        i = isort(nelement-ielement+1)

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
        ! 3b. Search neighbours of i (within some search radius D)
        !    Flag all searched particles
        !    Find j: rho_j is maximum, and |rj-origin| > |ri-origin|
        !    Also demand that change in velocity direction 
        !    between i and j is small (angcrit)
        !*****************************************************

        rhomax = -1.0e30
        jmax = -1
        do jelement=1,nelement

           ! Ignore all elements already tested
           if(spiralmember(jelement)/=0) cycle

           ! If element within linking length D, it's in spiral ispiral
           call calc_separation(i,jelement,rsep)
           if(rsep< D) then

              spiralmember(jelement)=ispiral
              call calc_origin_distance(jelement,rj,xcom)

              ! Calculate angle between velocity vectors of i and jpart
              call calc_position_angle(i,jelement,ispiral,pos_ang)
                 
              if(rj > ri .and. rho(jelement)>rhomax .and. pos_ang < angcrit) then
                 rhomax = rho(jelement)
                 jmax = jelement
                 rjmax = rj
              endif
           endif

        enddo


        ! If we have a successful match, add it to the spiral

        if(jmax>0) then                      
           call add_spiral_segment(ispiral, jmax)
           ! Now use this new position to find the next point
           i = jmax
           ri = rjmax
        endif

     enddo ! End loop to define spiral

     ! Write spiral data to file (only write spirals with more than two links)
     
     if(spirals(ispiral)%nseg>2) then
        call write_spiral_data(ispiral,eigenfile(ifile))
     else
        print*, 'Spiral ', ispiral, ' has insufficient segments'
        print*, 'Trying again'

        ! Remove all particles that were once part of this spiral from analysis

        do jelement=1,nelement
           if(spiralmember(jelement) == ispiral) spiralmember(jelement)=-1
        enddo
        ispiral = ispiral-1

     endif

  enddo  ! End loop over all particles

! Write spiral membership data to file

call write_spiralmember_data(eigenfile(ifile))

! Deallocate arrays ready for next loop
call deallocate_memory
deallocate(spiralmember)

enddo ! End loop over all files


end program spiralfind


