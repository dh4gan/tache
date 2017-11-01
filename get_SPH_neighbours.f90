subroutine get_SPH_neighbours(ifile)
!
! Subroutine either reads SPH neighbour data from file, or
! computes it via a regular grid
!

use tachedata
use sphdata, only: poten,iphase,xyzmh,isort,iorig
use sphneighbourdata

implicit none

integer,intent(in) :: ifile

integer :: j
real :: hmean
character(100) :: neighbourfile
logical :: existneigh

!****************************************************************
! 1. Test to see if neighbour file exists: if it does, no need to
! calculate neighbours
!*****************************************************************

neighbourfile ="neighbours_"//TRIM(filename(ifile))

inquire(file=neighbourfile,exist = existneigh)

if(existneigh.eqv..true.) then

   print*, 'Neighbour file',neighbourfile, ' found: reading'
   call read_neighbours(neighbourfile)

else
   !***************************************************************
   ! 2. If no neighbour file found, then we must generate the list
   !***************************************************************

   allocate(nneigh(nelement))
   allocate(neighb(nelement,neighmax))
   nneigh(:) = 0
   neighb(:,:) = 0

    ! Find maximum and minimum values for all coordinates
   rmax(:) = 0.0

   allocate(isort(nelement))
   allocate(iorig(nelement))

   do ielement=1,nelement
      ! rho(ielement) = xyzmh(4,ielement)/(xyzmh(5,ielement)**3) ! DEBUG LINE
      isort(ielement) = ielement
      iorig(ielement) = ielement
     
      if(iphase(ielement)==0) neigen = neigen +1
     
      do j=1,3
         IF(abs(xyzmh(j,ielement))> rmax(j)) xmax = abs(xyzmh(1,ielement))
      enddo

      hmean = hmean + xyzmh(5,ielement)
     
   enddo

   hmean = hmean/REAL(neigen)
  
   xmax = rmax(1)
   ymax = rmax(2)
   zmax = rmax(3)

   print*,'-----------------------------------------'
   print*, "Maximum values for spatial co-ordinates: "
   print*, "x: ", xmax
   print*, "y: ", ymax
   print*, "z: ", zmax
   
   print*, 'Mean smoothing length: ',hmean
  

   if(use_octree_grid=='o') THEN

      print*, "-----------------------------------------------"
      print*, 'Building octree'
      CALL make_octree(filename(ifile))
      
      !	Use octree to find neighbours    
  
      print*, "-----------------------------------------------"
      print*, 'Creating Neighbour Lists, nelement: ',nelement		
      CALL neighbours_octree(filename(ifile))
      
      print*, 'Neighbour lists created'
      print*, "-----------------------------------------------"			
      
     ELSE IF (use_octree_grid=='g') THEN

        print*, "-----------------------------------------------"
        print*, 'Building regular grid'
        CALL make_grid_sphdata

        !	Use grid to find neighbours    
  
        print*, "-----------------------------------------------"
        print*, 'Creating Neighbour Lists from grid, nelement: ',nelement		
        CALL neighbours_grid(filename(ifile))

        print*, 'Neighbour lists created'
        print*, "-----------------------------------------------"

     ELSE
        print*, 'Finding neighbours by brute force'
        CALL neighbours_brute(filename(ifile))
     ENDIF

     ! Write the neighbour data to file
     call write_neighbours(neighbourfile)

     if(tensorchoice=='tidal') then
        print*, 'Computing tidal tensor - require gravitational force'
        
        ! TODO - check if gravitational force file exists
    
        IF(use_octree_grid=='o' .and. grav_calc_choice=='o') THEN

           print*, 'Calculating Gravitational Force and Potential Using Octree'
           call calc_grav_tree

        ELSE IF(grav_calc_choice=='p' .and.allocated(poten)) THEN
           print*, 'Using potentials to calculate gravitational forces'
           call calc_grav_from_pot
           
        ELSE IF(grav_calc_choice=='b') THEN
           print*, 'Carrying out Brute Force gravity calculation'
           call calc_grav_brute
        ELSE
           print*, 'WARNING: Parameters not clear'
           print*, 'calculating gravity via brute force to be safe'
           print*, 'structure choice: ',use_octree_grid
           print*, 'gravity calculation choice: ', grav_calc_choice
           call calc_grav_brute
        endif

    
        !***********************************
        ! 4. Write to File
        !***********************************

        ! Write gravity data to file
        call wdump_grav(ifile)
     endif
  endif

end subroutine get_SPH_neighbours
