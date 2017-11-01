subroutine get_SPH_neighbours
!
! Subroutine either reads SPH neighbour data from file, or
! computes it via a regular grid
!


! First, test to see if neighbour file exists: if it does, no need to
! read or calculate octrees or calculate neighbours

neighbourfile ="neighbours_"//TRIM(filename(n))

inquire(file=neighbourfile,exist = existneigh)

if(existneigh.eqv..true.) then

   print*, 'Neighbour file',neighbourfile, ' found: reading'
   call read_neighbours(neighbourfile)

else
   ! If no neighbour file found, then we must generate the list
   
   allocate(nneigh(nelement))
   allocate(neighb(nelement,neighmax))
   nneigh(:) = 0
   neighb(:,:) = 0

   if(use_octree_grid=='o') THEN

      print*, "-----------------------------------------------"
      print*, 'Building octree'
      CALL make_octree(filename(n))
      
      !	Use octree to find neighbours    
  
      print*, "-----------------------------------------------"
      print*, 'Creating Neighbour Lists, nelement: ',nelement		
      CALL neighbours_octree(filename(n))
      
      print*, 'Neighbour lists created'
      print*, "-----------------------------------------------"			
      
     ELSE IF (use_octree_grid=='g') THEN

        print*, "-----------------------------------------------"
        print*, 'Building regular grid'
        CALL make_grid(filename(n),hmean,dgridmin)

        !	Use grid to find neighbours    
  
        print*, "-----------------------------------------------"
        print*, 'Creating Neighbour Lists from grid, nelement: ',nelement		
        CALL neighbours_grid(filename(n))

        print*, 'Neighbour lists created'
        print*, "-----------------------------------------------"

     ELSE
        print*, 'Finding neighbours by brute force'
        CALL neighbours_brute(filename(n))
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
        call wdump_grav(gravfile(n),potfile(n))
     endif
  endif

end subroutine get_SPH_neighbours
