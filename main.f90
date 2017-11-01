program sphNG_grav_plus_neighbours

  !	Code calculates gravitational forces, potentials and
  ! Computes neighbour lists for sphNG dumps

  use sphgravdata
  use treedata
  use sphkerneldata

  implicit none

  real, parameter :: pi = 3.141592653
  real, parameter :: twopi = 2.0*pi

  real :: dgridmin

  character(4) :: fileroot

  integer :: check,skip,i,n,ipart,start,finish, counter

  character(1) :: use_octree_grid, grav_calc_choice
  character(3),dimension(1000) :: num
  character(7),dimension(1000) :: filename,gravfile
  character(6),dimension(1000) :: potfile
  character(100) :: filedir, inputfile

  ! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"  
  print*, "     sphNG GRAVITY AND NEIGHBOURS CALCULATOR         "
  print*, "     Created by D. Forgan, 1st September 2014    "
  print*, "	                                       	  "
  print*, "-----------------------------------------------"
  print*, " "
  print*, " If SPH read fails, check the following:"
  print*, "    - File endianness"
  print*, "    - Default real size"
  print*, " "
  print*, " For the gfortran compiler, inserting or "
  print*, " removing the following flags should correct"
  print*, " the problem:"
  print*, "    -fconvert=swap"
  print*, "       (swaps endianness during read-in)"
  print*, "    -fdefault-real-8"
  print*, "       (sets default real to double precision)"
  print*, " "      
  print*, " "
  print*, " input parameters in ./sphNGgravplus.params"
  print*, "-----------------------------------------------"
  print*, " "

  OPEN(10,file='sphNGgravplus.params', status='old')
  READ(10,*) filedir ! File Directory
  READ(10,'(A4)') fileroot   ! File prefix
  READ(10,*) start  ! starting filenumber
  READ(10,*) finish ! Final filenumber
  READ(10,*) use_octree_grid ! Use octree (o) or grid (g) for neighbour search?
  READ(10,*) grav_calc_choice ! Use octree (o), brute force (b) or potential (p) to calc gravity forces?
  READ(10,*) dgridmin ! If using grid, this is the minimum grid length
  CLOSE(10)

  if (start>=1000) stop "Error! sphNG doesn't have filenumbers above 1000"  
  

  if (finish>=1000) then
     print*, "Warning! sphNG doesn't have filenumbers >1000"
     print*, "Setting finish=999"
     finish =999
  endif

  print*, 'sphNGgravplus.params successfully read '

  ! Write relevant filenames:
  do i = start,finish

     if (i <= 9) then
        write(num(i),'("00",I1)') i
     elseif (i <= 99) then
        write(num(i),'("0",I2)') i
     else
        write(num(i),'(I3)') i
     endif


     write(filename(i),'(A4,A3)') fileroot, num(i)
     write(gravfile(i),'("grav",A3)') num(i)
     write(potfile(i), '("pot",A3)') num(i)         
  enddo

  ! Load Kernel Tables for later interpolation

  print*, 'Loading Kernel Table Data'
  CALL ktable

  ! Loop over files
  DO n=start, finish


     !******************************************************************
     ! 1.  Read in SPH file
     !******************************************************************

     check =0
     skip = 0

     !	Read in SPH dump

     inputfile = TRIM(filedir)//filename(n)

     print*, 'File located at ', filedir
     print*, 'Address: ',inputfile

     print*, inputfile

     call rdump(inputfile,check,skip)

     ! Skip small dumps
     if(skip/=0) THEN
	print*, 'Skipping small dump'			
	cycle
     ENDIF

     if (check /= 0) then
        !	If file missing in series, skip it
	IF(n==start) THEN
           print*, 'ERROR: First file missing! Program aborted at sph read in'
           stop
        ELSE
           print*, 'Skipping missing dump ',filename(n)
           print*, "-----------------------------------------------"
           print*, " "
           cycle
	ENDIF
     endif


     allocate(isort(npart))
     allocate(iorig(npart))		


    !*************************************
    ! 2. Find neighbours for all particles
    !*************************************

     ! Find maximum and minimum values for all coordinates
     xmax = 0.0
     ymax = 0.0
     zmax = 0.0

     counter = 0
     hmean =0.0
     do ipart=1,npart		
	isort(ipart) = ipart
        iorig(ipart) = ipart

        IF(abs(xyzmh(1,ipart))> xmax) xmax = abs(xyzmh(1,ipart))
        IF(abs(xyzmh(2,ipart))> ymax) ymax = abs(xyzmh(2,ipart))
        IF(abs(xyzmh(3,ipart))> zmax) zmax = abs(xyzmh(3,ipart))

        if(iphase(ipart)==0) then
           hmean = hmean + xyzmh(5,ipart)
           counter = counter +1
        endif
     enddo

     hmean = hmean/real(counter)

     print*,  "-----------------------------------------------"
     print*, 'Mean smoothing length is ', hmean

    ! Set up neighbour lists

     allocate(nneigh(npart))
     allocate(neighb(npart,neighmax))
     nneigh(:) = 0
     neighb(:,:) = 0

     if(use_octree_grid=='o') THEN

        print*, "-----------------------------------------------"
        print*, 'Building octree'
        CALL make_octree(filename(n))

        !	Use octree to find neighbours    
  
        print*, "-----------------------------------------------"
        print*, 'Creating Neighbour Lists, npart: ',npart		
        CALL neighbours_octree(filename(n))

        print*, 'Neighbour lists created'
        print*, "-----------------------------------------------"			

     ELSE IF (use_octree_grid=='g') THEN

        print*, "-----------------------------------------------"
        print*, 'Building regular grid'
        CALL make_grid(filename(n),hmean,dgridmin)

        !	Use grid to find neighbours    
  
        print*, "-----------------------------------------------"
        print*, 'Creating Neighbour Lists from grid, npart: ',npart		
        CALL neighbours_grid(filename(n))

        print*, 'Neighbour lists created'
        print*, "-----------------------------------------------"

     ELSE
        print*, 'Finding neighbours by brute force'
        CALL neighbours_brute(filename(n))
     ENDIF
          
     !**********************************
     ! 3. Calculate Gravitational Forces
     !***********************************
          
    
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


    call deallocate_memory
 
  ENDDO
END program  sphNG_grav_plus_neighbours
