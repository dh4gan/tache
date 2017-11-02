PROGRAM link_clumps
  ! Program takes outputs from clumpfind
  ! Uses clump membership to identify the same clump over several timesteps, 
  ! and labels them correctly

 ! use clumpdata

  implicit none

 ! type(sphclump), allocatable :: clump(:) 

  integer :: start,finish,i,j,k,l,nelement, kprev, newclump
  integer,allocatable,dimension(:) :: nclumps

  integer, dimension(1) :: temp
  real,allocatable, dimension(:)  :: clumptot,newtot
  integer,allocatable, dimension(:,:) :: member,matches

  character(3),dimension(1000) :: num
  character(100), dimension(1000) :: cataloguefile
  character(100), dimension(1000) :: memberfile	
  character(7) :: file
  character(4) :: fileroot


  logical :: fileexist

  !	Read in user input			
  print*, " "
  print*, "-----------------------------------------------"  
  print*, "           CLUMP LINKING"
  print*, "     Created by D. Forgan, 7th July 2011    "
  print*, "	                            	  "
  print*, "-----------------------------------------------"
  print*, " "

  print*, "What is the initial input filename?"
  read*, file
  print*, 'How many particles does this dump have?'
  read*, nelement

  do 
     print*, "What is the final filenumber?"
     read*, finish
     if (finish < 1000) exit 
     print*, "Final filenumber too great - retry!"
  enddo

  write(fileroot,'(A4)') file(1:4)
  read(file(5:7),'(I3)') start

  ! Create filename strings to read in

  do i = start, finish
     if (i <= 9) then
        write(num(i),'("00",I1)') i                                 !(can use write to write into an array?)
     elseif (i <= 99) then
        write(num(i),'("0",I2)') i
     else
        write(num(i),'(I3)') i
     endif
     write(cataloguefile(i),'("catalogue_",A4,A3)')fileroot,num(i)
     cataloguefile(i) = TRIM(cataloguefile(i))                            ! (trim gets rid of any trailing blanks)
     write(memberfile(i),'("members_",A4,A3)') fileroot,num(i)
     !	 write(clumpsphfile(i), '("",A3)') num(i)
  enddo

  ! Open file to store linkages
  OPEN(20,file='clumplinks.dat', status='unknown')

  print*, " "
  print*, "-----------------------------------------------"

 ! allocate(clump(finish))
  allocate(nclumps(finish))                                                !(finish is just the number of array dimensions not the number of clumps)
  allocate(member(finish,nelement))

  !***************************************
  ! 1. Identify original clump memberships
  !***************************************

  OPEN(10,file=cataloguefile(start), status='old')                                       !(opening the first file to use as a reference and define each clump number)
  read(10,*) nclumps(start)
  CLOSE(10)

  OPEN(10,file=memberfile(start), status='old',form='unformatted')	
  read(10) (member(start,i), i=1,nelement)	
  CLOSE(10)

  !	Identify initial clump totals

  allocate(clumptot(nclumps(start)))
  DO i=1,nelement
     IF(member(start,i)==0) cycle
     clumptot(member(start,i)) = clumptot(member(start,i)) + 1
  ENDDO

  print*, 'Initial clump totals'
  print*, "-----------------------------------------------"

  DO j=1,nclumps(start)
     print*, j, clumptot(j)
  ENDDO

  print*, "-----------------------------------------------"

  kprev = start

  !**************************************************************************
  ! 2. Begin loop over all subsequent dumps, and compare against initial dump
  !**************************************************************************

  DO k=start+1,finish

     print*, 'Checking for files at dump no. ',k
     INQUIRE(file=cataloguefile(k), EXIST=fileexist)

     !	If file doesn't exist skip			
        IF(fileexist.eqv..FALSE.) cycle
        
   
     

     !	Read catalogue file - get total number of clumps
     OPEN(10,file=cataloguefile(k),status='old')
     READ(10,*) nclumps(k)			
     CLOSE(10)

     ! Read membership file
     OPEN(10,file=memberfile(k),status='old',form='unformatted')
     READ(10) (member(k,i), i=1,nelement)
     CLOSE(10)

     ! If there are less clumps here than in the initial file, print this to screen
     IF(nclumps(start)>nclumps(k)) THEN
	print*, 'Clumps missing ', kprev, k, nclumps(start), nclumps(k)
     ENDIF

     ! Keep track of new totals for each clump
     allocate(newtot(nclumps(k)))
     newtot(:) = 0.0

     DO i=1,nelement
	IF(member(k,i)==0) cycle		
	newtot(member(k,i)) = newtot(member(k,i)) + 1
     ENDDO
     print*, "-----------------------------------------------"		
     print*, 'New clump totals for dump ',k
     print*, "-----------------------------------------------"		
     DO j=1,nclumps(k)
	print*, j, newtot(j)
     ENDDO
     print*, "-----------------------------------------------"

     ! Compare previous dump's membership with new dump
     ! Create matrix of matches
     ! add 1 to matrix element (i,j) if particle is in clump i of dump kprev and clump j of dump k

     allocate(matches(nclumps(kprev),nclumps(k)))
     matches(:,:) = 0.0	

     !************************************************************************
     ! 3. Record all instances of particle matches in clumps from kprev and k
     !************************************************************************
     DO i=1,nelement
        IF(member(kprev,i)*member(k,i)/=0) THEN
           matches(member(kprev,i),member(k,i)) = matches(member(kprev,i),member(k,i)) + 1
        ENDIF
     ENDDO

     print*, 'Clump populations ready for comparison'

     DO l=1,nclumps(k)

        IF(l <=nclumps(kprev)) THEN 

           ! Find clump in new dump that contains the most matches to clump in kprev
           temp = maxloc(matches(:,l))
           newclump = temp(1)

           !	Remove its contributions to other clump counts (prevents multiple identifications)
           matches(newclump,l+1:nclumps(k)) = 0.0

           print*, k, l, newclump, maxval(matches(:,l))/REAL(clumptot(newclump))


           WRITE(20,*) num(k), l, newclump, maxval(matches(:,l))/REAL(clumptot(newclump))

        ELSE IF(l>nclumps(kprev)) THEN
           WRITE(20,*) num(k), l, l, 0.0
	ENDIF

        CALL FLUSH(20) 		 		

        !	End of loop over clumps

     ENDDO


     ! Pass current clump data to old clump variables in preparation for next loop
     kprev = k

     deallocate(clumptot)
     allocate(clumptot(nclumps(k)))

     clumptot(:) = newtot(:)

     deallocate(newtot,matches)

     print*, 'Run complete for dump ',num(k)

     !	End of loop over dumps				

  ENDDO

  close(20)

END PROGRAM link_clumps
