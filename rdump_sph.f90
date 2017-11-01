      subroutine rdump(filename,rcheck,skip)

!----------------------------------------------------------------------
! Subroutine to read in dump files.
! Updated from rdump_PJC.f, to be compatible with Fortan 90/95
!
! PJC 08/05/2008
!----------------------------------------------------------------------
      use sphgravdata
        
      implicit none
      integer(kind=4) :: int1,int2,itest1,itest2,itest3,itest4
      integer(kind=4) :: check,number,i,j
      integer(kind=8) :: number8
      integer,intent(inout) :: rcheck,skip
      real :: rtest1
      character(100),intent(in) :: filename

! Standard integers for read in:
      int1 = 690706
      int2 = 780806

! Open file to read in:
      print*, '1) SPH file    : ', filename
      open(10,file=filename,status='old',iostat=check,&
           form='unformatted',action='read')

! If open statement fails, return an error to main program:
      if (check /= 0) THEN
         print*, '      - SPH input file not found'
         rcheck = 1
         return
      endif
      print*, '      - Reading in SPH file'

      if (check /= 0) then
         print*, 'Error in rdump: file open failed'
         rcheck = 1
         return
      endif
!-------------------File header-----------------------------------------
! Check file header and return errors if necesary:
      read(10) itest1,rtest1,itest2,itest3,itest4
      if (itest1 /= int1) then
         print*, 'Error in rdump: ENDIANNESS wrong?'
         rcheck = 1
         return
      elseif (itest2 /= int2) then
         print*, 'Error in rdump: default integer size wrong'
         rcheck = 1
         return
      elseif (itest4 /= int1) then
         print*, 'Error in rdump: default real size wrong'
      endif
      read(10)  fileident

!-------------------Single values---------------------------------------
! Read default integers
		      
      read (10) number
      if (number /= 6) then
         print*, 'Error in rdump: No. of default integers wrong'
         rcheck = 1
         return
      endif
      read (10) npart,n1,n2,nreassign,naccrete,nkill

! Read in integers of size int*1, int*2, int*4, int*8
      do i = 1, 4
         read (10) number
      end do

!-------------------Real values-----------------------------------------
! Read default reals
      read(10) number
      if (number < 14) then
!         print*, 'Error in rdump: No. of default reals wrong'
	print*, 'Skipping small dump ',filename
	print*, ''
	skip = 1
         rcheck = 1
         return
      endif
      read(10) gt,dtmaxdp,gamma,rhozero,RK2,escap,tkin,tgrav, &
		      tterm,anglostx,anglosty,anglostz,specang,ptmassin			

! Read real*4s
      read(10) number

! Read real*8s
      read(10) number
      if (number /= 3) then
         print*, 'Error in rdump: No. of real*8s wrong'
         !rcheck = 1
         !return
      endif
      read(10) udist,umass,utime
	  
!-------------------Array header----------------------------------------
      read(10) number
      if (number /= 2) then
         print*, 'Error in rdump: No. of array types wrong'
         rcheck = 1
         return
      endif

! Read array type 1 header
      read(10) number8, (nums(i), i=1,8)
      if (number8 /= npart) then
         print*, 'Error in rdump: npart wrong'
         rcheck = 1
         return
      endif
      
! Read array type 2 header
      read(10) nptmass, (nums(i), i=1,8)                        


	print*, 'npart, nptmass ', npart, nptmass
			 
!-------------------Type 1 array values---------------------------------
! Allocate arrays as required
      allocate(iphase(npart))
      allocate(isteps(npart))
      allocate(xyzmh(5,npart))
      allocate(vxyzu(4,npart))
      allocate(rho(npart))
      allocate(dgrav(npart))
		      
	
! Read default integers!      integer,dimension(500001) :: testph

      read(10) (isteps(i), i=1,npart)
	
! Read integer*1s
      read(10) (iphase(i), i=1,npart)

	
! Read default reals
      do j = 1,5
         read(10) (xyzmh(j,i), i=1,npart)
      enddo
      do j = 1,4
         read(10) (vxyzu(j,i), i=1,npart)
      enddo
	
! Read real*4s
      read(10) (rho(i), i=1,npart)
      read(10) (dgrav(i), i=1,npart) 

	print*, 'SPH particles read'
			 
!-------------------Type 2 array values---------------------------------
! Allocate arrays as required

	allocate(listpm(nptmass))
	allocate(spinx(nptmass))
	allocate(spiny(nptmass))
	allocate(spinz(nptmass))
	allocate(angaddx(nptmass))
	allocate(angaddy(nptmass))
	allocate(angaddz(nptmass))
	allocate(spinadx(nptmass))
	allocate(spinady(nptmass))
	allocate(spinadz(nptmass))
			 

			 
! Read default integers
      read (10) (listpm(i), i=1,nptmass)

! Read default reals
      read (10) (spinx(i),i=1,nptmass)
      read (10) (spiny(i),i=1,nptmass)
      read (10) (spinz(i),i=1,nptmass)
      read (10) (angaddx(i),i=1,nptmass)
      read (10) (angaddy(i),i=1,nptmass)
      read (10) (angaddz(i),i=1,nptmass)
      read (10) (spinadx(i),i=1,nptmass)
      read (10) (spinady(i),i=1,nptmass)
      read (10) (spinadz(i),i=1,nptmass)

      print*, 'pointmasses read'
      close(10)

      print*, '      - SPH dump file correctly read in'
      print*, '      -',npart,'particles in total'
      print*, '      -',npart-naccrete-1,'gas particles active'
      print*, ' '

      return
     
    
      end subroutine rdump
