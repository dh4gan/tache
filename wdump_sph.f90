      subroutine wdump(filename,rcheck)

!----------------------------------------------------------------------
! Subroutine to write dump files.
! Updated from wdump.f, to be compatible with Fortan 90/95
!
! DHF 08/10/2009
!----------------------------------------------------------------------
      use sphdata
        
      implicit none
      integer(kind=4) :: int1,int2,i1
      integer(kind=4) :: check,number,i,j
      integer(kind=8) :: number8
      integer,intent(inout) :: rcheck
      real :: r1
      character(7),intent(in) :: filename

! Standard integers for write:
      int1 = 690706
      int2 = 780806
	  r1 = real(int2)
	  i1 = int1

		    
! Open file to write to:
      print*, 'Writing SPH file    : ', filename
      open(10,file=filename,status='unknown',iostat=check,&
           form='unformatted',action='write')

! If open statement fails, return an error to main program:
 
      if (check /= 0) then
         print*, 'Error in wdump: file open failed'
         rcheck = 1
         return
      endif
!-------------------File header-----------------------------------------
! Check file header and return errors if necesary:
      write(10) int1,r1,int2,i1,int1
	  fileident= 'FHydro_hybrid'
      write(10)  fileident

!-------------------Single values---------------------------------------
! Write default integers
	  number = 6	      
      write(10) number
	  write(10) npart,n1,n2,nreassign,naccrete,nkill

! Write in integers of size int*1, int*2, int*4, int*8
	  number = 0
      do i = 1, 4
         write(10) number
      end do

!-------------------Real values-----------------------------------------
! Write default reals
	  number = 14
      write(10) number
      write(10) gt,dtmaxdp,gamma,rhozero,RK2,escap,tkin,tgrav, &
		      tterm,anglostx,anglosty,anglostz,specang,ptmassin			

! Write real*4s
	  number=0
      write(10) number

! Write real*8s
	  number = 3
      write(10) number
      write(10) udist,umass,utime
      
!-------------------Array header----------------------------------------
	  number = 2
	  write(10) number
	   
! Write array type 1 header
	  number8 = npart
      nums(1) = 1
      nums(2) = 1
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 9
      nums(7) = 2
      nums(8) = 0

      write(10) number8, (nums(i), i=1,8)
            
! Write array type 2 header
	  number8 = nptmass
      nums(1) = 1
      nums(2) = 0
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 9
      nums(7) = 0
      nums(8) = 0

      write(10) number8, (nums(i), i=1,8)                        

!-------------------Type 1 array values---------------------------------
		     	
! Write default integers!      integer,dimension(500001) :: testph

      write(10) (isteps(i), i=1,npart)
	
! Write integer*1s
      write(10) (iphase(i), i=1,npart)

	
! Write default reals
      do j = 1,5
         write(10) (xyzmh(j,i), i=1,npart)
      enddo
      do j = 1,4
         write(10) (vxyzu(j,i), i=1,npart)
      enddo
	
! Write real*4s
      write(10) (rho(i), i=1,npart)
      write(10) (dgrav(i), i=1,npart)
			 	
!-------------------Type 2 array values---------------------------------


			 
! Write default integers
      write(10) (listpm(i), i=1,nptmass)

! Write default reals
      write (10) (spinx(i),i=1,nptmass)
      write (10) (spiny(i),i=1,nptmass)
      write (10) (spinz(i),i=1,nptmass)
      write (10) (angaddx(i),i=1,nptmass)
      write (10) (angaddy(i),i=1,nptmass)
      write (10) (angaddz(i),i=1,nptmass)
      write (10) (spinadx(i),i=1,nptmass)
      write (10) (spinady(i),i=1,nptmass)
      write (10) (spinadz(i),i=1,nptmass)

      close(10)

      print*, '      - SPH dump file written'
      print*, '      -',npart,'particles in total'
      print*, '      -',npart-naccrete-nptmass,'gas particles active'
      print*, ' '
	
	print*, '      - complete'
      return

      end subroutine wdump
