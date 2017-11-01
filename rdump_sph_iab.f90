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
      integer(kind=4) :: check,number,i,j, nptcnt, ic1,ic2,ic3, jj
      integer(kind=8) :: number8, nlistinactive, numberdummy
      integer, dimension(8) :: numssink,numsRT
      integer,allocatable,dimension(:) :: listinactive
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
         print*, 'Error in rdump: ENDIANNESS wrong? ',itest1,int1
         !rcheck = 1
         !return
      elseif (itest2 /= int2) then
         print*, 'Error in rdump: default integer size wrong'
         !rcheck = 1
         !return
      elseif (itest4 /= int1) then
         print*, 'Error in rdump: default real size wrong'
      endif
      read(10)  fileident

      print*, 'Ident: ',fileident
!-------------------Single values---------------------------------------
! Read default integers
		      
      read (10) number

      if(number==6) then
         read(10) npart,n1,n2,nreassign,naccrete,nkill
         nblocks = 1
      else if(number>=7) then
         read(10) npart,n1,n2,nreassign,naccrete,nkill,nblocks
      else
         print*, 'Error in rdump: No. of default integers wrong: ',number
         rcheck = 1
         return
      endif

      print*, 'npart, nblocks: ', npart, nblocks
! Read in integers of size int*1, int*2, int*4, int*8
      do i = 1, 6
         read (10) number
         print*, number
      end do

!-------------------Real values-----------------------------------------
! Read default reals

      read(10) gt,dtmaxdp,gamma,rhozero,RK2,escap,tkin,tgrav, &
		      tterm,anglostx,anglosty,anglostz,specang,ptmassin			

! Read real*4s
      read(10) number

! Read real*8s
      read(10) number
      if (number /= 3) then
         !print*, 'Error in rdump: No. of real*8s wrong'
         !rcheck = 1
         !return
      endif
      read(10) udist,umass,utime
	  
      write(*,*) 'Units - Distance, Mass, Time: ',udist, umass, utime
!-------------------Array header----------------------------------------
      read(10) number
      print*, 'Number of Array Types: ', number, number/nblocks
      number = number/nblocks
      if (number /= 2) then
         print*, 'Error in rdump: No. of array types wrong'
         rcheck = 1
         return
      endif


 icount = 0
 nptcnt = 0 

print*, 'Allocating Arrays'
allocate(xyzmh(5,npart))
allocate(vxyzu(4,npart))
allocate(iphase(npart))
allocate(isteps(npart))
allocate(rho(npart))
allocate(dgrav(npart))
allocate(listinactive(npart))

allocate(listpm(nptmax))
allocate(spinx(nptmax))
allocate(spiny(nptmax))
allocate(spinz(nptmax))
allocate(angaddx(nptmax))
allocate(angaddy(nptmax))
allocate(angaddz(nptmax))
allocate(spinadx(nptmax))
allocate(spinady(nptmax))
allocate(spinadz(nptmax))

print*, 'Arrays Allocated'

 DO jj = 1, nblocks

    READ (10, END=20) number8, (nums(i), i=1,8)
    print *, npart, nblocks,jj,number8,number8*nblocks
    npart = number8

    print *,'Reading block ',jj,npart, icount
    print *,'nums ',(nums(i), i=1,8)

    READ (10, END=20) number8, (numssink(i), i=1,8)
    nptmass = number8
    IF (number.EQ.3) THEN
       READ (10, END=20) number8, (numsRT(i), i=1,8)
    ENDIF

    print *,'numssink ',(numssink(i), i=1,8)
    print *,'numsRT ',(numsRT(i), i=1,8)

    ic1 = icount+1
    ic2 = icount+2
    ic3 = icount+3
    READ (10, END=20) (isteps(i), i=icount+1, & 
         icount+npart)
     print *,'isteps: ',isteps(ic1),isteps(ic2),isteps(ic3)
     print*, nums(1)
     IF (nums(1).GE.2) THEN
        READ(10) 
 !       READ (11, END=20) (nlistinactive, listinactive(i), &
  !               i=icount+1, icount+nlistinactive)
        print *, 'inactive: ', nlistinactive,listinactive(ic1)
     END IF
     READ (10, END=20) (iphase(i), i=icount+1, &
          icount+npart)


     !
     !--Read iunique
     !

     READ (10, END=20) numberdummy
     READ (10, END=20) (xyzmh(1,i), i=icount+1, icount+npart)
     READ (10, END=20) (xyzmh(2,i), i=icount+1, icount+npart)
     READ (10, END=20) (xyzmh(3,i), i=icount+1, icount+npart)
     READ (10, END=20) (xyzmh(4,i), i=icount+1, icount+npart)
     READ (10, END=20) (xyzmh(5,i), i=icount+1, icount+npart)

     print *,'reading velocities and u'

     READ (10, END=20) (vxyzu(1,i), i=icount+1, icount+npart)
     READ (10, END=20) (vxyzu(2,i), i=icount+1, icount+npart)
     READ (10, END=20) (vxyzu(3,i), i=icount+1, icount+npart)
     READ (10, END=20) (vxyzu(4,i), i=icount+1, icount+npart)

     print *,'reading density'

     READ (10, END=20) (rho(i), i=icount+1, &
          icount+npart)

     write (*,*) 'nums(7) is ',nums(7)

     DO j = 1, nums(7)-2
        READ (10, END=20) 
     END DO
     READ (10, END=20) (dgrav(i), i=icount+1, &
                 icount+npart)

     print*, 'Reading Sink Data'

     READ (10, END=20) (listpm(i),i=nptcnt+1,nptcnt+nptmass)

     DO i=nptcnt+1,nptcnt+nptmass
        listpm(i) = listpm(i)+icount
     END DO

     READ (10, END=20) (spinx(i),i=1,nptmass)
     READ (10, END=20) (spiny(i),i=1,nptmass)
     READ (10, END=20) (spinz(i),i=1,nptmass)
     READ (10, END=20) (angaddx(i),i=1,nptmass)
     READ (10, END=20) (angaddy(i),i=1,nptmass)
     READ (10, END=20) (angaddz(i),i=1,nptmass)
     READ (10, END=20) (spinadx(i),i=1,nptmass)
     READ (10, END=20) (spinady(i),i=1,nptmass)
     READ (10, END=20) (spinadz(i),i=1,nptmass)


     print *,'numssink(6)=',numssink(6)
     DO i = 1, numssink(6)-9
        READ (10, END=20)
     END DO
     DO i = 1, numssink(8)
        READ (10, END=20)
     END DO

!
!--rad trans
!
     IF (number.EQ.3) THEN
       print *,'Reading RT ',number,(numsRT(i),i=1,8)
!!$        IF (numsRT(3).EQ.1) THEN
!!$           READ (11, END=20) (nneigh(i), i=1, npart)
!!$        ENDIF
!!$        READ (11, END=20) (e(i), i=1, npart)
!!$        READ (11, END=20) (rkappa(i), i=1, npart)
!!$        READ (11, END=20) (cv(i), i=1, npart)
!!$        READ (11, END=20) (rlambda(i), i=1, npart)
!!$        READ (11, END=20) (edd(i), i=1, npart)
!!$        IF (numsRT(6).EQ.8) THEN
!!$           DO ii=1,3
!!$              READ (11, END=20) (force(ii,i), i=1, npart)
!!$           END DO
!!$        ENDIF
!!$        IF (numsRT(7).EQ.1) THEN
!!$           READ (11, END=20) (dlnTdlnP(i), i=1, npart)
!!$        ENDIF
!!$        IF (numsRT(7).EQ.13 .OR. numsRT(7).EQ.14) THEN
!!$           READ (11, END=20) (dlnTdlnP(i), i=1, npart)
!!$           IF (numsRT(7).EQ.11) THEN
!!$              READ (11, END=20) (adiabaticgradient(i), &
!!$                          i=1, npart)
!!$           ENDIF
!!$           DO ii=1,3
!!$              READ (11, END=20) (pressure(ii,i),  &
!!$                          i=1, npart)
!!$           END DO
!!$           DO ii=1,3
!!$              READ (11, END=20) (viscosity(ii,i), &
!!$                          i=1, npart)
!!$           END DO
!!$           DO ii=1,3
!!$              READ (11, END=20) (gravity(ii,i), 
!!$              &                          i=1, npart)
!!$           END DO
!!$           DO ii=1,3
!!$              READ (11, END=20) (radpres(ii,i), 
!!$              &                          i=1, npart)
!!$           END DO
!!$        ENDIF
     ENDIF
     
     icount = icount + npart
     nptcnt = nptcnt + nptmass
     !
     !--End block loop
     !
  END DO
  npart = icount
  nptmass = nptcnt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


close(10)

print*, 'npart, nptmass ', npart, nptmass
			 

      print*, '      - SPH dump file correctly read in'
      print*, '      -',npart,'particles in total'
      print*, '      -',npart-naccrete-1,'gas particles active'
      print*, ' '

20    CONTINUE
      return
     
    
      end subroutine rdump
