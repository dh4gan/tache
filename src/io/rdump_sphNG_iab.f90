! Written 26/2/15 by dh4gan
! Code adapted from rbin_v2 module written by William Lucas (St Andrews)
! Reads modern sphNG binary format

! rdump_sph is the principal subroutine, but it relies on several other 
! subroutines to read the header and the data in contiguous blocks


!IMPORTANT:
!  If compiled with ifort there may be issues with reading large files.
!  If on linux, unlimiting the stacksize first should help; on OSX, where this
!  isn't possible, compile with -heap-arrays 6400 to force creation of temporary
!  arrays on the heap instead. (This is just a number I've found to work :)
!  gfortran is by default more open to using the heap, so this isn't a problem.

!  If working on a multiple-block (MPI) file, choose 'contiguous' before use.
!  With .TRUE. the whole file is read and all data stored in one array.
!  With .FALSE. a single block is read and control passed back to the
!  controlling subroutine - any data already in the arrays is destroyed.

!  Ensure module is scoped correctly or data may be lost - if you don't
!  feel like being careful, use global scope ('USE rbin' in the main program).

SUBROUTINE rdump_sphNG_iab(filename,rcheck, skip)

use sphdata
use tachedata,only: nelement

implicit none 
integer :: i,check, rcheck,skip
character(100) :: filename

! No small dumps to read
skip =0
! Open file to read in:
      print*, '1) SPH file    : ', TRIM(filename)
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

      ! Allocate sink arrays

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

      call read_sphng_header
      call read_sphng_data

      if(nblocks>1) nelement = nelementtot
      close(10)

      print*, 'NELEMENT ', nelement
      allocate(isort(nelement))
      allocate(iorig(nelement))
      do i=1,nelement
         isort(i) = i
         iorig(i) = i
      enddo
return 
END SUBROUTINE rdump_sphNG_iab


! Subroutine reads the data file
! This contains flow control to either read the entire file as a single block
! or read multiple blocks in a contiguous fashion (MPI dumps)

  SUBROUTINE read_sphng_data

    use sphdata
    IMPLICIT NONE
    integer :: ii
    !If only 1 block exists or non-contiguous mode, only read once.
    IF (.NOT. contiguous .OR. nblocks == 1) THEN
       CALL read_sphng_block
       !If we want an entire MPI run, call multiple times.
    ELSE IF (contiguous) THEN
       DO ii = 1, nblocks
          CALL read_sphng_block
       END DO
    ELSE
       print *, "Unsure how to read data:"
       print *, "nblocks = ", nblocks
       IF (contiguous) THEN
          print *, "contiguous mode: ON"
       ELSE
          print *, "contiguous mode: OFF"
       END IF
       print *, "Exiting..."
       STOP
    END IF

    print*, 'sphNG file read complete'
    RETURN
  END SUBROUTINE read_sphng_data
  

  ! This subroutine reads the header data at the top of the binary file
  SUBROUTINE read_sphng_header

    use sphdata
    IMPLICIT NONE
  
    integer :: i !dummy variable
    integer :: ierr !error return from extract
 
    integer :: number
    
    character*16, dimension(128) :: tagsreal
    
    integer, parameter :: idimhead = 30
    real, dimension(idimhead) :: rheader
    
    print *, "Starting file read..."
  
    read(10)
    read(10) fileident
    print *, trim(fileident)
    IF (fileident(1:1).NE.'F') THEN
       print *, "File does not contain a full dump, exiting!"
       CLOSE(10)
       STOP
    END IF
    IF (fileident(2:2).EQ.'T') THEN
       tagged = .TRUE.
       print *, "File is tagged."
    ENDIF
    
    read(10) number
    IF (tagged) read(10) !skip tags here
    IF (number == 6) THEN
       read(10) nelementtot,n1,n2,nreassign,naccrete,nkill
       nblocks = 1
    ELSE
       read(10) nelementtot,n1,n2,nreassign,naccrete,nkill,nblocks
    END IF
    write (*,*) 'nelementtot = ',nelementtot,' nblocks ',nblocks
    allocate( nelementblocks(nblocks) )
    
    DO i = 1, 3
       read(10) 
    END DO
    read(10) number
    IF (number == 1) THEN
       IF (tagged) read(10)
       read(10) iuniquemax
    ELSE
       iuniquemax = nelementtot
    END IF
    print *, "iuniquemax =", iuniquemax
    
    read(10) number
    print *, "number real", number
    
    IF (tagged) THEN
       write(*,"(A,I3)") "Reading tags from file, number =", number
       read(10) tagsreal(1:number)
    ELSE
       print *, "Simulating tags for non-tagged file."
       tagsreal(1:30) = &
            (/'gt              ','dtmax           ','gamma           ', &
            'rhozero         ','RK2             ','escap           ', &
            'tkin            ','tgrav           ','tterm           ', &
            'anglostx        ','anglosty        ','anglostz        ', &
            'specang         ','ptmassin        ','tmag            ', &
            'Bextx           ','Bexty           ','Bextz           ', &
            'hzero           ','uzero_n2        ','hmass           ', &
            'gapfac          ','                ','sdprof          ', &
            'rorbit_orig     ','min_rplan       ','max_rplan       ', &
            'planetesimalmass','coremass_orig   ','coremass        '/)
    END IF
    
    read(10) (rheader(i),i=1,min(number,idimhead))
    
    CALL extract('gt',gt,rheader,tagsreal,number,ierr)
    IF (ierr.NE.0) STOP
    print *, "Time is", gt
    CALL extract('dtmax',dtmaxdp,rheader,tagsreal,number,ierr)
    IF (ierr.NE.0) STOP
    CALL extract('gamma',gamma,rheader,tagsreal,number,ierr)
    IF (ierr.NE.0) STOP
    CALL extract('rhozero',rhozero,rheader,tagsreal,number,ierr)
    IF (ierr.NE.0) STOP
    CALL extract('RK2',RK2,rheader,tagsreal,number,ierr)
    IF (ierr.NE.0) STOP
    CALL extract('escap',escap,rheader,tagsreal,number,ierr)
    CALL extract('tkin',tkin,rheader,tagsreal,number,ierr)
    CALL extract('tgrav',tgrav,rheader,tagsreal,number,ierr)
    CALL extract('tterm',tterm,rheader,tagsreal,number,ierr)
    CALL extract('anglostx',anglostx,rheader,tagsreal,number,ierr)
    CALL extract('anglosty',anglosty,rheader,tagsreal,number,ierr)
    CALL extract('anglostz',anglostz,rheader,tagsreal,number,ierr)
    CALL extract('specang',specang,rheader,tagsreal,number,ierr)
    CALL extract('ptmassin',ptmassin,rheader,tagsreal,number,ierr)
    IF (imhd == 1) THEN
       CALL extract('tmag',tmag,rheader,tagsreal,number,ierr)
       CALL extract('Bextx',Bextx,rheader,tagsreal,number,ierr)
       CALL extract('Bexty',Bexty,rheader,tagsreal,number,ierr)
       CALL extract('Bextz',Bextz,rheader,tagsreal,number,ierr)
       print *, 'External field found, Bext = ',Bextx,Bexty,Bextz
    END IF
    CALL extract('hzero',hzero,rheader,tagsreal,number,ierr)
    IF (iexf == 10 .AND. ierr /= 0) STOP      
    CALL extract('uzero_n2',uzero_n2,rheader,tagsreal,number,ierr)
    IF (n2 > 0) THEN
       print *, ' read u for surrounding medium = ',uzero_n2
    ENDIF
    CALL extract('hmass',hmass,rheader,tagsreal,number,ierr)
    CALL extract('gapfac',gapfac,rheader,tagsreal,number,ierr)
    CALL extract('sdprof',sdprof,rheader,tagsreal,number,ierr)
    IF (ierr /= 0) THEN
       sdprof = -0.5
       print *, 'Surface density pre vary (goes as r^-0.5)'
    ENDIF
    CALL extract('rorbit_orig',rorbit_orig,rheader,tagsreal,number,ierr)
    IF (imigrate > 0) rorbitmax = (rorbitmax - rorbit_orig)/pmrate
    CALL extract('min_rplan',min_rplan,rheader,tagsreal,number,ierr)
    CALL extract('max_rplan',max_rplan,rheader,tagsreal,number,ierr)
    IF (ierr /= 0) max_rplan = 1.0E6
    print *, 'Setting planetesimal boundaries; check values.'
    print *, 'r_min = ', min_rplan
    print *, 'r_max = ', max_rplan
    IF (.NOT.tagged) THEN
       print *, 'DO NOT USE CODE WITH DUMPS MADE FROM APRIL 2011'
       print *, ' AND BEFORE 16/05/2011. For these, reals 26, 27 are'
       print *, ' planetesimal radius and density respectively.'
    ENDIF
    CALL extract('planetesimalmass',planetesimalmass,rheader,tagsreal,number,ierr)
    print *, 'Planetesimal mass ', planetesimalmass
    CALL extract('coremass_orig',coremass_orig,rheader,tagsreal,number,ierr)
    print *, 'Core mass orig ', coremass_orig
    CALL extract('coremass',coremass,rheader,tagsreal,number,ierr)
    print *, 'Core mass running record ', coremass
    
  read(10) number
  read(10) number
  IF (number < 3) THEN
     print *, "Error in rbin, nreal8 too small in header section"
     STOP
  END IF
  IF (tagged) read(10) !skip tags
  IF (imhd == 1) THEN
     IF (number > 3) THEN
        read(10) udist, umass, utime, umagfd
     ELSE
        print *, 'WARNING: no mag field units in rdump'
        read(10) udist, umass, utime
        !umagfdi = umagfd
        !Will need to define B field units another way      
     END IF
  ELSE
     read(10) udist, umass, utime
  END IF

  !Old code to read general dump info
  !  read(10) time, dtmaxdp, gamma, rhozero, RK2, escap, tkin, tgrav, tterm, &
  !    anglostx, anglosty, anglostz, specang, ptmassin
  !  print *, "Time is", time
  !  read(10)
  !  read(10)
  
  !  read(10) udist, umass, utime
  !  print *, "Distance, mass, time units are:",udist,umass,utime,"in cgs"
  
  read(10) number
  numberarray = number/nblocks
  print *, "Array types", number, numberarray
  !IF (numberarray /= 2) THEN
  !  print *, "Expected 2 for no MHD & no RT - exiting"
  !  CLOSE(10)
  !  STOP
  !END IF
  
  icount = 0
  icountsink = 0
  
  !At this point we are onto the first data block
  iblock = 1
  
  RETURN
  END SUBROUTINE read_sphng_header
    
  !Reads in a single block of data.
  SUBROUTINE read_sphng_block
    use sphdata
    use tachedata, only: nelement

  IMPLICIT NONE
  
  character(len=16) :: tagi
  
  integer :: i !dummy integer
  integer*8 :: number8
  integer, dimension(8) :: numssink, numsRT
  integer :: ic1, ic2, ic3
  
  print *, " "
  print *, " "
  print *, "   Reading block ", iblock
  
  read(10) number8, nums(1:8)
  !print *, nelementtot, nblocks, iblock, number8, number8*nblocks
  nelement = number8
  nelementblocks(iblock) = nelement
  
  print *, "   Block contains ", nelement, "particles, current icount ", icount
  !print *, "nums ", nums(1:8)
  
  read(10) number8, numssink(1:8)
  nptmass = number8
  !print *, "numssink ", numssink(1:8)
  IF (numberarray == 3) THEN
    read(10) number8, numsRT(1:8)
   ! print *, "numsRT ", numsRT(1:8)
  END IF
  
  ic1 = icount + 1
  ic2 = icount + 2
  ic3 = icount + 3
  
  !Re/allocate as needed - isteps is taken to be representative of all arrays
  IF (.NOT. allocated(isteps)) THEN
    CALL allocate_arrays
  !Need to reallocate if we're on blockwise read or if we're starting to read a new file
  ELSE IF (.NOT. contiguous .OR. iblock == 1) THEN
    CALL reallocate_arrays
  END IF
  
  IF (tagged) read(10) !skip tag
  read(10) isteps(icount+1:icount+nelement)
  !print *, "isteps: ", isteps(ic1:ic3)
  IF (nums(1) >= 2) THEN
    IF (tagged) read(10) !skip tag
    read(10) !skip reading listinactive
  END IF
  IF (tagged) THEN
    read(10) tagi
    IF (trim(tagi) /= "iphase") THEN
      print*," WARNING: iphase does not match tag "//trim(tagi)
    END IF
  END IF
  read(10) iphase(icount+1:icount+nelement)
  !print *, "iphase: ", iphase(ic1:ic3)
  
  !In rdump, only if nums1(5) >= 5
  IF (tagged) read(10)
  read(10) iunique(icount+1:icount+nelement)
  
  !print *, "Reading xyzmh..."
  DO i = 1, 5
    IF (tagged) read(10) tagi !skip tag
    read(10) xyzmh(i,icount+1:icount+nelement)
  END DO
  !print *, "x: ", xyzmh(1,ic1:ic3)
  !print *, "h: ", xyzmh(5,ic1:ic3)
  
  !print *, "Reading vxyzu..."
  DO i = 1, 4
    IF (tagged) read(10) tagi !skip tag
    read(10) vxyzu(i,icount+1:icount+nelement)
  END DO
  !print *, "vx: ", vxyzu(1,ic1:ic3)
  !print *, "u: ", vxyzu(4,ic1:ic3)
  
  !print *, "Reading rho..."
  IF (tagged) read(10) tagi !skip tag
  read(10) alphaMM(icount+1:icount+nelement)
  rho(icount+1:icount+nelement) = alphaMM(icount+1:icount+nelement)
  !print *, "rho: ", rho(ic1:ic3)
  
  !print *, "nums(7) is ", nums(7)
  IF (igradh == 1) THEN
    IF (nums(7) >= 2) THEN
   !   print *, "Reading gradh..."
      IF (tagged) read(10) tagi !skip tag
    !  print *, "tag is ", trim(tagi)  
      read(10) gradh(icount+1:icount+nelement)
    !  print *, "gradh: ", gradh(ic1:ic3)
    END IF
    IF (nums(7) >= 3) THEN
     ! print *, "Reading gradhsoft..."
      IF (tagged) read(10) tagi !skip tag
      !print *, "tag is ", trim(tagi)  
      read(10) gradhsoft(icount+1:icount+nelement)
      !print *, "gradhsoft: ", poten(ic1:ic3)
    END IF
  ELSE
    !print *, "Reading dgrav..."
    IF (tagged) read(10) tagi !skip tag
    !print *, "tag is ", trim(tagi)
    read(10) dgrav(icount+1:icount+nelement)
    !print *, "dgrav: ", dgrav(ic1:ic3)  
  END IF 
  
  
!  print *, "nums(7) is ", nums(7)
!  IF (tagged) read(10) tagi !skip tag
!  DO i = 1, nums(7)-2
!    read(10)
!  END DO
 
  !print *, "Reading alphaMM..."
  IF (tagged) read(10) tagi !skip tag
  read(10) alphaMM(icount+1:icount+nelement)
  !print *, "alphaMM: ", alphaMM(ic1:ic3)

  ! Read the potential if it is here
  if(nums(7).ge.5) then
     print*, 'reading potential'
     if(tagged) read(10) tagi
     read(10) poten(icount+1:icount+nelement)
     print*, 'poten: ', poten(ic1:ic3)
  endif
  
  !If nptmass is zero, this may be incorrect (though it runs anyway).
  !It may be safer to read to a buffer and only store if nptmass > 0.
  print *, "   Reading sink particle data, nptmass: ", nptmass
  IF (tagged) read(10) tagi !skip tag
  read(10) listpm(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) spinx(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) spiny(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) spinz(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) angaddx(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) angaddy(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) angaddz(icountsink+1:icountsink+nptmass)
  IF (tagged) read(10) tagi !skip tag
  read(10) spinadx(icountsink+1:icountsink+nptmass)  
  IF (tagged) read(10) tagi !skip tag
  read(10) spinady(icountsink+1:icountsink+nptmass)  
  IF (tagged) read(10) tagi !skip tag
  read(10) spinadz(icountsink+1:icountsink+nptmass)
  
!I'm just not sure about this point onwards - if RT comes into things
! (and it probably won't) then check it. But that won't be the case.
  
  !print *, "numssink(6) = ", numssink(6)
  DO i = 1, numssink(6)-9
    read(10)
  END DO
  DO i = 1, numssink(8)
    read(10)
  END DO
  
  !Read RT data if it's there
  !NB Original indices were from 1 to nelement; I'm guessing they should be changed.
  !Also NB - is this old? Some parts don't work logically...
!!$  IF (numberarray == 3) THEN
!!$    print *, "Reading RT", numberarray, numsRT(1:8)
!!$    IF (.NOT.allocated(radneigh)) THEN
!!$      CALL allocate_arrays_RT
!!$      print*, 'RT arrays allocated'
!!$    ELSE IF (.NOT. contiguous .OR. iblock == 1) THEN       
!!$      CALL reallocate_arrays_RT
!!$      print*, 'RT arrays reallocated'
!!$    END IF
!!$    IF (numsRT(3) == 1) read(10) radneigh(icount+1:icount+nelement)
!!$    read(10) e(icount+1:icount+nelement)
!!$    read(10) rkappa(icount+1:icount+nelement)
!!$    read(10) cv(icount+1:icount+nelement)
!!$    read(10) rlambda(icount+1:icount+nelement)
!!$    read(10) edd(icount+1:icount+nelement)
!!$    IF (numsRT(6) == 8) THEN
!!$      DO i = 1, 3
!!$        read(10) force(i,icount+1:icount+nelement)
!!$      END DO
!!$    END IF
!!$    IF (numsRT(7) == 1) read(10) dlnTdlnP(icount+1:icount+nelement)
!!$    IF (numsRT(7) == 13 .OR. numsRT(7) == 14) THEN
!!$      read(10) dlnTdlnP(icount+1:icount+nelement)
!!$      IF (numsRT(7) == 11) read(10) adiabaticgradient(icount+1:icount+nelement)
!!$      DO i = 1, 3
!!$        read(10) pressure(i,icount+1:icount+nelement)
!!$      END DO
!!$      DO i = 1, 3
!!$        read(10) viscosity(i,icount+1:icount+nelement)
!!$      END DO
!!$      DO i = 1, 3
!!$        read(10) gravity(i,icount+1:icount+nelement)
!!$      END DO
!!$      DO i = 1, 3
!!$        read(10) radpres(i,icount+1:icount+nelement)
!!$      END DO
!!$    END IF
!!$  END IF
  
  IF (contiguous) THEN
    icount = icount + nelement
    icountsink = icountsink + icountsink
  ELSE
    icount = 0
    icountsink = 0
  END IF
  iblock = iblock + 1
  
  RETURN
  END SUBROUTINE read_sphng_block
    
  
!---------------------------------------------------------------
!These create and destroy all the module's arrays, except the sink arrays and nelementblocks.
  
  SUBROUTINE allocate_arrays

    use sphdata
    use tachedata, only: nelement

  IMPLICIT NONE
  integer :: nalloc
  IF (contiguous .AND. nblocks > 1) THEN
    nalloc = nelementtot
  ELSE
    nalloc = nelement
  END IF
  print *, "ALLOCATING HYDRO ARRAYS"
  allocate( isteps(nalloc), iphase(nalloc), iunique(nalloc) )
  allocate( xyzmh(5,nalloc), vxyzu(4,nalloc), alphaMM(nalloc), rho(nalloc) )
  allocate(poten(nalloc))
  IF (igradh == 1) THEN
    allocate( gradh(nalloc),gradhsoft(nalloc) )
  ELSE
    allocate( dgrav(nalloc) )
  END IF
  RETURN
  END SUBROUTINE allocate_arrays
  
  SUBROUTINE deallocate_arrays
    use sphdata
  IMPLICIT NONE
  print *, "DEALLOCATING HYDRO ARRAYS"
  deallocate( isteps, iphase, iunique )
  deallocate( xyzmh, vxyzu, alphaMM, rho )
  deallocate(poten)
  IF (igradh == 1) THEN
    deallocate( gradh, gradhsoft )
  ELSE
    deallocate( dgrav )
  END IF
  RETURN
  END SUBROUTINE deallocate_arrays
  
  SUBROUTINE reallocate_arrays
  IMPLICIT NONE
  CALL deallocate_arrays
  CALL allocate_arrays
  RETURN
  END SUBROUTINE reallocate_arrays
  
  SUBROUTINE allocate_arrays_RT
    use sphdata
    use tachedata, only:nelement
  IMPLICIT NONE
  integer :: nalloc
  IF (contiguous .AND. nblocks > 1) THEN
    nalloc = nelementtot
  ELSE
    nalloc = nelement
  END IF
  print *, "Allocating RT arrays"
  allocate( radneigh(nalloc) )
  allocate( e(nalloc), rkappa(nalloc), cv(nalloc), rlambda(nalloc), edd(nalloc) )
  allocate( force(3,nalloc) )
  allocate( dlnTdlnP(nalloc), adiabaticgradient(nalloc) )
  allocate( pressure(3,nalloc), viscosity(3,nalloc), gravity(3,nalloc), radpres(3,nalloc) )
  RETURN
  END SUBROUTINE allocate_arrays_RT
  
  SUBROUTINE deallocate_arrays_RT
    use sphdata
  IMPLICIT NONE
  print *, "Deallocating RT arrays"
  deallocate( radneigh )
  deallocate( e, rkappa, cv, rlambda, edd )
  deallocate( force )
  deallocate( dlnTdlnP, adiabaticgradient )
  deallocate( pressure, viscosity, gravity, radpres )
  RETURN
  END SUBROUTINE deallocate_arrays_RT
  
  SUBROUTINE reallocate_arrays_RT
  IMPLICIT NONE
  CALL deallocate_arrays_RT
  CALL allocate_arrays_RT
  RETURN
  END SUBROUTINE reallocate_arrays_RT
  
  !Deallocate all arrays held in this module, including nelementblocks.
  SUBROUTINE deallocate_all_arrays
    use sphdata
  IMPLICIT NONE
  print *, "Deallocating all arrays:"
  IF (allocated(isteps)) CALL deallocate_arrays
  IF (allocated(radneigh)) CALL deallocate_arrays_RT
  IF (allocated(nelementblocks)) deallocate(nelementblocks)
  RETURN
  END SUBROUTINE deallocate_all_arrays

