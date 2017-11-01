      subroutine rdump_grav(gravfile,potfile,rcheck)

      use sphdata

      implicit none
      integer :: ios,unit,rcheck,i
      integer(kind=1),allocatable,dimension(:) :: testph
      character(7) :: gravfile
      character(6) :: potfile

! Set unit value!
      unit = 1
      
! Open and read gravity file
      print*, '2a) Gravity file: ', gravfile

      allocate(gravxyz(3,nelement))
      allocate(testph(nelement)   )
      allocate(tphase(nelement)   )

      open(15,file=gravfile, form="unformatted",status="old", &
           action='read',iostat=ios)
      if(ios.ne.0) then
         print*, '      - ERROR! Gravity input file not found!'
         rcheck = 1
         return
      endif
      print*, '      - Reading in Gravity file'

      read(15,err=100) (gravxyz(1,i), i=1,nelement)
      read(15,err=100) (gravxyz(2,i), i=1,nelement)
      read(15,err=100) (gravxyz(3,i), i=1,nelement)
      read(15,err=100) (testph(i), i=1,nelement)

      do i=1,nelement
         tphase(i) = unit*testph(i)
      enddo

      close(15)
      deallocate(testph)

      print*, '      - Gravity file correctly read in'
      print*, ' '

! Open and read Potential file
      print*, '2b) Potential file: ', potfile

      allocate(poten(nelement))
      allocate(testph(nelement)   )
      allocate(uphase(nelement)   )

      open(15,file=potfile, form="unformatted",status="old", &
           action='read',iostat=ios)
      if(ios.ne.0) then
         print*, '      - ERROR! Potential input file not found!'
         rcheck = 1
         return
      endif
      print*, '      - Reading in Potential file'

      read(15,err=200) (poten(i), i=1,nelement)
!      read(15,err=200) (testph(i), i=1,nelement)

     do i=1,nelement
!         uphase(i) = unit*testph(i)
        uphase(i) = tphase(i)
      enddo
      close(15)
      deallocate(testph)

      print*, '      - Potential file correctly read in'
      print*, ' '
      return 
 100  print*, '      - Gravity file read failed'
      rcheck=1
      return
 200  print*, '      - Potential file read failed'
      rcheck = 1

      return
      end subroutine rdump_grav
