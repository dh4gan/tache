      subroutine wdump_grav(gravfile,potfile)

      use sphdata
      use tachedata, only: nelement

      implicit none
      integer :: i
      character(7) :: gravfile
      character(6) :: potfile

! Open and write gravity file

      print*, 'Gravitational Force file: ', gravfile

      open(15,file=gravfile, form="unformatted",status="unknown")

      print*, '      - Writing to file'

      write(15,err=100) (gravxyz(1,i), i=1,nelement)
      write(15,err=100) (gravxyz(2,i), i=1,nelement)
      write(15,err=100) (gravxyz(3,i), i=1,nelement)
      write(15,err=100) (iphase(i), i=1,nelement)

      close(15)

      print*, '      - Done'

! Open and read Potential file
      print*, 'Gravitational Potential file: ', potfile

      open(15,file=potfile, form="unformatted",status="unknown")

      print*, '      - Writing to file'

    write(15,err=200)(poten(i), i=1,nelement)
    write(15,err=200)(iphase(i), i=1,nelement)

      close(15)

      return
100 print*, "Gravitational Force File Write Failed"
      return
200 print*, "Gravitational Potential File Write Failed"

      return
      end subroutine wdump_grav
