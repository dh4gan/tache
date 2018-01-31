subroutine read_dump(ifile,skipdump)
!
! Subroutine controls the reading in of all data
! (SPH, grid or mesh)
!

use tachedata
implicit none

integer,intent(in) :: ifile
logical,intent(inout) :: skipdump
integer :: check,smallfile

skipdump = .false.

check =0
smallfile = 0

if(filetype=='SPH') then

   print*, 'Reading SPH file'
  
   if(fileformat=='sphNG_wkmr') then

      call rdump_sphNG_wkmr(filename(ifile),check,smallfile)

   else if(fileformat=='sphNG_iab') then
      call rdump_sphNG_iab(filename(ifile),check,smallfile)

   endif


 ! Skip small dumps
     if(smallfile/=0) THEN
        print*, 'Skipping small/incomplete dump'
        skipdump = .true.
     ENDIF

     if (check /= 0) then
        !	If file missing in series, skip it
        IF(ifile==1) THEN
           print*, 'ERROR: Program aborted at first file read in'
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
           stop
        ELSE
           print*, 'Skipping missing dump ',filename(ifile)
           skipdump = .true.
        ENDIF
     endif

endif

end subroutine read_dump
