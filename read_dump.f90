subroutine read_dump(filename,skipdump)
!
! Subroutine controls the reading in of all data
! (SPH, grid or mesh)
!

skipdump = .false.

check =0
smallfile = 0

if(filetype=='SPH') then

   if(fileformat=='sphNG_wkmr') then

      call rdump_sphNG_wkmr(filename,check,smallfile)

   else if(fileformat=='sphNG_iab') then
      call rdump_sphNG_iab(filename,check,smallfile)

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
           stop
        ELSE
           print*, 'Skipping missing dump ',filename(ifile)
           skipdump = .true.
        ENDIF
     endif

endif

end subroutine read_dump
