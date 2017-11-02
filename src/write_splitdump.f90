subroutine write_splitdump(ifile)

use tachedata
use sphdata

implicit none

integer, intent(in) :: ifile
integer :: i,iclass,ipart,check,counter,nparthold,nptmasshold
integer,allocatable,dimension(:) :: isorthold,iphasehold,istepshold
real,allocatable,dimension(:) :: rhohold,dgravhold
real,allocatable,dimension(:,:) :: xyzmhhold,vxyzuhold

real, dimension(3) :: eigensingle
character(100) :: classfile

!
! Splits up dump into four dumps according to element classification
! (Currently only writes SPH dumps in sphNG_wkmr format)
!

  prefixes(1) = "clusters"
  prefixes(2) = "filaments"
  prefixes(3) = "sheets"
  prefixes(4) = "voids"


!***********************************************************
! 1. Classify each particle according to its eigenvalues
!***********************************************************

print*, 'Classifying dump by eigenvalue: threshold=',threshold

class(:) = -1
classnum(:) = 0

do i=1,neigen
   ielement = eigenelement(i)
   eigensingle(:) = eigenvalues(:,i)
   CALL classify_by_eigenvalues(class(ielement), eigensingle,threshold)

enddo


if(filetype=='SPH') then

   !************************************************************
   ! 2a. Write classes to single file (for particle tracking)
   !************************************************************

   print*, 'Writing all particle classifications to file ', TRIM(memberfile(ifile))
   open(12,file=memberfile(ifile),form='unformatted')

   write(12) nelement
   write(12) gt
   write(12) (class(i),i=1,nelement)
   
   close(12)


   !*************************************************************
   ! 2b. Write to separate SPH dumps
   !*************************************************************

   print*, 'Now writing different classifications to separate SPH files'

   !	Write all SPH data to holding array

   nparthold = nelement
   nptmasshold = nptmass
 
   allocate(isorthold(nelement))
   allocate(iphasehold(nelement))
   allocate(istepshold(nelement))
   allocate(xyzmhhold(5,nelement))
   allocate(vxyzuhold(4,nelement))
   allocate(rhohold(nelement))
   allocate(dgravhold(nelement))
   
   !allocate(listpmhold(nptmass))
   !allocate(spinxhold(nptmass))
   !allocate(spinyhold(nptmass))
   !allocate(spinzhold(nptmass))
   !allocate(angaddxhold(nptmass))
   !allocate(angaddyhold(nptmass))
   !allocate(angaddzhold(nptmass))
   !allocate(spinadxhold(nptmass))
   !allocate(spinadyhold(nptmass))
   !allocate(spinadzhold(nptmass))
   
   
   isorthold(:) = isort(:)
   iphasehold(:) = iphase(:)		
   istepshold(:) = isteps(:)
   
   xyzmhhold(:,:) = xyzmh(:,:)
   vxyzuhold(:,:) = vxyzu(:,:)
   rhohold(:) = rho(:)
   dgravhold(:) = dgrav(:)
   
   ! listpmhold(:) = listpm(:)
   ! spinxhold(:) = spinx(:)
   ! spinyhold(:) = spiny(:)
   ! spinzhold(:) = spinz(:)
   
   ! angaddxhold(:) = angaddx(:)
   ! angaddyhold(:) = angaddy(:)
   ! angaddzhold(:) = angaddz(:)
   
   ! spinadxhold(:) = spinadx(:)
   ! spinadyhold(:) = spinady(:)
   ! spinadzhold(:) = spinadz(:)
   

   do iclass=1,nclasses
          
      classfile=trim(filename(ifile))//'_'//trim(prefixes(iclass))

      npart = classnum(iclass)
      
      print*, 'Classification ', prefixes(iclass), ' being written to ',TRIM(classfile)
      print*, 'Particle count ',npart
      
      nptmass=0	
      n1 = npart
      counter=0
      
      do ipart=1,nparthold
         
         IF(class(ipart)/=iclass) cycle
         
         counter = counter +1
         
         iphase(counter) = iphasehold(ipart)
         isteps(counter) = 0
         isort(counter) = counter
         
         xyzmh(:,counter) = xyzmhhold(:,ipart)
         vxyzu(:,counter) = vxyzuhold(:,ipart)
         rho(counter) = rhohold(ipart)
         dgrav(counter) = dgravhold(ipart)			
         
      enddo
      
      ! If requested, recalculate density (using same h)
      
      if(density_recalc=='y') call recalc_density
      
      !	Call wdump with these new arrays
      
      check =0
      call wdump(classfile, check)
      
   enddo

 print*, 'Classification files written'

 !**********************************************************
 ! 4. Deallocate arrays, ready for the next loop
 !**********************************************************
 
 deallocate(class)

endif

end subroutine write_splitdump
