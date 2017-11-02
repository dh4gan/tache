subroutine checksetup
!
! TODO - Subroutine checks the inputs are correct
!

use tachedata
implicit none

logical :: existlist


! Does the filename list exist?

inquire(file=listfile,exist= existlist)

if(existlist.eqv..false.) then
   print*, 'ERROR: List of filenames ',listfile,' not found'
   STOP
endif

! For now, code only does SPH data - enforce this
if(filetype/='SPH') then
   print*, 'WARNING: TACHE does not yet analyse ', filetype, ' data'
   print*, 'Halting program'
   STOP
endif

! Check fileformat is OK

! Check tensorchoice


! Is threshold a number?


! Is splitdump y/n?

if(splitdump/='y' .and. splitdump/='Y' .and. splitdump/='N' .and. splitdump/='n') then
   print*, 'WARNING: splitdump option must be (y/n)'
   print*, 'Assuming splitdump=n'
   splitdump='n'
endif


end subroutine checksetup
