subroutine checksetup
!
! Subroutine checks the inputs are correct
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

if(fileformat/='sphNG_wkmr'.and.fileformat/='sphNG_iab') then
   print*, 'Unrecognised SPH file format'
   print*, 'Halting program'
   STOP
endif


! Check tensorchoice
if(tensorchoice/='tidal'.and.tensorchoice/='velocity') then
   print*, 'WARNING: Unrecognised tensor choice ', tensorchoice
   print*, 'Assuming choice=velocity shear'
   tensorchoice='velocity'
endif

! Is splitdump y/n?

if(splitdumpchoice=='y' .or. splitdumpchoice=='Y') splitdump = .true.

if(splitdumpchoice/='y' .and. splitdumpchoice/='Y' .and. splitdumpchoice/='N' .and. splitdumpchoice/='n') then
   print*, 'WARNING: splitdump option must be (y/n)'
   print*, 'Assuming splitdump=n'
   splitdump=.false.
endif

if(splitdump.and.threshold<0.0) then
   print*, 'ERROR: threshold for classification must be positive: check tache.params'
   STOP
endif


end subroutine checksetup
