subroutine neighbours_grid(sphfile)
  ! Subroutine finds nearest neighbours in radius 2h, using the regular grid
  ! h already defined for all particles
  ! This counts accreted particles and pointmasses too

  use sphdata
  use sphneighbourdata
  IMPLICIT NONE

interface
   subroutine find_particles_in_range(xpart,ypart,zpart,range)
     real, intent(in) :: xpart,ypart,zpart,range
   end subroutine find_particles_in_range
end interface


  character(7) ::sphfile
  character(18) :: neighbours
  logical :: existneigh

  integer :: ielement,jelement, k

  real :: hi,hj, hmean, range, sep,percent,counter
  real, parameter :: tiny = 1.0e-34

  ! Check if a neighbourhood file currently exists. 
  print*, "-----------------------------------------------"
  print*, 'Searching for neighbourhood fileset'
  write(neighbours,'("neighbours_",A7)') sphfile


  ! Check if file is available
  INQUIRE(file=neighbours,exist=existneigh)

  IF(existneigh.eqv..true.) then
     print*, 'Neighbourhood file found: reading'

     call read_neighbours(neighbours)

  ENDIF



  !------------------------------------------------------------------------------

  IF(existneigh.eqv..false.) then
     print*, 'No satisfactory fileset found: creating neighbours from scratch'
     print*, "----------------------------------------------------------------"

     ! Loop over all particles

     percent = 0.0
     counter = 1.0

     do ielement = 1,nelement

        if(iphase(ielement)/=0) cycle
         ! printing the extent to which the process is complete
         percent = 100.0*REAL(ielement)/REAL(nelement)

        if(percent>counter)then
            print*, counter,'% complete'
            counter = counter +1.0
        endif


        ! Find all grid cells within tolerance*h of the particle

        range = tolerance*xyzmh(5,ielement)        
        hi = xyzmh(5,ielement)

        ! Create list of all cells in range of the particle
        ! And a list of all particles in these cells

        call find_particles_in_range(xyzmh(1,ielement),xyzmh(2,ielement),&
             xyzmh(3,ielement), range)

        ! Test for neighbourship all particles in the particle list

        !print*, ielement, ncellrange, nelementiclelist
        !$OMP PARALLEL &
        !$OMP shared(nelementiclelist,particlelist,iphase,hi,ielement)&
        !$OMP shared(xyzmh,nneigh,neighb) &
        !$OMP private(k,jelement,hj,hmean,sep)
        !$OMP DO SCHEDULE(runtime)
        do k = 1,nelementiclelist
           jelement = particlelist(k)
           !print*,ielement, k,nelementiclelist, jelement
           IF(ielement==jelement) cycle
           IF(iphase(jelement)/=0) cycle
           
           hj = xyzmh(5,jelement)
           
           if(ielement/=jelement) then

             hmean = (hi + hj)/2.0
        
              sep = (xyzmh(1,ielement) - xyzmh(1,jelement))**2 + &
                  (xyzmh(2,ielement) - xyzmh(2,jelement))**2 + &
                  (xyzmh(3,ielement) - xyzmh(3,jelement))**2 + tiny
              sep = sqrt(sep)
              !print*, xyzmh(:,ielement)
              !print*, xyzmh(:,jelement)
              !print*,ielement,jelement, member(ielement), member(jelement), hi, hj, sep/hmean
        
              !	if particle j in neighbour sphere, then add to neighbour list
              if(sep<2.0*hmean.and.nneigh(ielement)<neighmax) then
                 !$OMP CRITICAL
                 nneigh(ielement) = nneigh(ielement) + 1
                 neighb(ielement,nneigh(ielement)) = jelement
                 !$OMP END CRITICAL
              endif
           endif
           
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! End loop over particles in cell
        if(nneigh(ielement) < 20) then
           print*, 'nneigh for particle ',ielement,' is ',nneigh(ielement), neighmax
           print*, 'Particle type: ', iphase(ielement)
           print*, 'Number of cells/particles tested ', ncellrange, nelementiclelist
           !print*, 'xyzmh: ', xyzmh(:,ielement)
        endif
        
        
     enddo
     ! End of loop over all particles
  
     
     ! write neighbours to file.
     print*, 'Writing neighbours to file ', neighbours
     
     call write_neighbours(neighbours,tolerance)
     
  endif
  
  ! Calculate mean and standard deviation of neighbour counts
  meanneigh = sum(nneigh)/REAL(nelement)
  sdneigh = 0.0
  
  !$OMP PARALLEL &
  !$OMP shared(nneigh,meanneigh,nelement)&
  !$OMP private(ielement) &
  !$OMP reduction(+:sdneigh)
  !$OMP DO SCHEDULE(runtime)
  do ielement=1,nelement
     sdneigh = sdneigh+(nneigh(ielement)-meanneigh)**2
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  sdneigh = sqrt(sdneigh/REAL(nelement))
  
  print*, 'Mean neighbour number is ', meanneigh
  print*, 'Standard Deviation: ', sdneigh
  neighcrit = meanneigh-5.0*sdneigh
  
  print*, 'Clumps created only if neighbour number greater than ', neighcrit
  
  return
end subroutine neighbours_grid



SUBROUTINE find_particles_in_range(xpart,ypart,zpart,range)
! Subroutine finds all grid cells within a distance hpart from the co-ordinates
! xpart, ypart, zpart, and lists all particles in those cells

  use sphdata, only:nelement
  use sphneighbourdata

  IMPLICIT NONE

  real, intent(in) :: xpart,ypart,zpart,range
  integer :: icell, jcell, kcell, thiscell, istart,iend
  integer :: icellmin, icellmax, jcellmin, jcellmax, kcellmin, kcellmax
  integer :: nbinned, ielement, ientry

! Find cell indices containing x-range, x+range

  icellmin = int((xpart-range+xmax)/dgrid +1)
  icellmax = int((xpart+range+xmax)/dgrid +1)

  if(icellmin <1) icellmin = 1
  if(icellmax > ngridx) icellmax = ngridx

! Do same for y

  jcellmin = int((ypart-range+ymax)/dgrid +1)
  jcellmax = int((ypart+range+ymax)/dgrid +1)

  if(jcellmin <1) jcellmin = 1
  if(jcellmax > ngridy) jcellmax = ngridy

! Finally z

  kcellmin = int((zpart-range+zmax)/dgrid +1)
  kcellmax = int((zpart+range+zmax)/dgrid +1)

  if(kcellmin <1) kcellmin = 1
  if(kcellmax > ngridz) kcellmax = ngridz

! Create a list of cellIDs corresponding to these locations
! Count total particle number inside these cells

  nbinned = 0

  ncellrange = (icellmax-icellmin+1)*(jcellmax-jcellmin+1)*(kcellmax-kcellmin+1)

  if(allocated(cellist)) deallocate(cellist)
  allocate(cellist(ncellrange))

  ncellrange = 0

  do icell = icellmin, icellmax
     do jcell = jcellmin, jcellmax
        do kcell = kcellmin, kcellmax
           ncellrange = ncellrange +1
           call get_cellID(thiscell,icell,jcell,kcell,ngridx,ngridy,ngridz)
           cellist(ncellrange) = thiscell
           nbinned = nbinned + n_occ(thiscell)
        enddo
     enddo
  enddo

  if(allocated(particlelist)) deallocate(particlelist)
  allocate(particlelist(nbinned))

particlelist(:) = 0
nelementiclelist = 0

! Now find list of particles from each cell

do icell=1,ncellrange

   thiscell = cellist(icell)
   
   ! Find the beginning of the cell's members in the sorted particle list
   
   istart = 1

   ! Skip all the non-gas particles (find the first non-zero entry in isortcellpart)

   ientry = -10
   do while(ientry <0)

      ientry = isortcellpart(istart)
      istart= istart+1
   enddo

   do jcell=1,thiscell-1
      istart = istart + n_occ(jcell)
   enddo

   ! End is simply the number of particles in the cell
   iend = istart + n_occ(thiscell)-1
      
   if (iend> nelement) iend = nelement
   ! Now go through this section of the list and add particles
   do ielement=istart,iend
      nelementiclelist = nelementiclelist+1
      particlelist(nelementiclelist) = isortcellpart(ielement)      
   enddo

enddo


if(nelementiclelist/=nbinned) then
   print*, "WARNING! particlelist total doesn't match expected occupancy:",&
        nelementiclelist,nbinned
endif

!print*, 'Found ', nelementiclelist, ' particles in range from ',ncellrange, ' cells'

return

END SUBROUTINE find_particles_in_range
