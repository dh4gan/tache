SUBROUTINE make_grid(sphfile, hmean,dgridmin)
  ! Written 17/2/15 by dh4gan
  ! This subroutine makes a regular grid in x,y,z and 
  ! bins the particles accordingly

  ! x face = -xmax +(i-1)*dgrid and so on for y,z
  ! x centre = -xmax +(i-1)*dgrid +0.5*dgrid and so on for y,z
  
  ! Find a particle's cell, location (xi,yi,zi):
  
  ! icell = int((xi-xmin)/dgrid +1) = int((xi+xmax)/dgrid +1)
  ! (and so on for jcell, kcell (y,z)

  USE sphgravdata
  USE treedata
  
  IMPLICIT NONE  

  character(7) :: sphfile

  integer :: ielement, ix,iy,iz, icell,jcell,kcell
  integer :: thiscell, counter, celltot, ngas

  real :: hmean,dgridmin
  
  ! Set up grid in x,y,z

  ! Find maximum and minimum values for all coordinates
  xmax = 0.0
  ymax = 0.0
  zmax = 0.0

  allocate(isort(nelement))
  allocate(iorig(nelement))

  ngas = 0
  do ielement=1,nelement
     ! rho(ielement) = xyzmh(4,ielement)/(xyzmh(5,ielement)**3) ! DEBUG LINE
     isort(ielement) = ielement
     iorig(ielement) = ielement
     
     if(iphase(ielement)==0) ngas = ngas +1
     
     IF(abs(xyzmh(1,ielement))> xmax) xmax = abs(xyzmh(1,ielement))
     IF(abs(xyzmh(2,ielement))> ymax) ymax = abs(xyzmh(2,ielement))
     IF(abs(xyzmh(3,ielement))> zmax) zmax = abs(xyzmh(3,ielement))
     
  enddo
  
  print*,'-----------------------------------------'
  print*, "Maximum values for spatial co-ordinates: "
  print*, "x: ", xmax
  print*, "y: ", ymax
  print*, "z: ", zmax
  
  
  ! Grid spacing
  dgrid = 2.0*hmean
  
  if(dgrid< dgridmin) then
     print*, 'Calculated dgrid too small: ', hmean, dgrid
     dgrid = dgridmin
  endif

  ngridx = int(2.0*xmax/dgrid +1)
  ngridy = int(2.0*ymax/dgrid +1)
  ngridz = int(2.0*zmax/dgrid +1)

  ncells = ngridx*ngridy*ngridz

  print*, 'xmax, ymax, zmax: ', xmax, ymax, zmax
  print*, 'Setting up regular grid, spacing: ',dgrid
  print*, 'Grid dimensions:'
  print*, 'nx = ',ngridx
  print*, 'ny = ', ngridy
  print*, 'nz = ', ngridz
  print*, 'Total cell no.: ',ncells

!!$  ! Set up grid IDs
!!$
!!$  counter = 0
!!$
!!$  !allocate(cellID(ngridx,ngridy,ngridz))
!!$
!!$  do ix=1,ngridx
!!$     do iy = 1,ngridy
!!$        do iz = 1,ngridz
!!$           counter =counter +1
!!$           cellID(ix,iy,iz) = counter
!!$        enddo
!!$     enddo
!!$  enddo

  ! Set up grid variables

  allocate(n_occ(ncells))

  print*, 'Cell IDs and occupancy set up'

  n_occ(:) = 0

  ! Bin particles in grid
  !$OMP PARALLEL &
  !$OMP shared(nelement,iphase,xyzmh,cellID,n_occ,occ) &
  !$OMP private(ielement,icell,jcell,kcell,thiscell)
  !$OMP DO SCHEDULE(runtime)
  do ielement=1,nelement
     
     if(iphase(ielement)/=0) cycle

     icell = int((xyzmh(1,ielement)+xmax)/dgrid +1)
     jcell = int((xyzmh(2,ielement)+ymax)/dgrid +1)
     kcell = int((xyzmh(3,ielement)+zmax)/dgrid +1)

     if(icell*jcell*kcell==0) print*, xyzmh(1:3,ielement), icell, jcell, kcell
     !thiscell = cellID(icell,jcell,kcell)

     call get_cellID(thiscell,icell,jcell,kcell,ngridx,ngridy,ngridz)

     if(thiscell>ncells) then
        print*, 'Error: cellID exceeds cell total'
        print*, thiscell, ncells, icell,jcell,kcell
     endif

     !$OMP CRITICAL
     n_occ(thiscell) = n_occ(thiscell)+1
     !$OMP END CRITICAL

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  cellmax = maxval(n_occ)

  ! Sum n occ to confirm that all particles are binned

  celltot = sum(n_occ)

  print*, 'Maximum cell occupancy: ',cellmax

  if(celltot.ne.nelement) then
     print*, "WARNING: Binned total doesn't match particle total: ", celltot, nelement
  endif

  !allocate(occ(ncells,cellmax))
  allocate(member(nelement))

  print*, 'Cell occupancy allocation complete'
  n_occ(:) = 0
  !occ(:,:) = 0

  ! Bin particles in grid
  !$OMP PARALLEL &
  !$OMP shared(nelement,iphase,xyzmh,cellID,n_occ,occ) &
  !$OMP private(ielement,icell,jcell,kcell,thiscell)
  !$OMP DO SCHEDULE(runtime)
  do ielement=1,nelement
     
     if(iphase(ielement)/=0) cycle

     ngas = ngas +1
     icell = int((xyzmh(1,ielement)+xmax)/dgrid +1)
     jcell = int((xyzmh(2,ielement)+ymax)/dgrid +1)
     kcell = int((xyzmh(3,ielement)+zmax)/dgrid +1)

     if(icell*jcell*kcell==0) print*, xyzmh(1:3,ielement), icell, jcell, kcell
!     thiscell = cellID(icell,jcell,kcell)

     call get_cellID(thiscell,icell,jcell,kcell,ngridx,ngridy,ngridz)

     if(ielement==1) print*, 'particle ',ielement, ' in cell ', thiscell
     !$OMP CRITICAL
     n_occ(thiscell) = n_occ(thiscell)+1
!     occ(thiscell,n_occ(thiscell)) = ielement
     member(ielement) = thiscell
     !$OMP END CRITICAL
     if(ielement==1) print*, 'membership: ', member(ielement)
     
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
     
! Finally, create a sorted list of particle IDs so that we can access lists of
! cell members in cell order


print*, 'Sorting so that membership can be accessed in cell order'
allocate(isortcellpart(nelement))

do ielement=1,nelement
   isortcellpart(ielement) = ielement
   if(iphase(ielement)/=0)   isortcellpart(ielement) = -10
enddo

call sort2(nelement,member,isortcellpart,nelement)

print*, 'Sort complete'


return
end subroutine make_grid
