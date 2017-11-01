SUBROUTINE calc_grav_from_pot
  ! Subroutine calculates gravitational forces using particle potentials for all particles

  use sphgravdata
  use treedata
  use sphkerneldata

  implicit none

  real,dimension(3) :: dr

  integer :: ipart, jpart,k,ix
  real :: v, meanpot,sdpot, percent,counter
  real :: rhoj,grpm

  allocate(gravxyz(3,npart))
  gravxyz(:,:) = 0.0
  isoft = 0

  percent = 0.0
  counter = 1.0
  do ipart = 1,npart
    
     percent = REAL(ipart)/REAL(npart)*100.0

     if(percent>counter)then
        print*, counter,'% complete'
        counter = counter +1.0
     endif   

     !$OMP PARALLEL &
     !$OMP shared(ipart, npart,poten,gravxyz) &
     !$OMP private(jpart,ix) 
     !$OMP DO SCHEDULE(runtime)

     do k = 1, nneigh(ipart)

        jpart = neighb(ipart,k)
        
        ! Calculate gradient of SPH kernel
        hmean = 0.5*(xyzmh(5,ipart) + xyzmh(5,jpart))
        hmean21 = 1./(hmean*hmean)
        hmean41 = hmean21*hmean21
        pmassj = xyzmh(4,jpart)
        rhoj = rho(jpart)

        ! Separation of particles
        do ix = 1,3
            dr(ix) = xyzmh(ix,ipart) - xyzmh(ix,jpart)
        enddo

        rij2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

        v2 = rij2*hmean21
        v = rij/hmean

        index = v2/dvtable

	if(index<1) index=1

        dxx = v2 - index*dvtable
        index1 = index + 1

        IF (index1.GT.itable) then
           index = itable-1
           index1 = itable
        endif

        IF(index1 < 1 .or. index<1) THEN
           print*, jpart,ipart, index, itable, rij2, v2,dvtable, v2/dvtable, index1
        endif
        IF (v > radkernel) cycle
     
        dgrwdx = (grwij(index1) - grwij(index))/dvtable
        grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
        grpm = pmassj*grwtij

        ! Use gradient of kernel to calculate F_i = - grad_i phi (x,y,z)
        do ix = 1,3            
                gravxyz(ix,ipart) = gravxyz(ix,ipart) - grpm*dr(ix)*poten(jpart)/rhoj          
        enddo

     enddo
     !$OMP END DO
     !$OMP END PARALLEL

     !  Add contribution to the potential from pointmasses

     do k=1,nptmass

        jpart = listpm(k)

        print*, listpm
        if(jpart==ipart) cycle

        sep = (xyzmh(1,ipart) - xyzmh(1,jpart))**2 + &
             (xyzmh(2,ipart) - xyzmh(2,jpart))**2 +&
             (xyzmh(3,ipart) - xyzmh(3,jpart))**2

        sep =sqrt(sep)

        poten(ipart) = poten(ipart) + xyzmh(4,jpart)/sep

     enddo
     !print*, 'Particle ', ipart, ': Potential - ', poten(ipart), 'Force: ',gravi(:)
  ENDDO
  ! End loop over all particles

  ! Calculate mean and standard deviation of potential

  meanpot = sum(poten)/REAL(npart)
  sdpot = 0.0

  !$OMP PARALLEL &
  !$OMP shared(nneigh,meanneigh,npart)&
  !$OMP private(ipart) &
  !$OMP reduction(+:sdneigh)
  !$OMP DO SCHEDULE(runtime)
  do ipart=1,npart
     sdpot = sdpot+(poten(ipart)-meanpot)**2
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  sdpot = sqrt(sdpot/REAL(npart))

  print*, 'Mean potential is ', meanpot
  print*, 'Standard Deviation: ', sdpot


  return
END SUBROUTINE calc_grav_from_pot


