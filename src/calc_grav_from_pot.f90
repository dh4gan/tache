SUBROUTINE calc_grav_from_pot
  !*************************************************************
  ! Subroutine calculates gravitational forces of SPH particles
  ! using particle potential
  !*************************************************************

  use sphdata
  use sphneighbourdata
  use sphkerneldata
  use tachedata, only: nelement

  implicit none

  real,dimension(3) :: dr

  integer :: ielement, jpart,k,ix
  real :: v, meanpot,sdpot, percent,counter
  real :: rhoj,grpm

  allocate(gravxyz(3,nelement))
  gravxyz(:,:) = 0.0
  isoft = 0

  percent = 0.0
  counter = 1.0
  do ielement = 1,nelement
    
     percent = REAL(ielement)/REAL(nelement)*100.0

     if(percent>counter)then
        print*, counter,'% complete'
        counter = counter +1.0
     endif   

     ! Loop over all neighbours of ielement
     !$OMP PARALLEL &
     !$OMP shared(ielement, nelement,poten,gravxyz) &
     !$OMP private(jpart,ix) 
     !$OMP DO SCHEDULE(runtime)

     do k = 1, nneigh(ielement)

        jpart = neighb(ielement,k)
        
        ! Calculate gradient of SPH kernel
        hmean = 0.5*(xyzmh(5,ielement) + xyzmh(5,jpart))
        hmean21 = 1./(hmean*hmean)
        hmean41 = hmean21*hmean21
        pmassj = xyzmh(4,jpart)
        rhoj = rho(jpart)

        ! Separation of particles
        do ix = 1,3
            dr(ix) = xyzmh(ix,ielement) - xyzmh(ix,jpart)
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
           print*, jpart,ielement, index, itable, rij2, v2,dvtable, v2/dvtable, index1
        endif
        IF (v > radkernel) cycle
     
        dgrwdx = (grwij(index1) - grwij(index))/dvtable
        grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
        grpm = pmassj*grwtij

        ! Use gradient of kernel to calculate F_i = - grad_i phi (x,y,z)
        do ix = 1,3            
                gravxyz(ix,ielement) = gravxyz(ix,ielement) - grpm*dr(ix)*poten(jpart)/rhoj          
        enddo

     enddo
     !$OMP END DO
     !$OMP END PARALLEL

     !  Add contribution to the potential from pointmasses

     do k=1,nptmass

        jpart = listpm(k)

        print*, listpm
        if(jpart==ielement) cycle

        sep = (xyzmh(1,ielement) - xyzmh(1,jpart))**2 + &
             (xyzmh(2,ielement) - xyzmh(2,jpart))**2 +&
             (xyzmh(3,ielement) - xyzmh(3,jpart))**2

        sep =sqrt(sep)

        poten(ielement) = poten(ielement) + xyzmh(4,jpart)/sep

     enddo
     !print*, 'Particle ', ielement, ': Potential - ', poten(ielement), 'Force: ',gravi(:)
  ENDDO
  ! End loop over all particles

  ! Calculate mean and standard deviation of potential

  meanpot = sum(poten)/REAL(nelement)
  sdpot = 0.0

  !$OMP PARALLEL &
  !$OMP shared(nneigh,meanneigh,nelement)&
  !$OMP private(ielement) &
  !$OMP reduction(+:sdneigh)
  !$OMP DO SCHEDULE(runtime)
  do ielement=1,nelement
     sdpot = sdpot+(poten(ielement)-meanpot)**2
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  sdpot = sqrt(sdpot/REAL(nelement))

  print*, 'Mean potential is ', meanpot
  print*, 'Standard Deviation: ', sdpot


  return
END SUBROUTINE calc_grav_from_pot


