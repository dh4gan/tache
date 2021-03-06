SUBROUTINE calc_grav_brute
  ! Subroutine calculates gravitational forces and potentials for all particles
  ! It uses brute force O(N^2) to calculate forces and potentials 
  ! for particles outside the smoothing volume

  use sphdata
  use sphneighbourdata
  use sphkerneldata
  use tachedata, only:nelement

  implicit none

  real,dimension(3) :: gravi

  integer :: ielement, jpart,k,ix
  real :: meanpot,sdpot, percent,increment

  allocate(gravxyz(3,nelement))
  allocate(poten(nelement))

  gravxyz(:,:) = 0.0
  poten(:) = 0.0

  isoft = 0

  percent = 0.0
  increment = 1.0

  do ielement = 1,nelement

     call element_percent_complete(ielement,nelement,percent,increment)

     gravi(:) = 0.0
     poteni = 0.0

     ! Loop over all elements

     !$OMP PARALLEL &
     !$OMP shared(ielement, nelement,poten,gravxyz) &
     !$OMP private(jpart,poteni,gravi,ix) 
     !$OMP DO SCHEDULE(runtime)
     do jpart = 1,nelement

        if(ielement==jpart) cycle

        ! Compute the forces between ielement and jpart
        call particle_forces(ielement,jpart,poteni,gravi)

        poten(ielement) = poten(ielement) + poteni

        DO ix=1,3
           gravxyz(ix,ielement) = gravxyz(ix,ielement) + gravi(ix)
        ENDDO

     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     ! End of loop over jpart


     !  Add contribution to the potential from pointmasses

     do k=1,nptmass

        jpart = listpm(k)

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
END SUBROUTINE calc_grav_brute


