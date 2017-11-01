SUBROUTINE calc_grav_brute
  ! Subroutine calculates gravitational forces and potentials for all particles
  ! It uses brute force to calculate forces and potentials for particles outside
  ! the smoothing volume

  use sphgravdata
  use treedata
  use sphkerneldata

  implicit none

  real,dimension(3) :: gravi

  integer :: ipart, jpart,k,ix
  real :: meanpot,sdpot, percent,counter

  allocate(gravxyz(3,npart))
  allocate(poten(npart))

  gravxyz(:,:) = 0.0
  poten(:) = 0.0

  isoft = 0

  percent = 0.0
  counter = 1.0
  do ipart = 1,npart

     percent = REAL(ipart)/REAL(npart)*100.0

     if(percent>counter)then
        print*, counter,'% complete'
        counter = counter +1.0
     endif

     gravi(:) = 0.0
     poteni = 0.0

     !$OMP PARALLEL &
     !$OMP shared(ipart, npart,poten,gravxyz) &
     !$OMP private(jpart,poteni,gravi,ix) 
     !$OMP DO SCHEDULE(runtime)
     do jpart = 1,npart

        if(ipart==jpart) cycle

        call particle_forces(ipart,jpart,poteni,gravi)

        poten(ipart) = poten(ipart) + poteni

        DO ix=1,3
           gravxyz(ix,ipart) = gravxyz(ix,ipart) + gravi(ix)
        ENDDO

     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     ! End of loop over jpart



     !  Add contribution to the potential from pointmasses

     do k=1,nptmass

        jpart = listpm(k)

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
END SUBROUTINE calc_grav_brute


