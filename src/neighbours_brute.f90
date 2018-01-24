subroutine neighbours_brute
  ! Subroutine finds nearest neighbours in radius 2h, using the octree
  ! h already defined for all particles
  ! This counts accreted particles and pointmasses too
  use sphdata
  use sphneighbourdata
  use tachedata, only: nelement

  implicit none

 integer :: i,j

  real :: hi,hj, hmean, sep,percent
  real, parameter :: tiny = 1.0e-34

  print*, 'Beginning brute force search for neighbours ',nelement
  print*, "------------------------------------------------------------------------"

  !$OMP PARALLEL &
  !$OMP shared(nelement, iphase,xyzmh, nneigh,neighb) &
  !$OMP private(i,j,hi,hj,hmean,sep)
  !$OMP DO SCHEDULE(runtime)
  do i=1,nelement

      IF(iphase(i)<0) cycle
      percent = 100.0*REAL(i)/REAL(nelement)
      if(MOD(real(i), real(nelement)/10.0)<1.0) THEN
          write(*,'(I3," % complete")') INT(percent)
      ENDIF

      do j=1,nelement

          IF(i==j) cycle
          IF(iphase(j)<0) cycle
     
          hi = xyzmh(5,i)
          hj = xyzmh(5,j)
     
          if(i/=j) then
              hmean = (hi + hj)/2.0
        
              sep = (xyzmh(1,i) - xyzmh(1,j))**2 + &
                  (xyzmh(2,i) - xyzmh(2,j))**2 + &
                  (xyzmh(3,i) - xyzmh(3,j))**2 + tiny
              sep = sqrt(sep)
        
              !	if particle j in neighbour sphere, then add to neighbour list
              if(sep<2.0*hmean.and.nneigh(i)<neighmax) then
                 !$OMP CRITICAL
                  nneigh(i) = nneigh(i) + 1
                  neighb(i,nneigh(i)) = j
                  !$OMP END CRITICAL
              endif
          endif
     
      enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


  return
  end subroutine neighbours_brute
