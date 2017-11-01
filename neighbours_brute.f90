subroutine neighbours_brute(sphfile)
  ! Subroutine finds nearest neighbours in radius 2h, using the octree
  ! h already defined for all particles
  ! This counts accreted particles and pointmasses too
  use sphgravdata
  use treedata

 integer :: i,j

  real :: hi,hj, hmean, sep,percent
  real, parameter :: tiny = 1.0e-34

  character(7) :: sphfile
  character(18) :: neighbours

  logical :: existneigh

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

IF(existneigh.eqv..false.) then
  print*, 'Beginning brute force search for neighbours'
  print*, "------------------------------------------------------------------------"

  ! Parallelisable:
  ! Shared: iphase, xyzmh, nneigh, neighb, neighmax, tiny
  ! Private: i,j,hi,hj,hmean,sep
  !$OMP PARALLEL &
  !$OMP shared(npart, iphase,xyzmh, nneigh,neighb) &
  !$OMP private(i,j,hi,hj,hmean,sep)
  !$OMP DO SCHEDULE(runtime)
  do i=1,npart

      IF(iphase(i)<0) cycle
      percent = 100.0*REAL(i)/REAL(npart)
      if(MOD(real(i), real(npart)/10.0)<1.0) THEN
          write(*,'(I3," % complete")') INT(percent)
      ENDIF

      do j=1,npart

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


  ! Write the neighbour data to file
  call write_neighbours(neighbours)

  ! End loop over all particles
endif

  return
  end subroutine neighbours_brute
