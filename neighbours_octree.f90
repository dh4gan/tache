subroutine neighbours_octree(sphfile)
  ! Subroutine finds nearest neighbours in radius 2h, using the octree
  ! h already defined for all particles
  ! This counts accreted particles and pointmasses too
  use sphdata
  use sphneighbourdata
  IMPLICIT NONE
  integer :: i,ic,ipar,j,k,l,n_open
  integer :: n_test,failflag
  real,allocatable, dimension(:) :: hit
  real,allocatable,dimension(:,:) :: bb_cen

  real, dimension(3) :: r_choose
  character(7) ::sphfile
  character(18) :: neighbours
  logical :: existneigh

  real :: hi,hj, hmean	, r_sj, r_sn, sep,percent
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


     allocate(hit(n_node))

     !	Calculate centre of all AABBs

     allocate(bb_cen(3,n_node))	

    !$OMP PARALLEL &
    !$OMP shared(bb_cen, bbr_min, bbr_max,n_node) &
    !$OMP private(ic,k)
    !$OMP DO SCHEDULE(runtime)
     do ic=1,n_node			
        do k=1,3
           bb_cen(k,ic) = bbr_min(k,ic) + (bbr_max(k,ic)-bbr_min(k,ic))/2.0
        enddo
     enddo
    !$OMP END DO
    !$OMP END PARALLEL

     ! Now loop over all particles

     do j=1,nelement

         if(iphase(j)<0) cycle

         ! printing the extent to which the process is complete
         percent = 100.0*REAL(j)/REAL(nelement)

         if(MOD(real(j), real(nelement)/10.0)<1.0) THEN
             write(*,'(I3," % complete")') INT(percent)
         ENDIF

         ! Set all nodes to unhit
         hit(:) = 0

         n_test= 0

         ! The root node is always hit
         hit(1) = 1

         !	Make sure the particle's own leaf node is also hit
         hit(partbin(j)) = 1

         ipar = 0

         n_open = 8

         do while(ipar<n_node.and.n_open>0)

             !	Cycle through nodes until you find one that has been hit
             ipar = ipar+1

             IF(hit(ipar)==0) cycle
             IF(ipar> n_node) exit	!	Exit when ipar is too high

             !	Investigate children of this node

             IF(n_child(ipar)>0) then

                 do i=1,n_child(ipar)
                     ic = child(ipar,i)

                     ! **********************************************
                     !	Perform intersection tests with AABB
                     !	**********************************************

                     !	Is particle inside the AABB? Start by assuming it is

                     failflag = 0 ! If flag remains zero, particle is inside the AABB

                     do k=1,3
                         IF(xyzmh(k,j)>bbr_max(k,ic)) failflag=1
                         IF(xyzmh(k,j)<bbr_min(k,ic)) failflag=1
                     enddo

                     !	IF particle passes this test, proceed to next stage
                     IF(failflag==0) GOTO 10

                     !	If binning test fails, proceed with other tests

                     !	First calculate r_sn = distance from particle to AABB centre

                     r_sn = (xyzmh(1,j) - bb_cen(1,ic))**2 + &
                         (xyzmh(2,j) - bb_cen(2,ic))**2 + &
                         (xyzmh(3,j) - bb_cen(3,ic))**2 + tiny
                     r_sn = SQRT(r_sn)

                     !	If this distance is less than 2h, node is hit

                     IF(r_sn<tolerance*xyzmh(5,j)) THEN
                         failflag = 0
                         GOTO 10
                     ENDIF

                     !	If this distance is greater than 2h, identify nearest vertex of the AABB

                     ! r_choose determines closest face for each axis
                     ! Convergence of all 3 closest faces = closest vertex

                     do k=1,3
                         r_choose(k) = bbr_min(k,ic)

                         IF(ABS(xyzmh(k,j) - bbr_min(k,ic)) > ABS(xyzmh(k,j) - bbr_max(k,ic))) then
                             r_choose(k) = bbr_max(k,ic)
                         endif

                     enddo

                     !	Calculate r_sj = distance from particle to nearest vertex

                     r_sj = (xyzmh(1,j) - r_choose(1))**2 + &
                         (xyzmh(2,j) - r_choose(2))**2 + &
                         (xyzmh(3,j) - r_choose(3))**2 + tiny

                     r_sj = SQRT(r_sj)

                     !	If r_sj < 2h, then node is hit

                     IF(r_sj< tolerance*xyzmh(5,j)) then
                         failflag = 0
                     ENDIF

                 !	If intersect, identify node as hit

10               CONTINUE

                 IF(failflag==1) n_open=n_open-1 ! If node not hit, then close it and continue

                 IF(failflag==0) then ! If node is hit, then record the hit and open all its children for examination
                     hit(ic) = 1
                     n_open = n_open+n_child(ic)
                 endif

             enddo	!	End of loop over children nodes

            ! If current node has no children (i.e. it's a leaf node), then test its occupants for neighbours
         ELSE IF(n_child(ipar) ==0) then

             hj = xyzmh(5,j)

    ! Parallelise here - probably best compromise
    ! Doing full parallelisation doesn't seem to work

    !$OMP PARALLEL &
    !$OMP shared(n_occ,occ,iphase, hj,j,ipar,ic,neighb,nneigh) &
    !$OMP private(l,i,hmean, sep) &
    !$OMP reduction(+:n_test)
    !$OMP DO SCHEDULE(runtime)
             do l=1,n_occ(ipar)

                 i = occ(ipar,l)

                 IF(iphase(i)<0) cycle ! Skip accreted particles, etc

                 n_test = n_test +1

                 hi = xyzmh(5,i)

                 if(i/=j) then
                     hmean = (hi + hj)/2.0

                     sep = (xyzmh(1,i) - xyzmh(1,j))**2 + &
                         (xyzmh(2,i) - xyzmh(2,j))**2 + &
                         (xyzmh(3,i) - xyzmh(3,j))**2 + tiny
                     sep = sqrt(sep)


                     !	if particle j in neighbour sphere, then add to neighbour list
                     if(sep<2.0*hmean.and.nneigh(j)<neighmax) then
                        !$OMP CRITICAL
                         nneigh(j) = nneigh(j) + 1
                         neighb(j,nneigh(j)) = i
                         !$OMP END CRITICAL

                     endif
                 endif

             enddo
             !$OMP END DO
             !$OMP END PARALLEL

            ! End loop over leaf occupants
         endif

         ! Now that the node has been investigated fully, close it
         n_open = n_open-1
     enddo
     ! End loop over nodes

     if(nneigh(j) < 20) then
         print*, 'n_test for particle ',j,' is ',n_test, nneigh(j), neighmax
         print*, 'Number of nodes tested: ', int(sum(hit))
     endif

 enddo

 ! End loop over particles

   
 ! writing to the file.
 print*, 'Writing to file ', neighbours

 call write_neighbours(neighbours,tolerance)

 endif

 ! Calculate mean and standard deviation of neighbour counts
 meanneigh = sum(nneigh)/REAL(nelement)
 sdneigh = 0.0

!$OMP PARALLEL &
!$OMP shared(nneigh,meanneigh,nelement)&
!$OMP private(i) &
!$OMP reduction(+:sdneigh)
!$OMP DO SCHEDULE(runtime)
 do i=1,nelement
     sdneigh = sdneigh+(nneigh(i)-meanneigh)**2
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 sdneigh = sqrt(sdneigh/REAL(nelement))

 print*, 'Mean neighbour number is ', meanneigh
 print*, 'Standard Deviation: ', sdneigh
 neighcrit = meanneigh-5.0*sdneigh
     
 print*, 'Clumps created only if neighbour number greater than ', neighcrit



 return
 end subroutine neighbours_octree
