SUBROUTINE calc_grav_tree
  ! Subroutine calculates gravitational forces and potentials 
  ! for all particles
  ! It uses a pre-generated octree to approximate gravity at distances 
  ! defined by Multipole Acceptance Criterion (MAC)

    use sphdata
    use sphneighbourdata
    use sphkerneldata
    use tachedata, only: nelement

    implicit none

    real, dimension(3) :: gravi,drcell, drgeom
    real, parameter :: mactheta = 0.01

    integer, allocatable, dimension(:) :: used

    integer :: ielement, jpart,jcell,k,ix
    real :: meanpot,sdpot, percent,counter
    real :: sepcell,sepcell2,theta, dimcell, sepgeom

    logical :: fullcalc

    allocate(gravxyz(3,nelement))
    allocate(poten(nelement))

    gravxyz(:,:) = 0.0
    poten(:) = 0.0

    isoft = 0

    percent = 0.0
    counter = 1.0

    allocate(used(n_node))

    used(:) = 0

    ! Begin loop over particles
    do ielement = 1,nelement

       if(iphase(ielement)<0) cycle
        percent = REAL(ielement)/REAL(nelement)*100.0

        if(percent>counter)then
            print*, counter,'% complete'
            counter = counter +1.0
        endif

        gravi(:) = 0.0
        poteni = 0.0
        used(:) = 0
    
        ! Now loop over nodes of the tree
!$OMP PARALLEL &
!$OMP shared(ielement,n_node,m_node,dr_node,com_node,r_node, used,xyzmh)&
!$OMP shared(n_occ, occ,n_child, poten, gravxyz) &
!$OMP private(jcell,jpart,fullcalc,dimcell,drcell,drgeom,sepcell,sepcell2) &
!$OMP private(sepgeom,theta,poteni,gravi) 
!$OMP DO SCHEDULE(runtime)
        do jcell = 1,n_node
           
            ! If this cell has already been used as part of a calculation, skip it
            if(used(jcell)==1) cycle

            if(m_node(jcell)==0) cycle

            ! Begin by assuming full calculation not necessary for this node

            fullcalc = .false.

            ! Calculate distance between ielement and the centre of the cell (R).  If

            ! a) R < 3.0*ielement's smoothing length, or
            ! b) theta = R/D < mactheta (where D is maximum cell dimension)
            ! then do full calculation

            ! Maximum dimension of cell
            dimcell = maxval(dr_node(:,jcell))

            ! separation of particle from cell centre of mass
            do ix=1,3
                drcell(ix) = xyzmh(ix,ielement)-com_node(ix,jcell) ! Separation from cell COM                
                drgeom(ix) = xyzmh(ix,ielement) - r_node(ix,jcell) ! Separation from node geometric centre
            enddo

            sepcell2 = drcell(1)*drcell(1) + drcell(2)*drcell(2) + drcell(3)*drcell(3)
            sepcell = sqrt(sepcell2)

            sepgeom = sqrt(drgeom(1)*drgeom(1) + drgeom(2)*drgeom(2) + drgeom(3)*drgeom(3))

            ! Ratio of cell size to separation
            !if(sepcell-sepgeom>0.0) then
            !   theta = 2.0*dimcell/(sepcell-sepgeom)
            !else
            !   theta = 2.0*dimcell/amin1(sepcell,sepgeom)
            !endif

            theta = 2.0*dimcell/sepcell

            if(theta > mactheta) fullcalc = .true.
            if(sepcell < 3.0*xyzmh(5,ielement)) fullcalc = .true.

            ! If this cell contains particle ielement, full calculation must be done
            ! Also, if cell only contains one particle, full calculation must be done
            if(partbin(ielement)==jcell) fullcalc = .true.
            if(n_occ(jcell)==1) fullcalc = .true.
            if(n_child(jcell)==0) fullcalc = .true.
          
            ! If full calculation required (and not a leaf), then continue on to the next node
            ! Otherwise, calculate the forces from the particles in the leaf

            poteni = 0.0
            gravi(:) = 0.0

            if(fullcalc) then

                if(n_child(jcell)==0) then
                                  !print*, 'Particle ', ielement, ' cell, ',jcell, ' using full calc', m_node(jcell), sepcell, &
                                  !     dimcell, sepgeom, theta
                    do k=1,n_occ(jcell)
                        jpart = occ(jcell,k)
                        call particle_forces(ielement,jpart, poteni,gravi)
                    enddo
                endif

            else
                ! If not a full calculation, then use the node's properties

               if(m_node(jcell)>0.0) then

                 ! print*, 'Particle ', ielement, ' cell ', jcell, ' Using Tree ', m_node(jcell), sepcell, & 
                      ! dimcell, sepgeom, theta
                  do ix=1,3
                     gravi(ix) = gravi(ix) - m_node(jcell)*drcell(ix)/(sepcell2*sepcell)
                  enddo

                  poteni = poteni -m_node(jcell)/sepcell
                 
                  ! Mark all the children of this node as "used"
                
                  call mark_children(jcell,used)
                 
               endif
            endif

            ! Add calculated contributions to the potential and gravitational forces
            poten(ielement) = poten(ielement) + poteni

            DO ix=1,3
                gravxyz(ix,ielement) = gravxyz(ix,ielement) + gravi(ix)
            ENDDO            
          
         enddo
!$OMP END DO
!$OMP END PARALLEL

        ! End loop over nodes

         !print*, 'Particle: ',ielement, poten(ielement), gravxyz(:,ielement)
        !  Add contribution to the potential from pointmasses

        do k=1,nptmass

            jpart = listpm(k)

            if(jpart==ielement) cycle

            sep = (xyzmh(1,ielement) - xyzmh(1,jpart))**2 + &
                (xyzmh(2,ielement) - xyzmh(2,jpart))**2 +&
                (xyzmh(3,ielement) - xyzmh(3,jpart))**2

            sep =sqrt(sep)

            poten(ielement) = poten(ielement) - xyzmh(4,jpart)/sep

        enddo
!        print*, 'Particle ', ielement, ': Potential - ', poten(ielement), 'Force: ',gravxyz(:,ielement)
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
END SUBROUTINE calc_grav_tree


