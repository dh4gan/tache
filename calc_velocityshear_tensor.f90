SUBROUTINE calc_velocityshear_tensor(ipart, tensor)
! Subroutine calculates tidal tensor for particle ipart

use sphgravdata
use treedata
use sphkerneldata

real, dimension(3,3) :: tensor

integer :: ipart, j,k, imat, jmat

real, dimension(3) :: dr

real :: vmag

do k = 1, nneigh(ipart)

        j = neighb(ipart,k)
        
    ! Calculate gradient of SPH kernel
        hmean = 0.5*(xyzmh(5,ipart) + xyzmh(5,j))
        hmean21 = 1./(hmean*hmean)
        hmean41 = hmean21*hmean21
        pmassj = xyzmh(4,j)
        rhoj = rho(j)

    ! Separation of particles
        do jmat = 1,3
            dr(jmat) = xyzmh(jmat,ipart) - xyzmh(jmat,j)
        enddo

        rij2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

        v2 = rij2*hmean21
        v = rij/hmean

        index = v2/dvtable

        if(index <1) index=1

        dxx = v2 - index*dvtable
        index1 = index + 1

        IF (index1.GT.itable) then
           index = itable-1
           index1 = itable
        endif

        IF(index1 < 0) THEN
           print*, j,ipart, index, itable, rij2, v2,dvtable, v2/dvtable, index1
        endif
        IF (v > radkernel) cycle
     
        dgrwdx = (grwij(index1) - grwij(index))/dvtable
        grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
        grpm = pmassj*grwtij

        ! Loops for matrix elements i,j
        do imat = 1,3
            do jmat = 1,3
                tensor(imat,jmat) = tensor(imat,jmat) - (grpm/rhoj)*(dr(imat)*vxyzu(jmat,j) + dr(jmat)*vxyzu(imat,j))                
            enddo
        enddo
        ! End of matrix loop
enddo

! End of loop over nearest neighbours

! Make tensor dimensionless by scaling by velocity and smoothing length

vmag = sqrt(vxyzu(1,ipart)*vxyzu(1,ipart) + vxyzu(2,ipart)*vxyzu(2,ipart) + &
     vxyzu(3,ipart)*vxyzu(3,ipart))

tensor(:,:) = 0.5*tensor(:,:)*xyzmh(5,ipart)/vmag

!tensor(:,:) = 0.5*tensor(:,:)

return
END SUBROUTINE calc_velocityshear_tensor


