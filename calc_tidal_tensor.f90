SUBROUTINE calc_tidal_tensor(ielement, tidal)
! Subroutine calculates tidal tensor for particle ielement
! Tensor is scaled by its potential

use sphgravdata
use treedata
use sphkerneldata

real, dimension(3,3) :: tidal

integer :: ielement, j,k, imat, jmat

real, dimension(3) :: dr

do k = 1, nneigh(ielement)

        j = neighb(ielement,k)
        
    ! Calculate gradient of SPH kernel
        hmean = 0.5*(xyzmh(5,ielement) + xyzmh(5,j))
        hmean21 = 1./(hmean*hmean)
        hmean41 = hmean21*hmean21
        pmassj = xyzmh(4,j)
        rhoj = rho(j)

    ! Separation of particles
        do jmat = 1,3
            dr(jmat) = xyzmh(jmat,ielement) - xyzmh(jmat,j)
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
           print*, j,ielement, index, itable, rij2, v2,dvtable, v2/dvtable, index1
        endif
        IF (v > radkernel) cycle
     
        dgrwdx = (grwij(index1) - grwij(index))/dvtable
        grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
        grpm = pmassj*grwtij

        ! Loops for matrix elements i,j
        do imat = 1,3
            do jmat = 1,3
                tidal(imat,jmat) = tidal(imat,jmat) - grpm*dr(imat)*gravxyz(jmat,j)/rhoj
            enddo
        enddo
        ! End of matrix loop
enddo
! End of loop over nearest neighbours

! Tensor should be dimensionless - scale using potential and smoothing length
tidal(:,:) = tidal(:,:)*xyzmh(5,ielement)*xyzmh(5,ielement)/abs(poten(ielement))


return
END SUBROUTINE calc_tidal_tensor


