!------------------------------------------------------
!+
! Does brute force SPH density sum (assuming current h)
!+
!-------------------------------------------------------
subroutine recalc_density

use sphdata
use sphkerneldata
use tachedata,only: nelement

implicit none

integer :: ielement,jelement,k,nneigh
real,dimension(3) :: dr

print*, 'Recalculating density assuming current smoothing lengths'

rho(:) = 0.0
do ielement=1,nelement
   !print*, ielement

   nneigh = 0
   if(iphase(ielement)/=0) cycle
   
   do jelement=ielement,nelement 
      if(jelement==ielement) cycle
      
      ! Calculate gradient of SPH kernel
      hmean = 0.5*(xyzmh(5,ielement) + xyzmh(5,jelement))
      hmean21 = 1./(hmean*hmean)
      hmean31 = hmean21/hmean
      
      ! Separation of particles
      do k = 1,3
         dr(k) = xyzmh(k,ielement) - xyzmh(k,jelement)
      enddo
      
      rij2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
     
      v2 = rij2*hmean21
     
      
 
      if(v2>radkernel*radkernel) cycle

      nneigh = nneigh+1
      index = v2/dvtable
      
      if(index<1) index=1
      
      dxx = v2 - index*dvtable
      index1 = index + 1
      
      IF (index1.GT.itable) then
         index = itable-1
         index1 = itable
      endif
      
      dwdx = (wij(index1)-wij(index))/dvtable
      wtij = (wij(index) + dwdx*dxx)*hmean31
      rho(ielement) = rho(ielement) + xyzmh(4,jelement)*wtij
      rho(jelement) = rho(jelement) + xyzmh(4,ielement)*wtij

   enddo


   !print*, ielement, ' New density: ', rho(ielement), nneigh, rho(ielement)/rho_old
enddo
rho(:) = rho(:)*cnormk

end subroutine recalc_density
