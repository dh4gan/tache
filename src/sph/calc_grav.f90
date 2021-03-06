SUBROUTINE calc_grav
!
! Subroutine calculates gravitational forces and potentials for all elements
!

use sphdata
use sphneighbourdata
use sphkerneldata
use tachedata, only: nelement

implicit none

real,dimension(3) :: gravi,dr

integer :: ielement, j,k,ix
real :: v, meanpot,sdpot

gravxyz(:,:) = 0.0
poten(:) = 0.0

do ielement = 1,nelement

gravi(:) = 0.0
poteni = 0.0

do k = 1, nneigh(ielement)

   j = neighb(ielement,k)

   if(j==0) print*, 'AARGH, j=0', j, k,ielement, neighb(ielement,k)

   ! Separation of particles
   do ix = 1,3
      dr(ix) = xyzmh(ix,ielement) - xyzmh(ix,j)
   enddo

   rij2 = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
   rij = sqrt(rij2)
   rij1 = 1./rij
   pmassj = xyzmh(4,j)
   hi = xyzmh(5,ielement)
   hj = xyzmh(5,j)

   !
   !--Define mean h
   !
   IF (iphase(ielement).GE.1) THEN
      hmean = hj/2.0
   ELSEIF (iphase(j).GE.1) THEN
      hmean = hi/2.0
   ELSE
      hmean = 0.5*(hi + hj)
   ENDIF
   hmean21 = 1./(hmean*hmean)
   hmean41 = hmean21*hmean21
   
   v2 = rij2*hmean21
   v = rij/hmean
   
   index = v2/dvtable
   dxx = v2 - index*dvtable
   
   if(index<1) index=1
   index1 = index + 1
   
   IF (index1.GT.itable) index1 = itable
   
   
   if(index<1) print*, 'INDEX ERROR: ',index,v2/dvtable,v2,dvtable,hmean, rij2,hmean21
   IF (isoft.EQ.1) THEN
      rij2grav = dx*dx + dy*dy + dz*dz + psoft**2
      rijgrav = SQRT(rij2grav)
      fm = 1.0
      phi = - 1./rijgrav
      dphi = 0.0
   ELSEIF (isoft.EQ.0) THEN
      rij2grav = rij2
      rijgrav = rij
      IF (v.GE.radkernel) THEN
         fm = 1.0
         phi = -rij1
         dphi = 0.0
      ELSE
         dfmassdx = (fmass(index1) - fmass(index))/dvtable
         fm = (fmass(index) + dfmassdx*dxx)
         dfptdx = (fpoten(index1) - fpoten(index))/dvtable
         phi = (fpoten(index) + dfptdx*dxx)/hmean
         dpotdh = (dphidh(index1) - dphidh(index))/dvtable
         dphi = (dphidh(index) + dpotdh*dxx)*hmean21*dhmean
         IF (v.GT.part2kernel) THEN
            phi = phi + rij1*part2potenkernel
         ELSEIF (v.GT.part1kernel) THEN
            phi = phi + rij1*part1potenkernel
         ENDIF
      ENDIF
      
   ENDIF
   !
   !--Gravitational force calculation
   !
   xmasj = fm*pmassj/(rij2grav*rijgrav)
   
   DO ix=1,3
      gravi(ix) = gravi(ix) - xmasj*dr(ix)
   ENDDO
   
   poteni = poteni + phi*pmassj
   dphiti = dphiti + pmassj*dphi
   
   
   
enddo
! End of loop over nearest neighbours

poten(ielement) = poten(ielement) + poteni

DO ix=1,3
   gravxyz(ix,ielement) = gravxyz(ix,ielement) + gravi(ix)
ENDDO

!  Add contribution to the potential from pointmasses

do k=1,nptmass
   
   j = listpm(k)
   
   if(j==ielement) cycle
   
   sep = (xyzmh(1,ielement) - xyzmh(1,j))**2 + &
        (xyzmh(2,ielement) - xyzmh(2,j))**2 +&
        (xyzmh(3,ielement) - xyzmh(3,j))**2
   
   sep =sqrt(sep)
   
   poten(ielement) = poten(ielement) + xyzmh(4,j)/sep
   
enddo

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
END SUBROUTINE calc_grav


