SUBROUTINE particle_forces(ipart,jpart, potential, gravforce)
! Calculates contributions to the gravitational potential and force vector
! of ipart due to jpart

use sphgravdata
use sphkerneldata

implicit none

integer :: ipart, jpart, ix

real, dimension(3) :: dr
real, dimension(3), intent(inout):: gravforce
real, intent(inout) :: potential
real :: dph, v

if(ipart==jpart) return
if(iphase(ipart)<0) return
if(iphase(jpart)<0) return


 ! Separation of particles
              do ix = 1,3
                 dr(ix) = xyzmh(ix,ipart) - xyzmh(ix,jpart)
              enddo

              rij2 = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
              rij = sqrt(rij2)
              rij1 = 1./rij
              pmassj = xyzmh(4,jpart)
              hi = xyzmh(5,ipart)
              hj = xyzmh(5,jpart)

              !
              !--Define mean h
              !

              IF (iphase(ipart).GE.1) THEN
                 hmean = hj/2.0
              ELSEIF (iphase(jpart).GE.1) THEN
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

              IF (index1.GT.itable) then
                 index1 = itable            
                 index = index1-1
              endif

              IF (isoft.EQ.1) THEN
                 rij2grav = rij2 + psoft**2
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
                 gravforce(ix) = gravforce(ix) - xmasj*dr(ix)
              ENDDO

              potential = potential + phi*pmassj
              dphiti = dphiti + pmassj*dph


END SUBROUTINE particle_forces
