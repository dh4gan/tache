      subroutine ktable
!------------------------------------------------------------
!                                                           -
!  This subroutine builds a table for the various values    -
!     of the kernel, the gradient of the kernel, the mass   -
!     fraction and the potential energy.                    -
!     The entry is v**2.                                    -
!     (updated to F90 DHF 17/10/2009                        -
!------------------------------------------------------------

	 use sphkerneldata

	 real :: sum, v, v3, v4, v5, v6
	 real, parameter :: pi = 3.141596253d0
	 integer :: i1
	 
	 
!	Allocate arrays

	allocate(wij(itable))
	allocate(grwij(itable))
	allocate(fmass(itable))
	allocate(fpoten(itable))
	allocate(dphidh(itable))
	 
!
!--Maximum interaction length and step size
!
	  radkernel = 2.0
      part1kernel = 1.0
      part2kernel = 2.0
      v2max = radkernel*radkernel
      dvtable = v2max/itable
      i1 = part1kernel/dvtable +1
!
!--Build tables
!
!  a) v less than 1
!
      do i = 1, i1
	  
         v2 = i*dvtable
		 v = SQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
		 
         sum = 1. - 1.5*v2 + 0.75*v3
         wij(i) = sum
         sum = -3.*v + 2.25*v2
         grwij(i) = sum
         sum = 1.3333333333*v3 - 1.2*v5 + 0.5*v6
         fmass(i) = sum
         sum = 0.66666666666*v2 - 0.3*v4 + 0.1*v5 - 1.4
         fpoten(i) = sum
         sum = -1.4 + 2.*v2 - 1.5*v4 + 0.6*v5
         dphidh(i) = sum
      enddo
	  
	
!
!  b) v greater than 1
!
      do i = i1 + 1, itable
         v2 = i*dvtable
         v = SQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         dif2 = 2. - v
         sum = 0.25*dif2*dif2*dif2
         wij(i) = sum
         sum = -0.75*v2 + 3*v - 3.
         grwij(i) = sum
         sum = -0.16666666666*v6 + 1.2*v5 - 3.*v4 + 2.66666666666*v3 - &
              0.0666666666666
         fmass(i) = sum
         sum = -0.033333333333*v5 + 0.3*v4 - v3 + 1.3333333333*v2 - 1.6
         fpoten(i) = sum
         sum = -1.6 + 4.*v2 - 4.*v3 + 1.5*v4 - 0.2*v5
         dphidh(i) = sum
      enddo
	  
!
!--Normalisation constant
!
      cnormk = 1.0/pi
      selfnormkernel = 1.0
      part1potenkernel = 1.0/15.0
      part2potenkernel = 0.0
	
      return
	  
      end subroutine ktable
