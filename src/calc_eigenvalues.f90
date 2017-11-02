subroutine calc_eigenvalues 
!
! Subroutine calculates the tensor eigenvalues/eigenvectors
!

use tachedata
use sphdata, only:iphase

implicit none

integer, parameter :: it_max = 10
integer :: it_num,rot_num

real :: percentcount

real, dimension(3,3) :: tensor_element,eigenvec
real, dimension(3) :: eigen

print*, "Calculating Tensor Eigenvalues"

percentcount = 0.0

allocate(eigenvalues(3,nelement))
allocate(eigenvectors(3,3,nelement))

if(filetype=='SPH') then

   do ielement = 1,nelement

      call element_percent_complete(ielement,nelement,percentcount,10.0)
      
      tensor_element(:,:) = tensor(:,:,ielement)
      eigen(:) = 0.0
      eigenvalues(:,ielement) = 0.0
      
      if(iphase(ielement)/=0) cycle
      
      it_num = 0
      rot_num = 0
      
      call jacobi_eigenvalue(3,tensor_element, it_max, eigenvec, eigen,it_num, rot_num)
      
      eigenvalues(:,ielement) = eigen(:)
      eigenvectors(:,:,ielement) = eigenvec(:,:)
   enddo

print*, 'Eigenvalues calculated'
   
endif

end subroutine calc_eigenvalues
