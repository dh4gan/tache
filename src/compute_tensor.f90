subroutine compute_tensor(ifile)
!***********************************************************************
! Subroutine handles the computation of both tensors
!***********************************************************************

use tachedata
implicit none

integer :: ifile
real :: percentcount
real,dimension(3,3) :: tensor_element

print*, '-------------'
! If this is an SPH file, require neighbour lists to compute derivatives
if(filetype=='SPH') then

   call get_SPH_neighbours(ifile)

   allocate(tensor(3,3,nelement))

   if (tensorchoice=='tidal') then
      print*, "Calculating Tidal Tensor"
   else if(tensorchoice=='velocity') then      
      print*, "Calculating Velocity Shear Tensor"
   endif
   percentcount=0.0

   ! Loop over particles
   do ielement = 1,nelement

      call element_percent_complete(ielement,nelement,percentcount,10.0)
      
      tensor_element(:,:) =0.0

      ! Choose which tensor to compute
      if(tensorchoice=='velocity') then
         call calc_velocityshear_tensor(ielement,tensor_element)
      else if(tensorchoice=='tidal') then
         call calc_tidal_tensor(ielement,tensor_element)
      endif

      tensor(:,:,ielement) = tensor_element(:,:)
      
   enddo
   ! End of loop over particles
   
endif
print*, "Tensors calculated"
print*,'-------------------'
   
end subroutine compute_tensor
