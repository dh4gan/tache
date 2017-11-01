subroutine compute_tensor
!
! Subroutine handles the computation of derivatives that form the tensor
! to be classified
!

real, dimension(3,3) :: tensor_element
real, allocatable, dimension(:,:,:) :: tensor

! If this is an SPH file, require neighbour lists to compute derivatives
if(filetype=='SPH') then

   call get_SPH_neighbours

   allocate(tensor(3,3,nelement))

   if (tensorchoice=='tidal') then
      print*, "Calculating Tidal Tensor"
   else if(tensorchoice=='velocity') then      
      print*, "Calculating Velocity Shear Tensor"
   endif
   percentcount=0.0

   ! Loop over particles
   do ielement = 1,nelement

      call particle_percent_complete(ielement,nelement,percentcount,10.0)
      
      tensor_element(:,:) =0.0
      call calc_velocityshear_tensor(ielement,tensor_element)
      tensor(:,:,ielement) = tensor_element(:,:)
      
   enddo
   ! End of loop over particles
   
endif
print*, "Tensors calculated"
print*,'-------------------'
   
end subroutine compute_tensor
