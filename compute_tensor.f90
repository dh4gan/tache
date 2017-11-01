subroutine compute_tensor
!
! Subroutine handles the computation of derivatives that form the tensor
! to be classified
!


! If this is an SPH file, require neighbour lists to compute derivatives
if(filetype=='SPH') then

   call get_SPH_neighbours

   allocate(tensor(3,3,npart))

   if (tensorchoice=='tidal') then
      print*, "Calculating Tidal Tensor"
   else if(tensorchoice=='velocity') then      
      print*, "Calculating Velocity Shear Tensor"
   endif
   percentcount=0.0

   ! Loop over particles
   do ipart = 1,npart

      call particle_percent_complete(ipart,npart,percentcount,10.0)
      
      tensor_element(:,:) =0.0
      call calc_velocityshear_tensor(ipart,tensor_element)
      tensor(:,:,ipart) = tensor_element(:,:)
      
   enddo
   ! End of loop over particles
   
endif
print*, "Tensors calculated"
print*,'-------------------'
   
end subroutine compute_tensor
