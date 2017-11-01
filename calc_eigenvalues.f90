subroutine calc_eigenvalues 


print*, "Calculating Tensor Eigenvalues"

percentcount = 0.0

allocate(eigenvalues(3,npart))
allocate(eigenvectors(3,3,npart))

if(filetype=='SPH') then

   do ipart = 1,npart

      call particle_percent_complete(ipart,npart,percentcount,10.0)
      
      tensor(:,:) = velocityshear(:,:,ipart)
      eigen(:) = 0.0
      eigenvalues(:,ipart) = 0.0
      
      if(iphase(ipart)/=0) cycle
      
      it_num = 0
      rot_num = 0
      
      call jacobi_eigenvalue(3,tensor, it_max, eigenvec, eigen,it_num, rot_num)
      
      eigenvalues(:,ipart) = eigen(:)
      eigenvectors(:,:,ipart) = eigenvec(:,:)
   enddo

print*, 'Eigenvalues calculated'
   
endif

end subroutine calc_eigenvalues
