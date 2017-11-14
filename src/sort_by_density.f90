subroutine sort_by_density(nelement,isort)
  !
  ! subroutine sorts elements in order of increasing density 
  ! and performs a density
  !

  implicit none

  print*, 'Sorting elements by density'

  allocate(isort(nelement))
  allocate(tested(nelement))

  tested(:) = 0
     
   do ielement=1,nelement
      isort(ielement) = ielement
   enddo
     
   CALL sort2(npart,rho,isort,npart)

   print*, 'Particles sorted by Density'
   print*, "-----------------------------------------------"
   
end subroutine sort_by_density
