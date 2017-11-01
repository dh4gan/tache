subroutine write_eigendata

  integer,allocatable,dimension(:) :: eigenelement
  real, allocatable, dimension (:,:,:) :: eigenvecbin
  real, allocatable, dimension(:,:) :: eigenvalues, eigenbin
  real, allocatable, dimension(:) :: xbin, ybin, zbin
  

  print*, '----------------------'
  print*, 'Writing to file ', TRIM(eigenfile(n))

  ! Write data to file - for now, simple formatted file
  
  ! Write data to binary file
  
  ! Prepare special arrays for binary write

  allocate(xbin(ngas))
  allocate(ybin(ngas))
  allocate(zbin(ngas))
  allocate(eigenbin(3,ngas))
  allocate(eigenvecbin(3,3,ngas))
  allocate(eigenelement(ngas))

  counter =1
  do ielement=1,nelement
     if(iphase(ielement)/=0) cycle
     if(allocated(iunique)) then
        eigenelement(counter) = iunique(ielement)
     else
        eigenelement(counter) = ielement
     endif
     xbin(counter) = xyzmh(1,ielement)
     ybin(counter) = xyzmh(2,ielement)
     zbin(counter) = xyzmh(3,ielement)
     do k=1,3
        eigenbin(k,counter) = eigenvalues(k,ielement)
        do j=1,3
           eigenvecbin(j,k,counter) = eigenvectors(j,k,ielement)
        enddo
     enddo
     
     counter = counter +1
  enddo
  
  open(27,file=eigenfile(n), status='unknown',form='unformatted')
  write(27) ngas
  write(27) (eigenelement(i),i=1,ngas)
  write(27) (xbin(i), i=1,ngas)
  write(27) (ybin(i), i=1,ngas)
  write(27) (zbin(i), i=1,ngas)
  write(27) (eigenbin(1,i), i=1,ngas)
  write(27) (eigenbin(2,i), i=1,ngas)     
  write(27) (eigenbin(3,i), i=1,ngas)
  close(27)
  
  ! Now write the eigenvectors to file
  open(27,file=vectorfile(n),status='unknown', form='unformatted')
  write(27) ngas
  write(27) (eigenelement(i),i=1,ngas)
  write(27) (eigenvecbin(1,1:3,i),i=1,ngas)
  write(27) (eigenvecbin(2,1:3,i),i=1,ngas)
  write(27) (eigenvecbin(3,1:3,i),i=1,ngas)


  deallocate(xbin,ybin,zbin,eigenbin,eigenelement, eigenvecbin)



end subroutine write_eigendata
