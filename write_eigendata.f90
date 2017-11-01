subroutine write_eigendata

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
  allocate(eigenpart(ngas))

  counter =1
  do ipart=1,npart
     if(iphase(ipart)/=0) cycle
     if(allocated(iunique)) then
        eigenpart(counter) = iunique(ipart)
     else
        eigenpart(counter) = ipart
     endif
     xbin(counter) = xyzmh(1,ipart)
     ybin(counter) = xyzmh(2,ipart)
     zbin(counter) = xyzmh(3,ipart)
     do k=1,3
        eigenbin(k,counter) = eigenvalues(k,ipart)
        do j=1,3
           eigenvecbin(j,k,counter) = eigenvectors(j,k,ipart)
        enddo
     enddo
     
     counter = counter +1
  enddo
  
  open(27,file=eigenfile(n), status='unknown',form='unformatted')
  write(27) ngas
  write(27) (eigenpart(i),i=1,ngas)
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
  write(27) (eigenpart(i),i=1,ngas)
  write(27) (eigenvecbin(1,1:3,i),i=1,ngas)
  write(27) (eigenvecbin(2,1:3,i),i=1,ngas)
  write(27) (eigenvecbin(3,1:3,i),i=1,ngas)


  deallocate(xbin,ybin,zbin,eigenbin,eigenpart, eigenvecbin)



end subroutine write_eigendata
