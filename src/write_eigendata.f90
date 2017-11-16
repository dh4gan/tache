subroutine write_eigendata(ifile)

  !
  ! Subroutine writes all eigendata to binary files
  !

  use tachedata
  use sphdata, only: xyzmh,rho,iunique,iphase

  implicit none

  integer, intent(in) :: ifile

  integer :: i,j,k,counter

  real, allocatable, dimension (:,:,:) :: eigenvecbin
  real, allocatable, dimension(:,:) :: eigenbin
  real, allocatable, dimension(:) :: xbin, ybin, zbin, mbin, rhobin
  

  print*, '----------------------'
  print*, 'Writing eigenvalues to file ', TRIM(eigenfile(ifile))
  print*, 'Writing eigenvectors to file ', TRIM(vectorfile(ifile))

  ! Write data to file - for now, simple formatted file
  
  ! Write data to binary file
  
  ! Prepare special arrays for binary write

  allocate(xbin(neigen))
  allocate(ybin(neigen))
  allocate(zbin(neigen))
  allocate(rhobin(neigen))
  allocate(mbin(neigen))
  allocate(eigenbin(3,neigen))
  allocate(eigenvecbin(3,3,neigen))
  allocate(eigenelement(neigen))

  counter =1
  do ielement=1,nelement
     if(iphase(ielement)/=0) cycle
!     if(allocated(iunique)) then
!        eigenelement(counter) = iunique(ielement)
!     else
        eigenelement(counter) = ielement
!     endif
     xbin(counter) = xyzmh(1,ielement)
     ybin(counter) = xyzmh(2,ielement)
     zbin(counter) = xyzmh(3,ielement)
     rhobin(counter) = rho(ielement)
     mbin(counter) = xyzmh(4,ielement)

     do k=1,3
        eigenbin(k,counter) = eigenvalues(k,ielement)
        do j=1,3
           eigenvecbin(j,k,counter) = eigenvectors(j,k,ielement)
        enddo
     enddo
     
     counter = counter +1
  enddo

  open(27,file=eigenfile(ifile), status='unknown',form='unformatted')
  write(27) neigen
  write(27) (eigenelement(i),i=1,neigen)
  write(27) (xbin(i), i=1,neigen)
  write(27) (ybin(i), i=1,neigen)
  write(27) (zbin(i), i=1,neigen)
  write(27) (eigenbin(1,i), i=1,neigen)
  write(27) (eigenbin(2,i), i=1,neigen)     
  write(27) (eigenbin(3,i), i=1,neigen)
  write(27) (rhobin(i), i =1,neigen)
  write(27) (mbin(i),i=1,neigen)
  close(27)
  
  ! Now write the eigenvectors to file
  open(27,file=vectorfile(ifile),status='unknown', form='unformatted')
  write(27) neigen
  write(27) (eigenelement(i),i=1,neigen)
  write(27) (eigenvecbin(1,1:3,i),i=1,neigen)
  write(27) (eigenvecbin(2,1:3,i),i=1,neigen)
  write(27) (eigenvecbin(3,1:3,i),i=1,neigen)


  deallocate(xbin,ybin,zbin,eigenbin, eigenvecbin)
  if(splitdump.eqv..false.) deallocate(eigenelement) 

end subroutine write_eigendata
