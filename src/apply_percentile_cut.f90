subroutine apply_percentile_cut
!
! Removes all but the top x% in density from the analysis
!

implicit none 

  ! Only consider the top x percent

   print '(a,f5.1,a)', 'Only considering the top ',xpercentile,' density percentile'

   ipercentile = int(xpercentile*real(nelement)/100.0)

   print'(a,I7,a)', 'This constitutes ', ipercentile, ' elements'
   
   do ielement= ipercentile,nelement
      i = isort(nelement-ielement+1)
      spiralmember(i)=-1
   enddo




end subroutine apply_percentile_cut
