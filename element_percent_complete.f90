  subroutine element_percent_complete(i, ntot,previouspercent,increment)
    !
    ! Calculates (and displays) % completion at element i of total ntotal
    !

    
    implicit none
    integer :: i,ntot
    real :: percent, previouspercent,increment

    if(previouspercent==0.0) previouspercent = increment
    percent = (real(i)/real(ntot))*100.0
    if(percent > previouspercent) then
       write(*,"(F4.0,A)") previouspercent, "% complete"
       previouspercent = previouspercent + 10.0
    endif

    return
  end subroutine element_percent_complete
