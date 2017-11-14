subroutine shift_to_centre_of_mass(xcom,vcom)

real,dimension(3) :: xcom,vcom

xcom(:) = 0.0
vcom(:) = 0.0
totmass = 0.0

do ipart=1,npart
   if(iphase(ipart)/=0) cycle

   do k=1,3
      xcom(k) = xcom(k) + xyzmh(k,ipart)*xyzmh(4,ipart)
      vcom(k) = vcom(k) + vxyzu(k,ipart)*xyzmh(4,ipart)
      totmass = totmass + xyzmh(4,ipart)
   enddo
   
enddo

if(totmass>1.0e-30) then
   xcom(:) = xcom(:)/totmass
   vcom(:) = vcom(:)/totmass
else
   xcom(:) = 0.0
   vcom(:) = 0.0
endif

print*, 'Centre of Mass: Position', xcom
print*, 'Centre of Mass: Velocity', vcom

xyzmh(1:3,ipart) = xyzmh(1:3,ipart) - xcom(1:3)
vxyzu(1:3,ipart) = vxyzu(1:3,ipart) - vcom(1:3)



end subroutine shift_to_centre_of_mass
