SUBROUTINE make_octree(sphfile)
  !******************************************************
  ! Builds an octree from the SPH dataset
  ! (used for neighbour finding and gravity calculations)
  !*******************************************************

  USE sphdata	
  USE sphneighbourdata
  use tachedata,only: nelement

  IMPLICIT NONE  

  integer :: ipar,n_open,binlog,ic,ielement,j,l
  integer :: flag,testnode,maxcheck, n_leaf
  integer,allocatable, dimension(:) :: ibin,counter, is_leaf
  character(7) :: sphfile, filecheck,filecheck2
  character(14):: nodefile
character(18)::partbinfile

  logical :: existnode,partbinexist

  !	Check if a nodes file currently exists
  print*, "-----------------------------------------------"
  print*, 'Searching for node fileset'
  write(nodefile,'("octree_",A7)') sphfile
  write(partbinfile,'("octree_bin",A7)') sphfile	

  INQUIRE(file=nodefile,exist=existnode)
  INQUIRE(file=partbinfile,exist=partbinexist)

  existnode = existnode.and.partbinexist	

  IF(existnode.eqv..true.) then
      print*, 'nodes fileset found: reading'
      open(2,file=nodefile,form='unformatted')

      read(2) filecheck
      read(2) n_node,maxcheck


      !	If nodemax incorrect, then file cannot be used
      IF(nodemax/=maxcheck) then
          print*, 'nodemax mismatch: node file invalid'
          print*, 'Require nodemax of ', nodemax, ', find nodemax of ',maxcheck
          print*, 'Generating new octree'
          existnode = .false.
          GOTO 100
      else
          print*, 'node file valid'

      endif

      !	allocate tree variables

      allocate(parent(n_node))
      allocate(child(n_node,8))
      allocate(n_child(n_node))
      allocate(n_occ(n_node))
      allocate(m_node(n_node))
      allocate(occ(n_node,nodemax))
      allocate(ibin(n_node))
      allocate(partbin(nelement))

      allocate(r_node(3,n_node))
      allocate(com_node(3,n_node))
      allocate(dr_node(3,n_node))
      allocate(bbr_min(3,n_node))
      allocate(bbr_max(3,n_node))


      do ic=1,n_node
          read(2) l, parent(ic),n_child(ic), &
              r_node(1,ic),r_node(2,ic),r_node(3,ic), &
              dr_node(1,ic),dr_node(2,ic),dr_node(3,ic), &
              bbr_min(1,ic), bbr_max(1,ic),bbr_min(2,ic), &
              bbr_max(2,ic),bbr_min(3,ic), bbr_max(3,ic), &
              n_occ(ic),m_node(ic), &
              com_node(1,ic), com_node(2,ic), com_node(3,ic)
              
      enddo
      close(2)

      print*, 'Nodes file read'

      !	Must now create children catalogue

      print*, 'Creating child catalogues'

      allocate(counter(n_node))

      counter = 0

      do ic=2,n_node
          IF(parent(ic)/=0) THEN
              counter(parent(ic)) = counter(parent(ic))+1
              child(parent(ic),counter(parent(ic))) = ic
          ENDIF
      enddo
      print*, 'Children catalogued'

      deallocate(counter)

      open(2,file=partbinfile,form='unformatted')
      read(2) filecheck2
      do ielement = 1, nelement
          read(2) l, iphase(ielement), partbin(ielement)
      enddo
      close(2)

      print*, 'Partbin file read'

      print*, 'Updating occupancy catalogue'
      print*, "-----------------------------------------------"
      ibin(:) = 0

      do ielement = 1,nelement
          ic = partbin(ielement)
          ibin(ic) = ibin(ic) + 1
          if(ibin(ic)<=nodemax) occ(ic,ibin(ic)) = ielement
      enddo

  endif

100 CONTINUE


  IF(existnode.eqv..false.) then
     !	Estimate required number of nodes

     print*, 'No satisfactory node fileset found: creating octree from scratch'
     print*, "-----------------------------------------------"
     n_node = INT(nelement/nodemax)	
     testnode =1

     do while (n_node>testnode)
	testnode = testnode*8
     enddo

     n_node = testnode		

999 continue ! Returns here if n_node insufficient
	
     WRITE(*,*) 'Nodes, nodemax, particles:'
     WRITE(*,*) n_node, nodemax, nelement		
     n_leaf = 0

     !	Now create tree - allocate variables

     allocate(parent(n_node))
     allocate(child(n_node,8))		
     allocate(n_child(n_node))
     allocate(n_occ(n_node))
     allocate(m_node(n_node))
     allocate(com_node(3,n_node))
     allocate(occ(n_node,nodemax))
     allocate(partbin(nelement))
     allocate(ibin(n_node))				
     allocate(r_node(3,n_node))
     allocate(dr_node(3,n_node))
     allocate(bbr_min(3,n_node))
     allocate(bbr_max(3,n_node))

     allocate(is_leaf(n_node))

     parent(:) = 0
     child(:,:) = 0
     n_child(:) = 0
     n_occ(:) = 0
     occ(:,:) = 0
     ibin(:) = 0
     r_node(:,:) = 0.0
     com_node(:,:) = 0.0
     dr_node(:,:) = 0.0
     bbr_min(:,:) = 0.0
     bbr_max(:,:) = 0.0
     n_child(:) = 0
     is_leaf(:) = 0
     m_node(:) = 0.0

     !	Algorithm as follows:


     !	do WHILE (nodes still available)
     !	i) Define current parent node (set ipar)
     !	do Loop over nodes
     !	ii) create new children for this node
     !	iii) Bin particles into these new nodes
     !	end do Loop over nodes
     !	iv) Find child nodes that are open
     !	v) Rebin closed nodes to find occupants
     !	end do WHILE

     !	Define root node


     parent(1) = 0
     n_child(1) = 8
     n_occ(1) = nelement

     partbin(:) = 1

     do j=1,3
	r_node(j,1) = 0.0d0			
     enddo

     dr_node(1,1) = xmax
     dr_node(2,1) = ymax
     dr_node(3,1) = zmax

     bbr_min(1,1) = -xmax
     bbr_min(2,1) = -ymax
     bbr_min(3,1) = -zmax

     bbr_max(1,1) = xmax
     bbr_max(2,1) = ymax
     bbr_max(3,1) = zmax

     ipar = 0
     n_open = 1
     inode = 1

     !	do WHILE(open nodes are still available)
     !	(Cannot exceed the maximum number of nodes also)

     do WHILE(inode<n_node.and.ipar<n_node.and.n_open>0)

        !	i) Define Current Parent Node (set ipar)

	ipar = ipar+1

	if(n_child(ipar)==0) cycle

 !	Opening this node ==> remove one from n_open
	n_open = n_open-1

 !	If number of nodes insufficient, then start again			
	if(inode+8>n_node.and.n_open/=0) then
           WRITE(*,*) 'n_nodes insufficient: recalculating'
           n_node = n_node*8

           deallocate(parent,child,n_child)
           deallocate(n_occ, occ, m_node)
           deallocate(ibin, partbin)
           deallocate(r_node,dr_node, com_node)
           deallocate(bbr_min,bbr_max)
           deallocate(is_leaf)



           GOTO 999
 	endif



  !	Open the child nodes
	n_open = n_open+8


	binlog = 0
 !	do Loop over nodes
	do ic=inode+1,inode+8

    !	ii) Create children for this node

    !	Define lineage
           parent(ic) = ipar
           child(ipar,ic-inode) = ic

           !	Define node size (constant over given generation)

           do j=1,3
              dr_node(j,ic) = dr_node(j,ipar)/2.0d0
           enddo
 	enddo

  !	Define each octant using different permutation

  !	Octant 1
	r_node(1,inode+1) = r_node(1,ipar) + dr_node(1,inode+1)
	r_node(2,inode+1) = r_node(2,ipar) + dr_node(2,inode+1)
	r_node(3,inode+1) = r_node(3,ipar) + dr_node(3,inode+1)

 !	Octant 2
	r_node(1,inode+2) = r_node(1,ipar) - dr_node(1,inode+2)
	r_node(2,inode+2) = r_node(2,ipar) + dr_node(2,inode+2)
	r_node(3,inode+2) = r_node(3,ipar) + dr_node(3,inode+2)

 !	Octant 3
	r_node(1,inode+3) = r_node(1,ipar) + dr_node(1,inode+3)
	r_node(2,inode+3) = r_node(2,ipar) - dr_node(2,inode+3)
	r_node(3,inode+3) = r_node(3,ipar) + dr_node(3,inode+3)

 !	Octant 4
	r_node(1,inode+4) = r_node(1,ipar) + dr_node(1,inode+4)
	r_node(2,inode+4) = r_node(2,ipar) + dr_node(2,inode+4)
	r_node(3,inode+4) = r_node(3,ipar) - dr_node(3,inode+4)

 !	Octant 5
	r_node(1,inode+5) = r_node(1,ipar) - dr_node(1,inode+5)
	r_node(2,inode+5) = r_node(2,ipar) - dr_node(2,inode+5)
	r_node(3,inode+5) = r_node(3,ipar) + dr_node(3,inode+5)

 !	Octant 6
	r_node(1,inode+6) = r_node(1,ipar) + dr_node(1,inode+6)
	r_node(2,inode+6) = r_node(2,ipar) - dr_node(2,inode+6)
	r_node(3,inode+6) = r_node(3,ipar) - dr_node(3,inode+6)

 !	Octant 7
	r_node(1,inode+7) = r_node(1,ipar) - dr_node(1,inode+7)
	r_node(2,inode+7) = r_node(2,ipar) + dr_node(2,inode+7)
	r_node(3,inode+7) = r_node(3,ipar) - dr_node(3,inode+7)

 !	Octant 8
	r_node(1,inode+8) = r_node(1,ipar) - dr_node(1,inode+8)
	r_node(2,inode+8) = r_node(2,ipar) - dr_node(2,inode+8)
	r_node(3,inode+8) = r_node(3,ipar) - dr_node(3,inode+8)


 !	iii) Bin particles into these new nodes

	do ic = inode+1, inode+8
           n_occ(ic) = 0
           !	Define minimum AABB (essentially a zero box)

           do j=1,3
              bbr_min(j,ic) = r_node(j,ic)-0.01d0*dr_node(j,ic)
              bbr_max(j,ic) =r_node(j,ic)+0.01d0*dr_node(j,ic)
           enddo

           !	Now begin binning: Record first (nodemax) particles binned, count the rest

           ibin(ic) = 0


           !$OMP PARALLEL &
           !$OMP shared(nelement, iphase, partbin, n_occ, xyzmh) &
           !$OMP shared(ic, ipar, occ, binlog, ibin) &
           !$OMP shared(bbr_min, bbr_max) &
           !$OMP private(ielement,j, flag)
           !$OMP DO SCHEDULE(runtime)

           do ielement = 1,nelement
               if(iphase(ielement) < 0) cycle ! Skip all accreted particles
               flag=0

               !	If particle not in parent node, then no need to test
               if(partbin(ielement)/=ipar) cycle

               !	If all particles in parent node binned, no need to test
               ! any further
               if(binlog==n_occ(ipar)) cycle

               !	Flag particle if outside cell bounds
               do j=1,3
                   if(ABS(xyzmh(j,ielement)-r_node(j,ic))>dr_node(j,ic)) flag=1
               enddo

               !	If particle in cell bounds, then increase n_occ
               if(flag==0) then

                  !$OMP CRITICAL

                   ibin(ic) = ibin(ic) +1
                   binlog = binlog +1
                   if(ibin(ic)<=nodemax) then
                       occ(ic,ibin(ic)) = ielement
                   endif
                   n_occ(ic) = n_occ(ic) +1
                   ! Update total mass of node, and centre of mass
                   m_node(ic) = m_node(ic)+xyzmh(4,ielement)

                   do j=1,3
                      com_node(j,ic) = com_node(j,ic) + xyzmh(j,ielement)*xyzmh(4,ielement)
                   enddo                  

                   partbin(ielement) = ic

                    IF(partbin(ielement)==0) then
                    print*, "Huh? ", ielement, iphase(ielement), ic, binlog, n_occ(ipar), n_occ(ic)
                    endif
                   !	Reset AABB if necessary
                   do j=1,3

                       if(xyzmh(j,ielement)-2.0d0*xyzmh(5,ielement)<bbr_min(j,ic)) then
                           bbr_min(j,ic) = xyzmh(j,ielement)-2.0d0*xyzmh(5,ielement)
                       endif

                       if(xyzmh(j,ielement)+2.0d0*xyzmh(5,ielement)>bbr_max(j,ic)) then
                           bbr_max(j,ic) = xyzmh(j,ielement)+2.0d0*xyzmh(5,ielement)
                       endif

                   enddo
                   !$OMP END CRITICAL 

               endif
              !	End binning for particle ielement

           enddo
           !$OMP END DO
           !$OMP END PARALLEL

                  !	End binning for cell ic
       enddo

    do ic=inode+1,inode+8
        !	iv) Find child nodes that are not leaf nodes (allow them to have children)
        if(n_occ(ic)>nodemax) then
            n_child(ic)=8
        endif

        !	v) Close nodes that satisfy leaf conditions
        if(n_occ(ic)<=nodemax) then

            n_child(ic)=0
            n_leaf = n_leaf +1
            is_leaf(ic)=1
            n_open = n_open-1

        endif

    enddo

    inode = inode+8

!	WRITE(*,*), ipar,n_open

enddo
!	end do WHILE

! Finish centre of mass calculation

do ic=1,n_node
   if(m_node(ic)>0.0) then
      com_node(:,ic) = com_node(:,ic)/m_node(ic)
   endif
enddo


!	Write data to file

print*, 'Tree construction complete: total number of nodes ', inode
print*, 'Total number of leaf nodes: ', n_leaf
WRITE(*, '("Writing tree to files ",A20,"   ",A20)') partbinfile, nodefile
OPEN (2,file=nodefile, form='unformatted')
write(2) sphfile
write(2) inode,nodemax
do ic=1,inode

    WRITE(2) ic, parent(ic),n_child(ic), &
        r_node(1,ic),r_node(2,ic),r_node(3,ic), &
        dr_node(1,ic),dr_node(2,ic),dr_node(3,ic), &
        bbr_min(1,ic), bbr_max(1,ic),bbr_min(2,ic), &
        bbr_max(2,ic),bbr_min(3,ic), bbr_max(3,ic), &
        n_occ(ic), m_node(ic), &
        com_node(1,ic), com_node(2,ic), com_node(3,ic)
enddo
CLOSE(2)

open(2,file=partbinfile,form='unformatted')
write(2) sphfile
do ielement = 1, nelement
    write(2) ielement, iphase(ielement), partbin(ielement)
enddo
close(2)
endif

deallocate(ibin)

return
end subroutine make_octree
