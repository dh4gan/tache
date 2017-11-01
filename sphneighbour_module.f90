module sphneighbourdata

  !-----------------------------------------------------------------------
  ! Data module for saving data on SPH neighbours
  ! DHF 01/09/2009
  !-----------------------------------------------------------------------
  
  implicit none
  save
 
  !-------------------Single values---------------------------------------
  
  ! Integers
  
  integer ::  inode,n_node
  real :: hmin,xmax,ymax,zmax
  real:: meanneigh, sdneigh, neighcrit
  integer, parameter :: nodemax = 10
  integer, parameter :: neighmax = 150
  real,parameter:: tolerance = 5.0 ! Search radius (in smoothing lengths) for neighbours using octree (minimum value is 2.0)
  
  
  !----------------------Arrays---------------------------------------
  
  integer, allocatable, dimension(:) :: member
  integer, allocatable, dimension(:) :: ray,n_occ,parent,n_child,nneigh,partbin
  integer, allocatable, dimension(:,:) :: child, occ, neighb		
  real,allocatable, dimension(:) :: b,t_min, m_node
  real, allocatable, dimension(:,:) :: r_node, dr_node,com_node,bbr_min, bbr_max
  
 ! Data for regular grids only (occ, n_occ, xmax,ymax,zmax shared with octree)
  
  integer :: ncells, ngridx,ngridy,ngridz, cellmax, nelementiclelist, ncellrange
  integer, allocatable,dimension(:) :: particlelist, cellist,isortcellpart

  real :: dgrid
  

contains

  subroutine get_cellID(cellID, icell,jcell,kcell,nx,ny,nz)
    ! Simple routine to get a unique ID for every cell in a 3D Cartesian grid
    
    integer :: cellID, icell, jcell,kcell
    integer :: nx,ny,nz
    
    if(icell >nx .or. jcell >ny .or. kcell>nz) then
       print*, 'Error in get_cellID! cell index exceeds grid limits'
       print*, 'Cell: ',icell,jcell,kcell
       print*, 'Grid Limits: ', nx,ny,nz
    endif
    
    if(nz>1) then
       cellID = icell + nx*(jcell-1) + (kcell-1)*nx*ny
    else
       cellID = icell + (jcell-1)*nx
    endif
    
    return
  end subroutine get_cellID
  
end module sphneighbourdata
