#####################################################
###                                               ###
###         Makefile for sph_tidal_tensor         ###
###  				            	  ###
###                                               ###
###         Duncan H. Forgan (8/8/2014)           ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables:
FC     = gfortran

# For big-endian files use these flags
#FFLAGS = -O3 -frecord-marker=4 -fconvert=swap -fdefault-real-8 -Wunused -fbounds-check
 
# For little-endian files generated use these flags
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8 -Wunused -fbounds-check


# Create object files:
#%.o: %.f
#	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = sphgravmodule.f90 tree_module.f90 sphkernelmodule.f90 main.f90 \
	calc_velocityshear_tensor.f90 jacobi_eigenvalue.f90 deallocate_memory.f90 \
	ktable.f90 make_octree.f90 \
	neighbours_octree.f90 neighbours_brute.f90 \
	rdump_sph_iab.f90 rdump_grav.f90 read_neighbours.f90 write_neighbours.f90


OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: sph_velocityshear_tensor

sph_velocityshear_tensor:  $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)


# Clean statements:
clean: 
	\rm *.o *.mod sph_velocityshear_tensor

# End Makefile
