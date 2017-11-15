#####################################################
###                                               ###
###         	Makefile for TACHE	          ###
###  				            	  ###
###                                               ###
###         Duncan H. Forgan (1/11/2017)          ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables:
FC     = gfortran

# For big-endian files use these flags
#FFLAGS = -O3 -frecord-marker=4 -fconvert=swap -fdefault-real-8 -Wunused -fbounds-check
 
# For little-endian files use these flags
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8 -Wunused -fbounds-check


# Create object files:
#%.o: %.f
#	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

#
# Sources for the main tache program
#
SRCTACHE = sphdata_module.f90 sphkernel_module.f90 \
	sphneighbour_module.f90 tachedata_module.f90  main.f90 \
	calc_eigenvalues.f90 calc_grav_brute.f90 calc_grav.f90 \
	calc_grav_from_pot.f90 calc_grav_tree.f90 \
	calc_tidal_tensor.f90 calc_velocityshear_tensor.f90 \
	checksetup.f90 classify_by_eigenvalue.f90 \
	compute_tensor.f90 deallocate_memory.f90 \
	element_percent_complete.f90 \
	get_SPH_neighbours.f90 initial.f90 jacobi_eigenvalue.f90 \
	ktable.f90 neighbours_octree.f90 make_grid_sphdata.f90 \
	make_octree.f90 mark_children.f90 neighbours_brute.f90 \
	neighbours_grid.f90 particle_forces.f90 read_dump.f90 \
	rdump_grav.f90 	rdump_sphNG_wkmr.f90 \
	rdump_sphNG_iab.f90  recalc_density.f90 \
	read_neighbours.f90 sort2.f90 wdump_grav.f90 \
	wdump_sph.f90 write_eigendata.f90 \
	write_splitdump.f90 write_neighbours.f90

OBJTACHE    = $(SRCTACHE:.f90=.o)

#
# Sources for the spiralfind program
#

SRCSPIRAL = spiraldata_module.f90  \
	 spiralfind.f90  \
	classify_by_eigenvalue.f90  \
	element_percent_complete.f90 \
	extract_spiral_elements.f90 \
	read_eigenvalue_file.f90 \
	sort2.f90 

OBJSPIRAL = $(SRCSPIRAL:.f90=.o)


# Create executable files:
build: tache

tache:  $(OBJTACHE)
	$(FC) $(FFLAGS) -o $@ $(OBJTACHE)

spiralfind: $(OBJSPIRAL)
	    $(FC) $(FFLAGS) -o $@ $(OBJSPIRAL)

# Clean statements:
clean: 
	\rm *.o *.mod tache spiralfind

# End Makefile
