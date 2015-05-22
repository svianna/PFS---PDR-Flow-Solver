#stl reader

FC=gfortran

#DEBUG
#FFLAGS=-c -fbounds-check
#RUN
FFLAGS=-c -O5
FDEP=common_block.o

FDEP=common_block.o

FSOURCES=main_EULER.o common_block.o read_data.o inguess.o read_porosity.o check_grid.o check_porosity.o set_timestep.o crude_guess.o set_others.o inout_bc.o flux_ma.o wall_bc.o flux_energy.o flux_mom.o sum_fluxes.o smooth.o check_conv.o output.o read_parameters.o

all: MAIN_EULER

MAIN_EULER: $(FSOURCES)
	$(FC) $(FSOURCES) -o MAIN_EULER 

main_EULER.o: main_EULER.f90 $(FDEP)
	$(FC) $(FFLAGS) main_EULER.f90

common_block.o: common_block.f90 $(FDEP)
	$(FC) $(FFLAGS) common_block.f90

read_data.o: read_data.f90
	$(FC) $(FFLAGS) read_data.f90

inguess.o: inguess.f90
	$(FC) $(FFLAGS) inguess.f90

check_porosity.o: check_porosity.f90
	$(FC) $(FFLAGS) check_porosity.f90

read_porosity.o: read_porosity.f90
	$(FC) $(FFLAGS) read_porosity.f90

check_grid.o: check_grid.f90
	$(FC) $(FFLAGS) check_grid.f90

crude_guess.o: crude_guess.f90
	$(FC) $(FFLAGS) crude_guess.f90

set_timestep.o: set_timestep.f90
	$(FC) $(FFLAGS) set_timestep.f90

set_others.o: set_others.f90
	$(FC) $(FFLAGS) set_others.f90

inout_bc.o: inout_bc.f90
	$(FC) $(FFLAGS) inout_bc.f90

flux_ma.o: flux_ma.f90
	$(FC) $(FFLAGS) flux_ma.f90

wall_bc.o: wall_bc.f90
	$(FC) $(FFLAGS) wall_bc.f90

flux_energy.o: flux_energy.f90
	$(FC) $(FFLAGS) flux_energy.f90

flux_mom.o: flux_mom.f90
	$(FC) $(FFLAGS) flux_mom.f90

sum_fluxes.o: sum_fluxes.f90
	$(FC) $(FFLAGS) sum_fluxes.f90

smooth.o: smooth.f90
	$(FC) $(FFLAGS) smooth.f90

check_conv.o: check_conv.f90
	$(FC) $(FFLAGS) check_conv.f90

output.o: output.f90
	$(FC) $(FFLAGS) output.f90

read_parameters.o: read_parameters.f90
	$(FC) $(FFLAGS) read_parameters.f90

clean:
	rm -rf *.o
