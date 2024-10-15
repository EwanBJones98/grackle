/***********************************************************************
/
/ Example executable using libgrackle
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <sys/time.h>

#include <grackle.h>

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16

void initialise_structs(double initial_redshift,
                        double field_size,
                        chemistry_data *my_chemistry_data,
                        chemistry_data_storage *my_rates, 
                        grackle_field_data *my_fields,
                        code_units *my_units)
{

  double tiny_number = 1e-20;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  my_units->comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units->density_units = 1.67e-24;
  my_units->length_units = 1.0;
  my_units->time_units = 1.0e12;
  my_units->a_units = 1.0; // units for the expansion factor
  // Set expansion factor to 1 for non-cosmological simulation.
  my_units->a_value = 1. / (1. + initial_redshift) / my_units->a_units;
  set_velocity_units(my_units);

  // Second, create initialise chemistry to default value.
  my_chemistry_data = malloc(sizeof(my_chemistry_data));
  set_default_chemistry_parameters(my_chemistry_data);
  my_chemistry_data->use_grackle = 1;            // chemistry on
  my_chemistry_data->with_radiative_cooling = 1; // cooling on
  my_chemistry_data->primordial_chemistry = 3;   // molecular network with H, He, D
  my_chemistry_data->dust_chemistry = 1;         // dust processes
  my_chemistry_data->metal_cooling = 1;          // metal cooling on
  my_chemistry_data->UVbackground = 1;           // UV background on
  my_chemistry_data->grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5"; // data file

  // Finally, initialize the rates.
  if (local_initialize_chemistry_data(my_chemistry_data, my_rates, my_units) == 0)
    fprintf(stderr, "Error in local_initialize_chemistry!\n");

  // Initialise struct for storing grackle field data
  if (gr_initialize_field_data(my_fields) == 0)
    fprintf(stderr, "Error in gr_initialize_field_data!\n");

  my_fields->density         = malloc(field_size * sizeof(gr_float));
  my_fields->internal_energy = malloc(field_size * sizeof(gr_float));
  my_fields->x_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields->y_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields->z_velocity      = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 1
  my_fields->HI_density      = malloc(field_size * sizeof(gr_float));
  my_fields->HII_density     = malloc(field_size * sizeof(gr_float));
  my_fields->HeI_density     = malloc(field_size * sizeof(gr_float));
  my_fields->HeII_density    = malloc(field_size * sizeof(gr_float));
  my_fields->HeIII_density   = malloc(field_size * sizeof(gr_float));
  my_fields->e_density       = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 2
  my_fields->HM_density      = malloc(field_size * sizeof(gr_float));
  my_fields->H2I_density     = malloc(field_size * sizeof(gr_float));
  my_fields->H2II_density    = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 3
  my_fields->DI_density      = malloc(field_size * sizeof(gr_float));
  my_fields->DII_density     = malloc(field_size * sizeof(gr_float));
  my_fields->HDI_density     = malloc(field_size * sizeof(gr_float));
  // for metal_cooling = 1
  my_fields->metal_density   = malloc(field_size * sizeof(gr_float));

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields->volumetric_heating_rate = malloc(field_size * sizeof(gr_float));
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields->specific_heating_rate = malloc(field_size * sizeof(gr_float));

  // radiative transfer ionization / dissociation rate fields (provide in units [1/s])
  my_fields->RT_HI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_HeI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_HeII_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_H2_dissociation_rate = malloc(field_size * sizeof(gr_float));
  // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
  my_fields->RT_heating_rate = malloc(field_size * sizeof(gr_float));

  // set temperature units
  double temperature_units = get_temperature_units(my_units);

  for (int i=0; i<field_size; i++)
  {
    my_fields->density[i] = 1.0;
    my_fields->HI_density[i] = my_chemistry_data->HydrogenFractionByMass * my_fields->density[i];
    my_fields->HII_density[i] = tiny_number * my_fields->density[i];
    my_fields->HM_density[i] = tiny_number * my_fields->density[i];
    my_fields->HeI_density[i] = (1.0 - my_chemistry_data->HydrogenFractionByMass) *
                                  my_fields->density[i];
    my_fields->HeII_density[i] = tiny_number * my_fields->density[i];
    my_fields->HeIII_density[i] = tiny_number * my_fields->density[i];
    my_fields->H2I_density[i] = tiny_number * my_fields->density[i];
    my_fields->H2II_density[i] = tiny_number * my_fields->density[i];
    my_fields->DI_density[i] = 2.0 * 3.4e-5 * my_fields->density[i];
    my_fields->DII_density[i] = tiny_number * my_fields->density[i];
    my_fields->HDI_density[i] = tiny_number * my_fields->density[i];
    my_fields->e_density[i] = tiny_number * my_fields->density[i];
    // solar metallicity
    my_fields->metal_density[i] = my_chemistry_data->SolarMetalFractionByMass *
                                  my_fields->density[i];

    my_fields->x_velocity[i] = 0.0;
    my_fields->y_velocity[i] = 0.0;
    my_fields->z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
    my_fields->internal_energy[i] = 1000. / temperature_units;

    my_fields->volumetric_heating_rate[i] = 0.0;
    my_fields->specific_heating_rate[i] = 0.0;

    my_fields->RT_HI_ionization_rate[i] = 0.0;
    my_fields->RT_HeI_ionization_rate[i] = 0.0;
    my_fields->RT_HeII_ionization_rate[i] = 0.0;
    my_fields->RT_H2_dissociation_rate[i] = 0.0;
    my_fields->RT_heating_rate[i] = 0.0;
  }
}

double run_test(int grid_rank, int *grid_dimensions)
{
  // Initialise timing struct
  struct timeval time_start, time_end;

  // Initialise grackle structs
  chemistry_data *my_data;
  chemistry_data_storage my_rates;
  grackle_field_data my_fields;
  code_units my_units;

  // Set initial redshift (for internal units).
  double initial_redshift = 0.;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  my_fields.grid_rank = grid_rank;
  my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  my_fields.grid_dimension = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_start = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_end = malloc(my_fields.grid_rank * sizeof(int));

  int field_size = 1;
  for (int i=0; i<my_fields.grid_rank; i++)
    field_size *= grid_dimensions[i];

  for (int i = 0; i < my_fields.grid_rank; i++)
  {
    my_fields.grid_dimension[i] = grid_dimensions[i]; // the active dimension not including ghost zones.
    my_fields.grid_start[i] = 0;
    my_fields.grid_end[i] = grid_dimensions[i] - 1;
  }

  initialise_structs(initial_redshift, field_size, my_data, &my_rates, &my_fields, &my_units);

  /*********************************************************************
  / Calling the chemistry solver
  *********************************************************************/

  // Evolving the chemistry.
  // some timestep
  double dt = 3.15e7 * 1e6 / my_units.time_units;

  gettimeofday(&time_start, NULL);
  if (local_solve_chemistry(my_data, &my_rates, &my_units, &my_fields, dt) == 0)
  {
    fprintf(stderr, "Error in solve_chemistry.\n");
  }
  gettimeofday(&time_end, NULL);

  double elapsed_time = (time_end.tv_sec - time_end.tv_sec) + (time_end.tv_usec - time_end.tv_usec)*1e-6;

  return elapsed_time;
}

int main(int argc, char *argv[])
{
  int log2_field_size = 19;
  int num_timing_iter = 1000;

  int field_size = pow(2, log2_field_size);

  FILE *fptr;
  fptr = fopen("test_results.txt", "w");

  if (fptr == NULL)
  {
      printf("Error opening file!\n");
      exit(0);
  }
  fprintf(fptr, "#i size, j size, k size, number of iterations, mean time (s)\n");

  int i_size, j_size, k_size;
  int log2_i_size, log2_j_size, log2_k_size;

  double mean_test_runtime;
  for (log2_i_size = 0; log2_i_size <= log2_field_size; log2_i_size)
  {
    log2_j_size = log2_field_size - log2_i_size;
    log2_k_size = 1;

    i_size = (int) pow(2, log2_i_size);
    j_size = (int) pow(2, log2_j_size);
    k_size = (int) pow(2, log2_k_size);

    int grid_dimensions[] = {i_size, j_size, k_size};

    mean_test_runtime = 0.;
    for (int iter=0; iter < num_timing_iter; iter++)
      mean_test_runtime += run_test(3, grid_dimensions) / num_timing_iter;

    fprintf(fptr, "%d, %d, %d, %d, %f\n", i_size, j_size, k_size, num_timing_iter, mean_test_runtime);
  }

  fclose(fptr);


  return 0;
}
