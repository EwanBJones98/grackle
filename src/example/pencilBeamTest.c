/**************************************************
 * 
 * Testing code to see the optimal array dimensions
 * for iterating over a 3D grid
 * 
 * Based on c_example.c
 * 
 * Author: Ewan Jones
 * Date: August 2023
 * 
***************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <grackle.h>

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16

//? Function to write stats from supplied array to file.
// Currently finds average, max and min.
void write_stats_from_array(double* arr, int arr_len, int* arr_dims, int num_iter, char* functionName)
{
    FILE *fptr;
    int write_header = 0;
    // See if the file already exists, if not then create it and write the table header.
    if ((fptr = fopen("pencilBeamTimings.txt", "r")) == NULL){
        write_header = 1;
    }
    // Open file to write results to.
    if ((fptr = fopen("pencilBeamTimings.txt", "a")) == NULL){
        fprintf(stderr, "Error opening file!\n");
        exit(1);
    }
    // Write header if file was just created.
    if (write_header){
        fprintf(fptr, "%-30s\t%-25s\t%-12s\t%-10s\t%-10s\t%-10s\t%-10s\t%-20s\t%-12s\t%-12s\n",
        "Function", "Number of Iterations", "Array Size", "x-length", "y-length", "z-length",
        "Mean Time", "Standard Deviation", "Minimum Time", "Maximum Time");
    }

    double x, average_x=0, max_x=0, min_x=0;
    // Calculate average, max and min x of array.
    for (int i = 0; i < arr_len; i++){
        // Get x from array.
        x = arr[i];
        // See if x is the max.
        if (max_x == 0 || x > max_x){
            max_x = x;
        }
        // See if x is the min.
        if (min_x == 0 || x < min_x){
            min_x = x;
        }
        // Add x to calculation of average.
        average_x += x;
    }
    // Complete calculation of average x.
    average_x /= arr_len;

    // Calculate sample standard deviation of x.
    double stddev_x = 0;
    for (int i = 0; i < arr_len; i++){
        stddev_x += pow(x - average_x, 2);
    }
    stddev_x = sqrt(stddev_x / (arr_len - 1));

    // Write results to the table.
    fprintf(fptr, "%-30s\t%-25d\t%-12d\t%-10d\t%-10d\t%-10d\t%-10f\t%-20f\t%-14f\t%-14f\n",
            functionName, num_iter, arr_dims[0]*arr_dims[1]*arr_dims[2], arr_dims[0], arr_dims[1],
            arr_dims[2], average_x, stddev_x, min_x, max_x);
}

//? Function to set units of given code_units struct
int setup_units(code_units* grackle_units, double initial_redshift)
{
    // First, set up the units system.
    // These are conversions from code units to cgs.
    grackle_units->comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
    grackle_units->density_units = 1.67e-24;
    grackle_units->length_units = 1.0;
    grackle_units->time_units = 1.0e12;
    grackle_units->a_units = 1.0; // units for the expansion factor
    // Set expansion factor to 1 for non-cosmological simulation.
    grackle_units->a_value = 1. / (1. + initial_redshift) / grackle_units->a_units;
    set_velocity_units(grackle_units);

    return 1;
}

//? Function to allocate memory for chemistry fields
int initialize_chemistry_fields(int total_field_size, chemistry_data* my_grackle_data,
                                grackle_field_data* grackle_fields, code_units* grackle_units)
{
    grackle_fields->density         = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->internal_energy = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->x_velocity      = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->y_velocity      = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->z_velocity      = malloc(total_field_size * sizeof(gr_float));
    // for primordial_chemistry >= 1
    grackle_fields->HI_density      = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->HII_density     = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->HeI_density     = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->HeII_density    = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->HeIII_density   = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->e_density       = malloc(total_field_size * sizeof(gr_float));
    // for primordial_chemistry >= 2
    grackle_fields->HM_density      = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->H2I_density          = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->H2II_density    = malloc(total_field_size * sizeof(gr_float));
    // for primordial_chemistry >= 3
    grackle_fields->DI_density      = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->DII_density     = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->HDI_density     = malloc(total_field_size * sizeof(gr_float));
    // for metal_cooling = 1
    grackle_fields->metal_density   = malloc(total_field_size * sizeof(gr_float));

    // volumetric heating rate (provide in units [erg s^-1 cm^-3])
    grackle_fields->volumetric_heating_rate = malloc(total_field_size * sizeof(gr_float));
    // specific heating rate (provide in units [egs s^-1 g^-1]
    grackle_fields->specific_heating_rate = malloc(total_field_size * sizeof(gr_float));

    // radiative transfer ionization / dissociation rate fields (provide in units [1/s])
    grackle_fields->RT_HI_ionization_rate = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->RT_HeI_ionization_rate = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->RT_HeII_ionization_rate = malloc(total_field_size * sizeof(gr_float));
    grackle_fields->RT_H2_dissociation_rate = malloc(total_field_size * sizeof(gr_float));
    // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
    grackle_fields->RT_heating_rate = malloc(total_field_size * sizeof(gr_float));

    double tiny_number = 1e-20;
    for (int i = 0; i < total_field_size; i++) {
        grackle_fields->density[i] = 1.0;
        grackle_fields->HI_density[i] = my_grackle_data->HydrogenFractionByMass * grackle_fields->density[i];
        grackle_fields->HII_density[i] = tiny_number * grackle_fields->density[i];
        grackle_fields->HM_density[i] = tiny_number * grackle_fields->density[i];
        grackle_fields->HeI_density[i] = (1.0 - my_grackle_data->HydrogenFractionByMass) *
        grackle_fields->density[i];
        grackle_fields->HeII_density[i] = tiny_number * grackle_fields->density[i];
        grackle_fields->HeIII_density[i] = tiny_number * grackle_fields->density[i];
        grackle_fields->H2I_density[i] = tiny_number * grackle_fields->density[i];
        grackle_fields->H2II_density[i] = tiny_number * grackle_fields->density[i];
        grackle_fields->DI_density[i] = 2.0 * 3.4e-5 * grackle_fields->density[i];
        grackle_fields->DII_density[i] = tiny_number * grackle_fields->density[i];
        grackle_fields->HDI_density[i] = tiny_number * grackle_fields->density[i];
        grackle_fields->e_density[i] = tiny_number * grackle_fields->density[i];
        // solar metallicity
        grackle_fields->metal_density[i] = my_grackle_data->SolarMetalFractionByMass *
        grackle_fields->density[i];

        grackle_fields->x_velocity[i] = 0.0;
        grackle_fields->y_velocity[i] = 0.0;
        grackle_fields->z_velocity[i] = 0.0;

        // initilize internal energy (here 1000 K for no reason)
        grackle_fields->internal_energy[i] = 1000. / get_temperature_units(grackle_units);

        grackle_fields->volumetric_heating_rate[i] = 0.0;
        grackle_fields->specific_heating_rate[i] = 0.0;

        grackle_fields->RT_HI_ionization_rate[i] = 0.0;
        grackle_fields->RT_HeI_ionization_rate[i] = 0.0;
        grackle_fields->RT_HeII_ionization_rate[i] = 0.0;
        grackle_fields->RT_H2_dissociation_rate[i] = 0.0;
        grackle_fields->RT_heating_rate[i] = 0.0;
    }

    return 1;
}


int main(int argc, int *argv[])
{   
    //? Read in command-line arguments
    int gridDimX, gridDimY, gridDimZ;
    if (argc==4){
        gridDimX = (int) strtol(argv[1], NULL, 10);
        gridDimY = (int) strtol(argv[2], NULL, 10);
        gridDimZ = (int) strtol(argv[3], NULL, 10);
    } else {
        fprintf(stderr, "Three arguments must be supplied: gridDimX, gridDimY, gridDimZ\nPlease try again.\n");
        return 0;
    }

    //!-------------------------------------------------------------------
    //! Initialising the chemistry structures and parameters.
    //!-------------------------------------------------------------------

    //? Set up unit system
    code_units my_units;
    double initial_redshift = 0.;
    if (!setup_units(&my_units, initial_redshift)){
        fprintf(stderr, "setup_units failed!");
        return 0;
    }

    //? Set up chemistry parameters
    // Create struct to store them and set them to their default values.
    chemistry_data *my_grackle_data;
    my_grackle_data = malloc(sizeof(chemistry_data));
    // Set them to default
    if (set_default_chemistry_parameters(my_grackle_data) == 0) {
        fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
        return 0;
    }
    // Change parameters as required.
    my_grackle_data->use_grackle = 1;            // chemistry on
    my_grackle_data->with_radiative_cooling = 1; // cooling on
    my_grackle_data->primordial_chemistry = 3;   // molecular network with H, He, D
    my_grackle_data->dust_chemistry = 1;         // dust processes
    my_grackle_data->metal_cooling = 1;          // metal cooling on
    my_grackle_data->UVbackground = 1;           // UV background on
    my_grackle_data->grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5"; // data file

    //? Finally, initialize the chemistry object.
    // Create struct to store chemical rates
    chemistry_data_storage my_rates;    
    if (local_initialize_chemistry_data(my_grackle_data, &my_rates, &my_units) == 0) {
        fprintf(stderr, "Error in initialize_chemistry_data.\n");
        return 0;
    }

    //? Create struct to store the grackle field data
    grackle_field_data my_fields;

    //? Set up particle grid
    // Scalars
    my_fields.grid_rank = 3; // Number of grid dimensions
    my_fields.grid_dx = 0.0; // This is only ised for the H2 self-shielding approximation
    // Arrays: first assign memory
    my_fields.grid_dimension = malloc(my_fields.grid_rank * sizeof(int)); // Number of particles in each dimension
    my_fields.grid_start = malloc(my_fields.grid_rank * sizeof(int)); // Starting index in each dimension
    my_fields.grid_end = malloc(my_fields.grid_rank * sizeof(int)); // Ending index in each dimension    
    // Now assign values
    int field_size[3] = {gridDimX,gridDimY,gridDimZ}; // Number of particles in each dimension
    int total_field_size = 1; // Total number of particles to be updated in loop below
    for (int i=0; i < my_fields.grid_rank; i++){
        my_fields.grid_dimension[i] = field_size[i];
        my_fields.grid_start[i] = 0;
        my_fields.grid_end[i] = field_size[i] - 1;
        total_field_size *= field_size[i];
    } 

    //? Allocate memory for chemistry fields and set them to initial values
    if (!(initialize_chemistry_fields(total_field_size, my_grackle_data, &my_fields, &my_units))){
        fprintf(stderr, "Error in initialize_chemistry_fields.\n");
        return 0;
    }


    //!-------------------------------------------------------------------
    //! Calling the chemistry solving routines.
    //! Select which ones you want to run with the flags below.
    //!-------------------------------------------------------------------

    //? Flags to select which computations you want to run.
    int run_calc_cooling_time, run_calc_temperature,
        run_calc_pressure, run_calc_gamma, run_calc_dust_temperature;
    run_calc_cooling_time       = 1;
    run_calc_temperature        = 0;
    run_calc_pressure           = 0;
    run_calc_gamma              = 0;
    run_calc_dust_temperature   = 0;

    //? Set up benchmarking settings.
    int number_of_iterations = 100;
    double elapsed_times[number_of_iterations], average_elapsed_time;
    clock_t start_time, end_time;

    //? Run chemistry solver, this has to be done.
    // Set some timestep to evolve over.
    double dt = 3.15e7 * 1e6 / my_units.time_units;
    // Run calculation.
    if (!local_solve_chemistry(my_grackle_data, &my_rates, &my_units, &my_fields, dt)){
        fprintf(stderr, "Error in solve_chemistry.\n");
        return 0;
    }

    //? Run cooling time calculation.
    if (run_calc_cooling_time){
        // Allocate array to store results in.
        gr_float *cooling_time;
        cooling_time = malloc(total_field_size * sizeof(gr_float));
        // Run and time calculation.
        for (int i=0; i < number_of_iterations; i++){
            start_time = clock();
            if (!local_calculate_cooling_time(my_grackle_data, &my_rates, &my_units, &my_fields,
                                        cooling_time)){
            fprintf(stderr, "Error in calculate_cooling_time.\n");
            return 0;
            }
            end_time = clock();
            elapsed_times[i] = ((double) end_time - start_time) / CLOCKS_PER_SEC;
        }
        // Write stats from this test into file
        write_stats_from_array(elapsed_times, number_of_iterations, field_size, number_of_iterations, "local_calculate_cooling_time");
    }
    
    //? Run temperature calculation.
    if (run_calc_temperature){
        // Allocate array to store results in.
        gr_float *temperature;
        temperature = malloc(total_field_size * sizeof(gr_float));
        // Run calculation.
        if (local_calculate_temperature(my_grackle_data, &my_rates, &my_units, &my_fields,
                                    temperature) == 0) {
            fprintf(stderr, "Error in calculate_temperature.\n");
            return 0;
        }
    }
    
    //? Run pressure calculation.
    if (run_calc_pressure){
        // Create array to store results in.
        gr_float *pressure;
        pressure = malloc(total_field_size * sizeof(gr_float));
        // Run calculation.
        if (!local_calculate_pressure(my_grackle_data, &my_rates, &my_units, &my_fields,
                                pressure)) {
            fprintf(stderr, "Error in calculate_pressure.\n");
            return 0;
        }
    }
    
    //? Run gamma calculation.
    if (run_calc_gamma){
        // Create array to store results in.
        gr_float *gamma;
        gamma = malloc(total_field_size * sizeof(gr_float));
        // Run calculation.
        if (!local_calculate_gamma(my_grackle_data, &my_rates, &my_units, &my_fields,
                            gamma)) {
            fprintf(stderr, "Error in calculate_gamma.\n");
            return 0;
        }
    }

    //? Run dust temperature calculation.
    if (run_calc_dust_temperature){
        // Create array to store results in.
        gr_float *dust_temperature;
        dust_temperature = malloc(total_field_size * sizeof(gr_float));
        // Run calculation.
        if (!local_calculate_dust_temperature(my_grackle_data, &my_rates, &my_units, &my_fields,
                            dust_temperature)) {
            fprintf(stderr, "Error in calculate_dust_temperature.\n");
            return 0;
        }
    }

    return 1;
}