#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "phys_constants.h"

// Function prototype for calc_tdust_1d
int calc_tdust_1d(double *tdust, double *tgas, double *nh, double *gasgr,
                    double *isrf, int *itmask, double trad, int j, int k,
                    chemistry_data_storage *my_rates,
                    grackle_field_data *my_fields);

/*****************************************
*_________ FUNCTION CALC_TDUST_3D_________
* 
*  AUTHORS:
*   Original Fortran function written by
*   Britton Smith (July 2011).
*
*   Ported to C by Ewan Jones (Jan 2022).
*
*   PURPOSE:
!    Calculate dust heat balance to get the dust temperature
*   
******************************************/
int calc_tdust_3d(gr_float *gas_temperature, gr_float *dust_temperature, double *temperature_units,
                    double *co_length_units, double *co_density_units, grackle_field_data *my_fields,
                    chemistry_data *my_chemistry, chemistry_data_storage *my_rates, 
                    code_units *my_units)
{
    /*============== INIT ==============*/
    int in = my_fields->grid_dimension;
    int itmask[in], indixe[in];
    double myisrf[in], nh[in], tgas[in], tdust[in], logtem[in], t1[in], t2[in], tdef[in], gasgr[in];
    /*==================================*/

    //* Set log values of start and dend of lookup tables
    double logtem0 = log(my_chemistry->TemperatureStart);
    double logtem9 = log(my_chemistry->TemperatureEnd);
    double dlogtem = ( log(my_chemistry->TemperatureEnd) - log(my_chemistry->TemperatureStart) )
                        / ( (float) my_chemistry->NumberOfTemperatureBins );

    //* Set units
    double tbase1   = my_units->time_units;
        // co_length_units is [x]*a = [x]*[a]*a'
    double xbase1   = *co_length_units / (my_units->a_value * my_units->a_units);
        // co_density units is [dens]/a^3 = [dens]/([a]*a')^3
    double dbase1   = *co_density_units * pow(my_units->a_value * my_units->a_units, 3);
    double coolunit = ( pow(my_units->a_units,5) * pow(xbase1,2) * pow(mh,2) ) /
                        ( pow(tbase1,3) * dbase1 );
    double zr       = 1.0 / (my_units->a_value*my_units->a_units) - 1.0;

    //* Set CMB temperature
    double trad = 2.73 * (1.0 + zr);

    //* Loop over zones, and do enture i-column in one go
    //! Not sure if I need the +1's here as it may be from fortran indexing
    int dk = my_fields->grid_end+2 - my_fields->grid_start+2 + 1;
    int dj = my_fields->grid_end+1 - my_fields->grid_start+1 + 1;

    // Parallelise the k and j loops with OpenMP.
    // Flat j and k loops for better parallelism.
    #ifdef _OPENMP
    //! This opemMP code is commented out in the fortran code, ask Britton
    //! and if its needed copy from the fortran version (not the below, its
    //! just a placeholder!)
    //#pragma omp parallel for default(None) private()
    #endif
        for (int t = 0; t < dk*dj; t++) {
            //! Again not sure if you need the +1 here, don't think you do 
            //! as its used directly as an index
            int k = t/dj + my_fields->grid_start+2; // + 1;
            int j = t%dj + my_fields->grid_start+1; // + 1;
            
            int is = my_fields->grid_start;
            int ie = my_fields->grid_end;
            for (int i = is; i <= ie; i++) {
                //Set itmask to all true
                itmask[i] = TRUE;

                // Compute interstellar radiation field
                if (my_chemistry->use_isrf_field > 0) {
                    myisrf[i] = my_fields->isrf_habing[i,j,k];
                } else {
                    myisrf[i] = my_chemistry->interstellar_radiation_field;
                }

                // Compute hydrogen number density
                nh[i] = my_fields->HI_density[i,j,k] + my_fields->HII_density[i,j,k];
                if (my_chemistry->primordial_chemistry > 1) {
                    nh[i] = nh[i] + my_fields->H2I_density[i,j,k]
                            + my_fields->H2II_density[i,j,k];
                }

                // We have not converted to proper, so use comoving density units
                nh[i] = nh[i] * *co_density_units / mh;

                // Compute log temperature and truncate if above/below table max/min
                tgas[i]   = gas_temperature[i,j,k];
                logtem[i] = log(tgas[i]);
                logtem[i] = max(logtem[i], logtem0);
                logtem[i] = min(logtem[i], logtem9);

                // Compute index into the table and precompute parts of linear interp
                //! Need to check the zero indexing on this too
                indixe[i] = min(my_chemistry->NumberOfTemperatureBins - 1,
                                max(1, (int) ((logtem[i] - logtem0)/dlogtem) + 1));
                t1[i] = logtem0 + (indixe[i] - 1)*dlogtem;
                t2[i] = logtem0 + (indixe[i])*dlogtem;
                tdef[i] = (logtem[i] - t1[i]) / (t2[i] - t1[i]);

                // Lookup calues and do a linear temperature in log(T)
                // Convert back to cgs
                gasgr[i] = my_rates->gas_grain[indixe[i]] + tdef[i]*(
                            my_rates->gas_grain[indixe[i]+1] -
                            my_rates->gas_grain[indixe[i]]);
                gasgr[i] = gasgr[i] * my_chemistry->local_dust_to_gas_ratio * coolunit / mh;
            }

            //* Compute dust temperature in a slice
            calc_tdust_1d(tdust, tgas, nh, gasgr, myisrf, itmask, trad, j, k,
                            my_rates, my_fields);

            //* Copy slice values back to grid
            for (int i = is; i <= ie; i++) {
                dust_temperature[i,j,k] = tdust[i];
            }
        }

    // End of function
    return SUCCESS;
}