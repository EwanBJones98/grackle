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

/*****************************************
*_________ FUNCTION CALC_KAPPA_GR_________
* 
*  AUTHORS:
*   Original Fortran function written by
*   Britton Smith (Sep 2011).
*
*   Ported to C by Ewan Jones (Jan 2022).
*
*   PURPOSE:
!    Calculate grain plank mean opacity
!
!  INPUTS:
!     in       - i dimension of 3D fields
!
!     tdust    - dust temperature
!
!     is,ie    - start and end indices of active region (zero based)
!
!     itmask   - iteration mask
!
!     t_subl   - grain sublimation temperature
!
!  OUTPUTS:
!     kgr      - opacities
*   
******************************************/
void calc_kappa_gr(double tdust[], double kgr[], int itmask[], int in, int is,
                    int ie, double t_subl)
{
    /*============== INIT ==============*/
    //* Defining constants
    double kgr1 = 4.0e-4, kgr200 = 16.0;
    /*==================================*/

    for (int i = is; i <= ie; i++) {
        if (itmask[i]) {
        // Temperature dependence from Dopcke et al. (2011).
        // Normalised to Omukai (2000).
            if (tdust[i] < 200.0) {
                kgr[i] = kgr1 * pow(tdust[i], 2);
            } else if (tdust[i] < t_subl) {
                kgr[i] = kgr200;
            } else {
                kgr[i] = max(tiny, kgr200 * pow(tdust[i]/1.5e3, -12));
            }
        }
    }
}

/*****************************************
*_________ FUNCTION CALC_GR_BALANCE_________
* 
*  AUTHORS:
*   Original Fortran function written by
*   Britton Smith (Sep 2019).
*
*   Ported to C by Ewan Jones (Jan 2022).
*
*   PURPOSE:
!    Calculate grain heating/cooling balance
!
!  INPUTS:
!     in       - i dimension of 3D fields
!
!     tdust    - dust temperature
!     tgas     - gas temperature
!     kgr      - grain opacity
!     trad4    - CMB temperature to 4th power
!     gasgr    - gas/grain heat transfer rate
!     gamma_isrf - heating from interstellar radiation field
!     nh       - hydrogen number density
!
!     is,ie    - start and end indices of active region (zero based)
!
!     itmask   - iteration mask
!
!
!  OUTPUTS:
!     sol      - heating/cooling balance (heating - cooling)
*   
******************************************/
void calc_gr_balance (double tdust[], double tgas[], double kgr[],
                    double trad4, double _gamma_isrf[], double nh[],
                    int itmask[], double sol[], int in, int is, int ie,
                    chemistry_data_storage *my_rates)
{
    /*============== INIT ==============*/
    //* Defining constants
    double radf = 4.0 * sigma_sb;
    /*==================================*/

    for (int i = is; i <= ie; i++) {
        if (itmask[i]) {
            sol[i] = _gamma_isrf[i] + radf * kgr[i] *
                        (trad4 - pow(tdust[i],4)) +
                        (my_rates->gas_grain[i] * nh[i] *
                        (tgas[i] - tdust[i]));
        }
    }
}

/*****************************************
*_________ FUNCTION CALC_TDUST_1D_________
* 
*  AUTHORS:
*   Original Fortran function written by
*   Britton Smith (Feb 2011).
*
*   Ported to C by Ewan Jones (Jan 2022).
*
*   PURPOSE:
*    Calculate dust temperature.
*
*   INPUTS:
*     tdust    - dust temperature
*     tgas     - gas temperature

*     nh       - H number density

*     isrf     - interstellar radiation field in Habing units
*     trad     - CMB temperature
*
*     j,k      - indices of 1D slice
*
*     itmask   - iteration mask
*   
******************************************/
int calc_tdust_1d(double *tdust, double *tgas, double *nh, double isrf,
                    int *itmask, double trad, int j, int k,
                    chemistry_data_storage *my_rates,
                    grackle_field_data *my_fields)
{
    /*============== INIT ==============*/
    //* Defining constants
    double t_subl = 1.5e3; // Grain sublimation temperature
    double radf = 4.0 * sigma_sb, kgr1 = 4.0e-4, tol = 1.0e-5, bi_tol = 1.0e-3,
            minpert = 1.0e-10;
    int itmax = 50, bi_itmax = 30;
    int in = my_fields->grid_dimension, is = my_fields->grid_dimension+1,
            ie = my_fields->grid_dimension+2;
    double _gamma_isrf[in]; // Internal gamma_isrf array

    //* Local variables
    int i, iter, c_done, c_total, nm_done;
    double pert_i = 1.0e-3;
    double trad4 = pow( max(1.0, trad), 4 );
    // Set total cells for calculation
    c_done = 0, nm_done = 0;
    c_total = ie - is + 1;

    //* Local slice variables
    double kgr[in], kgrplus[in], sol[in], solplus[in], slope[in], tdplus[in],
            tdustnow[in], tdustold[in], pert[in], bi_t_mid[in], bi_t_high[in];
    int nm_itmask[in], bi_itmask[in];
    /*==================================*/

    //* Set local iteration mask and initial guess
    for (int i = is; i <= ie; i++) {
        if (itmask[i]) {
            _gamma_isrf[i] = isrf * my_rates->gamma_isrf;
        }
    }

    for (int i = is; i <= ie; i++) {
        nm_itmask[i] = itmask[i];
        bi_itmask[i] = itmask[i];
        if (nm_itmask[1]) {
            if (trad >= tgas[i]) {
                tdustnow[i] = trad;
                nm_itmask[i] = FALSE;
                bi_itmask[i] = FALSE;
                c_done ++;
                nm_done ++;
            } else if (tgas[i] > t_subl) {
            //  Use bisection if T_gas > grain sublimation temperature.
            nm_itmask[i] = FALSE;
            nm_done ++;
            } else {
                tdustnow[i] = max(trad, pow(_gamma_isrf[i] / radf / kgr1, 0.17));
                pert[i] = pert_i;
            }
        } else {
            c_done ++;
            nm_done ++;
        }
    }

    //* Iterate to convergence with Newton's method
    for (int iter = 0; iter <= itmax; iter++) {
    
        // Loop over slice
        for (int i = is; i <= ie; i++) {
            if (nm_itmask[i]) {
                tdplus[i] = max(1.0e-3, (1.0 + pert[i])*tdustnow[i]);
            }
        }

        // Calculate grain opacities
        calc_kappa_gr(tdustnow, kgr, nm_itmask, in, is, ie, t_subl);
        calc_kappa_gr(tdplus, kgrplus, nm_itmask, in, is, ie, t_subl);

        // Calculate heating/cooling balance
        calc_gr_balance(tdustnow, tgas, kgr, trad4, _gamma_isrf, nh,
                        nm_itmask, sol, in, is, ie, my_rates);
        calc_gr_balance(tdplus, tgas, kgrplus, trad4, _gamma_isrf, nh,
                        nm_itmask, solplus, in, is, ie, my_rates);

        for (int i = is; i <= ie; i++) {
            if (nm_itmask[i]) {
            // Use Newton's method to solve for Tdust
                slope[i] = (solplus[i] - sol[i]) /
                            (pert[i] * tdustnow[i]);
                
                tdustold[i] = tdustnow[i];
                tdustnow[i] = tdustnow[i] - (sol[i] / slope[i]);

                pert[i] = max(min(pert[i], 0.5 * abs(tdustnow[i]
                                - tdustold[i]) / tdustnow[i]),
                                minpert);

                // If negative solution calculated, give up and wait for bisection step.
                if (tdustnow[i] < trad) {
                    nm_itmask[i] = FALSE;
                    nm_done ++;
                } else if (abs(sol[i] / solplus[i]) < tol) {
                // Check for convergence of solution.
                    nm_itmask[i] = FALSE;
                    c_done ++;
                    bi_itmask[i] = FALSE;
                    nm_done ++;
                }
            }
        } // End of loop over slice

        // Check for all cells converged
        if (c_done >= c_total) goto label_converged;

        // Check for all cells done with Newon method.
        // This includes attempts where a negative solution was found.
        if (nm_done >= c_total) goto label_nmdone;
    } // End of loop for Newton's method

    label_nmdone:
    
    // If iteration count exceeded, try once more with bisection
    if (c_done < c_total) {
        for (int i = is; i <= ie; i++){
            if (bi_itmask[i]) {
                tdustnow[i] = trad;
                bi_t_high[i] = tgas[i];
            }
        }

        for (int iter = 1; iter <= bi_itmax; iter++) {
            for (int i = is; i = ie; i++) {
                if (bi_itmask[i]) {
                    bi_t_mid[i] = 0.5 * (tdustnow[i] + bi_t_high[i]);
                    if (iter == 1) {
                        bi_t_mid[i] = min(bi_t_mid[i], t_subl);
                    }
                }
            }

            calc_kappa_gr(bi_t_mid, kgr, bi_itmask, in, is, ie, t_subl);
            
            calc_gr_balance(bi_t_mid, tgas, kgr, trad4, _gamma_isrf, nh,
                            bi_itmask, sol, in, is, ie, my_rates);

            for (int i = is; i <= ie; i++) {
                if (bi_itmask[i]) {
                    if (sol[i] > 0.0) {
                        tdustnow[i] = bi_t_mid[i];
                    } else {
                        bi_t_high[i] = bi_t_mid[i];
                    }

                    if ((abs(bi_t_high[i] - tdustnow[i]) / tdustnow[i])
                                                            <= bi_tol) {
                        bi_itmask[i] = FALSE;
                        c_done ++;
                    }

                    // Check for all cells converged
                    if (c_done >= c_total) goto label_converged;
                }
            }
        }
        // If iteration count exceeded with bistection, end of the line.
        if (iter > itmax) {
#ifdef _OPENMP
#pragma omp critical
{
#endif
            fpritnf(stderr, 'calc_tdust_1d failed using both methods for %d cells.',
                        c_total - c_done);
#ifdef _OPENMP
}
#endif
        }
    }

    label_converged:

    // Copy values back to thrown slice
    for (int i = is; i <= ie; i++) {
        if (itmask[i]) {

            // Check for bad solutions
            if (tdustnow[i] < 0.0) {
#ifdef _OPENMP
#pragma omp critical
{
#endif
                fprintf(stderr, "calc_tdust_1d Newton method - T_dust < 0: i = %d\
                         j = %d k = %d nh = %g t_gas = %g t_rad = %g t_dust = %g",
                         i, j, k, nh[i], tgas[i], trad, tdustnow[i]);
#ifdef _OPENMP
}
#endif
            }

            tdust[i] = tdustnow[i];
        }
    }
    
    //End of function
    return SUCCESS;
}
