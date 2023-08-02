#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
  
int main () {

  #ifndef _OPENMP
    fprintf(stdout, "Not compiled with OpenMP! Quitting...");
    exit(1);
  #endif

  // Initialise values of x and y
  double x, y;
  x=0;

  // Check gpu is activated
  int gpu_active;
  #pragma omp target map(tofrom:gpu_active)
  {
    gpu_active = omp_is_initial_device();
  }
  fprintf(stdout, "GPU active = %d, (0 means gpu IS being used)\n", gpu_active);

  // Print number of threads being used on GPU
  int nteams, nthreads; 
  #pragma omp target teams map(tofrom:nteams,nthreads)
  {
    nteams = omp_get_num_teams();
    #pragma omp parallel
    #pragma omp single
    nthreads = omp_get_num_threads();
  }
  fprintf(stdout, "Number of teams = %d.\nNumber of threads = %d.\n", nteams, nthreads);

  double t1, t2, elapsed_time;
  t1 = omp_get_wtime();
  #pragma omp target teams distribute parallel for map(tofrom:x) reduction(+:x)
  for (int i=0; i<1000000; i++) {
    y = 0.5;
    for (int j=0; j<1000; j++) {
      y=(y+sqrtf(y)/2);
      x=x+1-y;
    }
  }
  t2 = omp_get_wtime();
  elapsed_time = t2 - t1;
  fprintf(stdout, "Our final result is: %f\n", x);
  fprintf(stdout, "Computation took %g seconds!\n",elapsed_time);
  return 0;
}
