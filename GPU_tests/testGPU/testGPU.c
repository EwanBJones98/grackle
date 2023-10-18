#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_thread_num() {return 0;}
static int omp_get_num_threads() {return 1;}
static double omp_get_wtime() {return (clock() / CLOCKS_PER_SEC);}
#endif

// Function to write timing statistics to file.
void write_timings(int matrixSidelength, double timetaken, int numiter)
{
    FILE *fptr;
    int write_header = 0;
    // See if the file already exists, if not then create it and write the table header.
    if ((fptr = fopen("timings.txt", "r")) == NULL){
        write_header = 1;
    }
    // Open file to write results to.
    if ((fptr = fopen("timings.txt", "a")) == NULL){
        fprintf(stdout, "Error opening file!\n");
        exit(1);
    }
    // Write header if file was just created.
    if (write_header){
        fprintf(fptr, "MatrixSideLength \t NumIters \t NumThreads \t OpenMP \t TotalTime \t AverageTime \n\n\n\n");
    }

    double averageTimetaken;
    // Calculate average, max and min x of array.
    // Complete calculation of average x.
    averageTimetaken = timetaken / numiter;

    // Find OpenMP mode being used.
    char *omp_mode = malloc(20 * sizeof(char));
    #ifdef GPU
            omp_mode = "GPU";
    #elif CPU
        omp_mode = "CPU";
    #elif SERIAL
        omp_mode = "off";
    #endif

    // Get number of threads being used.
    int num_threads;
    #ifdef CPU
        #pragma omp parallel
        {
            #pragma omp single
            {
                num_threads = omp_get_num_threads();
            }
        }
    #elif GPU
        #pragma omp target teams map(tofrom:num_threads)
        {
            #pragma omp parallel
            #pragma omp single
            num_threads = omp_get_num_threads() * omp_get_num_teams();
        }
    #elif SERIAL
        num_threads = 1;
    #endif

    // Write results to the file.
    fprintf(fptr, "%d | %d | %d | %s | %g | %g\n", matrixSidelength, numiter, num_threads, omp_mode, timetaken, averageTimetaken);
    fclose(fptr);
}


int main(int argc, char *argv[])
{
    // Number of iterations to time
    int numiter = 1000000000;

    // Allocate space for matrices
    int N=10000;
    double* x = (double*) malloc(N * N * sizeof(double));
    double* y = (double*) malloc(N * N * sizeof(double));
    double* z = (double*) malloc(N * N * sizeof(double));

    // Initialise matrices x and y
    for (int i = 0; i < N; i++){
        for (int j=0; j < N; j++){
            x[(i*N) + j] = i;
            y[(i*N) + j] = i*2;
        }
    }

    // Time 1000 iterations of the matrices x and y
    double start_time = omp_get_wtime();
    #ifdef GPU
        #pragma omp target enter data map(to:x[0:N*N],y[0:N*N],z[0:N*N],N,numiter)
        #pragma omp target teams distribute parallel for collapse(4)
    #elif CPU
        #pragma omp parallel for default (none) shared(x,y,z,N,numiter) collapse(4)
    #endif
    for (int iter=0; iter < numiter; iter++){
        for (int i = 0; i < N; i++){
            for (int j=0; j < N; j++){
                for (int k=0; k < N; k++){
                    z[(i*N) + j] += x[(i*N) + k] * y[(k*N) + j];
                }
            }
        }
    }
    #ifdef GPU
        #pragma omp target exit data map(from:z[0:N*N])
    #endif
    double elapsed_time = omp_get_wtime() - start_time;

    // Write results to file
    write_timings(N, elapsed_time, numiter);

    return 1;
}