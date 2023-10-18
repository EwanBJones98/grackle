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

// ####################################
//
// Author: Ewan Jones
//
// This code is just to test if I can get GPU parallelism up and working on Cuillin.
// The code itself will compute the inner product of two NxN matrices on the GPU.
//
// ####################################
    

// Function to multiply two NxN matrices, A and B, together to get result C.
// Then set A=C, B=A and compute C again. Do this num_iter, times.
double multiply_matrix(int *A, int *B, int N, int num_iter)
{
    // Create output matrix -- all elements set to zero.
    int* C = (int*) malloc(N * N * sizeof(int));
    for (int i = 0; i < N*N; i++){
        C[i] = 0;
    }
    
    // Populate output matrix elements as dot-product of both input matrices.
    double time_start = omp_get_wtime();
    #ifdef OMP_GPU
        #pragma omp target enter data map(to:A[0,N*N],B[0,N*N],C[0,N*N],N)
    #endif
    #pragma omp target
    {
    for (int iter = 0; iter < num_iter; iter++){
        // Multiply A and B to get C.
        #ifdef OMP_GPU
            #pragma omp teams distribute parallel for collapse(2)
        #elif OMP_CPU
            #pragma omp parallel for default(none) shared(A,B,C,N) collapse(2)
        #endif
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                for (int k = 0; k < N; k++){
                    C[(i*N) + j] += A[(i*N) + k] * B[(k*N) + j];
                }
            }
        }
    }
    }
    #ifdef OMP_GPU
        #pragma omp target exit data map(from:C[0,N*N])
    #endif
    double time_total = omp_get_wtime() - time_start;

    // Return timings
    return time_total;
}

// Function to write stats from given array to a file for analysis.
// Currently finds average, max and min.
void write_stats_from_array(double x, int matrix_size, int num_iter)
{
    FILE *fptr;
    int write_header = 0;
    // See if the file already exists, if not then create it and write the table header.
    if ((fptr = fopen("parallel_performance.txt", "r")) == NULL){
        write_header = 1;
    }
    // Open file to write results to.
    if ((fptr = fopen("parallel_performance.txt", "a")) == NULL){
        fprintf(stdout, "Error opening file!\n");
        exit(1);
    }
    // Write header if file was just created.
    if (write_header){
        fprintf(fptr, "MatrixSize \t NumIters \t NumThreads \t OpenMP \t TotalTime \t AverageSingleTime \n\n\n\n");
    }

    double average_x;
    // Calculate average, max and min x of array.
    // Complete calculation of average x.
    average_x = x / num_iter;

    // Find OpenMP mode being used.
    char *omp_mode = malloc(20 * sizeof(char));
    #ifdef OMP_GPU
            omp_mode = "GPU";
    #elif OMP_CPU
        omp_mode = "CPU";
    #else
        omp_mode = "off";
    #endif

    // Get number of threads being used.
    int num_threads;
    #ifdef OMP_CPU
        #pragma omp parallel
        {
            #pragma omp single
            {
                num_threads = omp_get_num_threads();
            }
        }
    #elif OMP_GPU
        #pragma omp target teams map(tofrom:num_threads)
        {
            #pragma omp parallel
            #pragma omp single
            num_threads = omp_get_num_threads() * omp_get_num_teams();
        }
    #endif

    // Write results to the table.
    fprintf(fptr, "%d \t %d \t %d \t %s \t %g \t %g \n", matrix_size, num_iter, num_threads, omp_mode, x, average_x);
}

int main(int argc, char *argv[])
{
    // Print message indicating parallelisation scheme.
    #ifdef OMP_GPU
        fprintf(stdout, "Running with GPU offload.\n");
        // Check gpu is activated
        int gpu_active;
        #pragma omp target map(tofrom:gpu_active)
        {
            gpu_active = omp_is_initial_device();
        }
        fprintf(stdout, "GPU active = %d, (0 means gpu IS being used)\n", gpu_active);
    #elif OMP_CPU
        fprintf(stdout, "Running with CPU threading.\n");
    #else
        fprintf(stdout, "Running serial code.\n");
    #endif

    // Get command line arguments from the user.
    int matrix_size, num_iter;
    #ifdef _OPENMP
        // Must pass three arguments in from command line -- matrix size, number of iterations and number of threads.
        if (argc != 4){
            fprintf(stderr, "You must pass the matrix size, number of iterations and number of threads. Only %d arguments supplied!", argc);
        return 0;
        }
        omp_set_num_threads(atoi(argv[3]));
        //#ifdef OMP_CPU
            //omp_set_num_threads(atoi(argv[3]));
        //#endif
    #else
        // Must pass two arguments in from command line, matrix size and number of iterations.
        if (argc != 3){
            fprintf(stderr, "You must pass the matrix size and number of iterations. Only %d arguments supplied!", argc);
            return 0;
        }
        
    #endif
    matrix_size = atoi(argv[1]);
    num_iter = atoi(argv[2]);
    fprintf(stdout, "-- finished reading command line arguments --\n");

    // Run the matrix multiplication code many times to get many timings.
    int *matrix_1 = (int*) malloc(matrix_size * matrix_size * sizeof(int));
    int *matrix_2 = (int*) malloc(matrix_size * matrix_size * sizeof(int));
    for (int i = 0; i < matrix_size*matrix_size; i++){
        matrix_1[i] = i*i;
        matrix_2[i] = i*7;
    }
    fprintf(stdout, "-- finished creating initial matrices --\n");

    // Calculate the dot product of both matrices and time the calculation.
    double totalTime = multiply_matrix(matrix_1, matrix_2, matrix_size, num_iter);

    // Print min, max and average of timings.
    write_stats_from_array(totalTime, matrix_size, num_iter);

    return 1;
}