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
// Date:   August 2022
//
// This code is just to test if I can get GPU parallelism up and working on Cuillin.
// The code itself will compute the inner product of two NxN matrices on the GPU.
//
// ####################################
    
// Function to populate NxN matrix with random numbers in the interval [0,X], where X is integer.
int *rand_matrix(int N, int X)
{
    // Create empty NxN array.
    int* A = (int*) malloc(N * N * sizeof(int));
    
    // Populate array with random numbers.
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            A[(i*N) + j] = rand() % X;
        }
    }

    // Return pointer to the array.
    return A;   
}

// Function to multiply two NxN matrices together.
int *multiply_matrix(int *A, int *B, int N)
{
    // Create output matrix -- all elements set to zero.
    int* C = (int*) malloc(N * N * sizeof(int));
    for (int i = 0; i < N*N; i++){
        C[i] = 0;
    }
   
    // Populate output matrix elements as dot-product of both input matrices.
    #ifdef OMP_GPU
        #pragma omp target teams distribute parallel for collapse(2) map(to:A,B,N), map(tofrom:C)
    #elif OMP_CPU
        #pragma omp parallel for default(none) shared(A,B,C,N) collapse(2)
    #endif
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                C[(i*N) + j] += A[(i*N) + j] * B[(i*N) + j];
            }
        }
    }

    // Return output matrix.
    return C;
}

// Function to print a matrix to the terminal such that it is human-readable.
void print_matrix(int* A, int N)
{
    
    // Print matrix elements separated by an indent.
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            fprintf(stdout, "%d\t", A[(i*N) + j]);
        }
        fprintf(stdout, "\n");
    }
}

// Function which creates, multiplies and prints to terminal the calculation of the dot-product of two square matrices.
// The function will return the time taken to multiply the matrices together.
double dot_rand_matrix(int matrix_size, int max_rand)
{

    // Create two NxN matrices populated with random integers.
    int *matrix_1 = rand_matrix(matrix_size, max_rand);
    int *matrix_2 = rand_matrix(matrix_size, max_rand);

    double time_start, time_end, time_elapsed;
    // Begin timer.
    time_start = omp_get_wtime();

    // Calculate the dot product of both matrices and time the calculation.
    int *matrix_3 = multiply_matrix(matrix_1, matrix_2, matrix_size);

    // Calculate elapsed time.
    time_end = omp_get_wtime();
    time_elapsed = time_end - time_start;

   return time_elapsed;
}

// Function to write stats from given array to a file for analysis.
// Currently finds average, max and min.
void write_stats_from_array(double* arr, int arr_len, int matrix_size, int num_iter)
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
        fprintf(fptr, "MatrixSize \t NumIters \t NumThreads \t OpenMP \t MinTime \t MaxTime \t MeanTime \n\n\n\n");
    }

    double average_x, max_x=0, min_x=0;
    // Calculate average, max and min x of array.
    for (int i = 0; i < arr_len; i++){
        double x;
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

    // Find OpenMP mode being used.
    char *omp_mode = malloc(20 * sizeof(char));
    #ifdef _OPENMP
        #ifdef GPU
            omp_mode = "GPU";
        #endif
        omp_mode = "CPU";
    #else
        omp_mode = "off";
    #endif

    // Get number of threads being used.
    int num_threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
        }
    }

    // Write results to the table.
    fprintf(fptr, "%d \t %d \t %d \t %s \t %g \t %g \t %g \n", matrix_size, num_iter, num_threads, omp_mode, min_x, max_x, average_x);
}

int main(int argc, char *argv[])
{
    // Print message indicating parallelisation scheme.
    #ifdef OMP_GPU
        fprintf(stdout, "Running with GPU offload.\n");
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
    #else
        // Must pass two arguments in from command line, matrix size and maximum random number generated.
        if (argc != 3){
            fprintf(stderr, "You must pass the matrix size and number of iterations. Only %d arguments supplied!", argc);
            return 0;
        }
        
    #endif
    matrix_size = atoi(argv[1]);
    num_iter = atoi(argv[2]);
    

    // Initialisation for random number generation.
    int max_rand = 10; 
    srand(time(NULL));

    // Run the matrix multiplication code 100 times to get many timings.
    double times[num_iter];
    int ten_percent = (int) num_iter/10;
    for (int i = 0; i < num_iter; i++){
        times[i] = dot_rand_matrix(matrix_size, max_rand);
        if (i % ten_percent == 0){
            fprintf(stderr, "%d%% done...\n", i*10/ten_percent);
        }
    }
    fprintf(stderr, "100%% done...\n");

    // Print min, max and average of timings.
    write_stats_from_array(times, num_iter, matrix_size, num_iter);

    return 1;
}