#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Declare global variables
double **F; // F matrix
double **U; // U matrix
double **U_old; // Temporary U matrix

//Declare blocks for continuous memory allocation
double *F_block;
double *U_block;
double *U_old_block;

int N; // Grid size
double Delta; // Grid spacing
int k_max; // Number of iterations
double threshold; // Threshold limit

void init()
{
    // Variables used for loops
    int i, j;

    // Calculate grid spacing
    Delta = 2.0 / ((double) N - 1.0);
    // Make room for boundaries
    N += 2;

    // Allocate memory for arrays in the heap
    F = (double **) malloc(N * sizeof(double **));
    U = (double **) malloc(N * sizeof(double **));
    U_old = (double **) malloc(N * sizeof(double **));
    F_block = (double *) calloc(N * N, sizeof(double));
    U_block = (double *) calloc(N * N, sizeof(double));
    U_old_block = (double *) calloc(N * N, sizeof(double));

    for (i = 0; i < N; i++)
    {
        F[i] = &F_block[i * N];
        U[i] = &U_block[i * N];
        U_old[i] = &U_old_block[i * N];
    }

    // Declare relative coordinates
    double x = -1.0;
    double y = -1.0;
    double x_lower = 0.0;
    double x_upper = 1.0 / 3.0;
    double y_lower = -2.0 / 3.0;
    double y_upper = -1.0 / 3.0;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            // Place radiator for F in the right place
            if (x <= x_upper && x >= x_lower && y <= y_upper && y >= y_lower)
            {
                // Set radiator value to 200 degrees
                F[i][j] = 200.0;
            }
            // Place temperature for walls
            if (i == (N - 1) || i == 0 || j == (N - 1))
            {
                // Set temperature to 20 degrees for 3 of the walls
                U[i][j]     = 20.0;
                U_old[i][j] = 20.0;
            }
            // Move relative coordinates by one unit of grid spacing
            y += Delta;
        }
        // Move relative coordinates by one unit of grid spacing
        x += Delta;
        y = -1.0;
    }
}

void deinit()
{
    free(F);
    free(U);
    free(U_old);
}

void print_matrix(double **M)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            // Swap indecies to show correct x and y-axes
            printf("%.2f\t", M[j][i]);
        }
        printf("\n");
    }
}

void jacobi()
{
    // Declare variables
    double one_fourth = 0.25;
    double d = 0;
    int k = 0, i, j;
    double **temp; // Create temporary pointer
    double old_val;

    // Run do-while loop
    do
    {
        // Run through entire matrix
        for (i = 1; i < N - 1; i++)
        {
            for (j = 1; j < N - 1; j++)
            {
                // Save old value
                old_val = U_old[i][j];
                // Calculate new value from surrounding points
                U_old[i][j] = (U[i][j - 1] + U[i][j + 1] + U[i - 1][j] + U[i + 1][j] + (Delta * Delta * F[i][j])) * one_fourth;
                // Calculate difference between old and new value
                old_val -= U_old[i][j];
                // Take absolute value of difference
                d += (old_val < 0 ? -old_val : old_val);
            }
        }
        // Take average of difference between old and new values
        d /= ((double) N) * ((double) N);

        // Swap pointers for U and U_old
        temp = U;
        U = U_old;
        U_old = temp;
        k += 1;
    }
    while (k < k_max && d > threshold);
}

void gs()
{
    // Declare variables
    double one_fourth = 0.25;
    double d = 0;
    int k = 0, i, j;
    double old_val;

    // Run do-while loop
    do
    {
        // Run through entire matrix
        for (i = 1; i < N - 1; i++)
        {
            for (j = 1; j < N - 1; j++)
            {
                // Save old value
                old_val = U[i][j];
                // Calculate new value from surrounding points
                U[i][j] = (U[i][j - 1] + U[i][j + 1] + U[i - 1][j] + U[i + 1][j] + (Delta * Delta * F[i][j])) * one_fourth;
                // Calculate difference between old and new value
                old_val -= U[i][j];
                // Take absolute value of difference
                d += (old_val < 0 ? -old_val : old_val);
            }
        }
        // Take average of difference between old and new values
        d /= ((double) N) * ((double) N);

        k += 1;
    }
    while (k < k_max && d > threshold);
}

void omp_jacobi()
{
    // Declare variables
    double one_fourth = 0.25;
    double d = 0;
    int k = 0, i, j;
    double **temp; // Create temporary pointer
    double old_val;

    // Run do-while loop
    do
    {
        // Run through entire matrix
        #pragma omp parallel
        {
            #pragma omp for private(i,j,old_val) reduction(+:d)
            for (i = 1; i < N - 1; i++)
            {
                for (j = 1; j < N - 1; j++)
                {
                    // Save old value
                    old_val = U_old[i][j];
                    // Calculate new value from surrounding points
                    U_old[i][j] = (U[i][j - 1] + U[i][j + 1] + U[i - 1][j] + U[i + 1][j] + (Delta * Delta * F[i][j])) * one_fourth;
                    // Calculate difference between old and new value
                    old_val -= U_old[i][j];
                    // Take absolute value of difference
                    d += (old_val < 0 ? -old_val : old_val);
                }
            }
        } /* end omp parallel */
        // Take average of difference between old and new values
        d /= ((double) N) * ((double) N);

        // Swap pointers for U and U_old
        temp = U;
        U = U_old;
        U_old = temp;
        k += 1;
    }
    while (k < k_max && d > threshold);
}

int main(int argc, char *argv[])
{
    // Set default values
    N = 30;           // Grid size
    k_max = 1000;     // Number of iterations
    threshold = 0.01; // Threshold limit

    // Set values according to program arguments
    if (argc > 2)
    {
        N = atoi(argv[2]);
    }
    if (argc > 3)
    {
        k_max = atoi(argv[3]);
    }
    if (argc > 4)
    {
        threshold = atof(argv[4]);
    }

    // Initialize values
    init();

    // Run algorithm specified by program argument
    if (!strcmp(argv[1], "jac"))
    {
        jacobi();
    }
    else if (!strcmp(argv[1], "gs"))
    {
        gs();
    }
    else if (!strcmp(argv[1], "omp"))
    {
        omp_jacobi();
    }
    else
    {
        printf("Did not give a method as argument, try adding 'jac', 'gs', 'omp'\n");
    }
    print_matrix(U);

    // Deinitialize to prevent memory leaks
    deinit();
    return 0;
}
