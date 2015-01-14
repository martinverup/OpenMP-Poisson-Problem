#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENPM
#include <omp.h>
#endif

void init(int n, double delta, double **f, double **u, double **u_out)
{

    double x = -1.0;
    double y = -1.0;
    int i, j;
    double x_lower = 0.0;
    double x_upper = 1.0 / 3.0;
    double y_lower = -2.0 / 3.0;
    double y_upper = -1.0 / 3.0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            f[i][j] = 0.0;
            u[i][j] = 0.0;
            u_out[i][j] = 0.0;
            if (x <= x_upper && x >= x_lower && y <= y_upper && y >= y_lower)
            {
                f[i][j] = 200.0;
            }

            if (i == (n - 1) || i == 0 || j == (n - 1))
            {
                u[i][j] = 20.0;
                u_out[i][j] = 20.0;
            }
            y += delta;
        }
        x += delta;
        y = -1.0;
    }
}

void print_matrix(int n, double **f)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%.1f\t", f[j][i]);
        }
        printf("\n");
    }
}

double one_jacobi_step(int n, double one_fourth, double delta2, double **u, double **f, double **u_new)
{
    double sum = 0;
    int i, j;
    for (i = 1; i < n - 1; i++)
    {
        for (j = 1; j < n - 1; j++)
        {
            double temp = u_new[i][j];
            u_new[i][j] = (u[i][j - 1] + u[i][j + 1] + u[i - 1][j] + u[i + 1][j] + (delta2 * f[i][j])) * one_fourth;

            temp -= u_new[i][j];
            sum += (temp < 0 ? -temp : temp);
        }
    }
    return sum /= ((double) n) * ((double) n);
}

void jacobi(int n, int k_max, double delta2, double **u, double **f, double **u_out, double threshold)
{

    double one_fourth = 1.0 / 4.0, d = 0;
    int k = 0;
    do   //running k times
    {
        d = one_jacobi_step(n, one_fourth, delta2, u, f, u_out);
        //copying back to the initializing array
        double **temp = u;
        u = u_out;
        u_out = temp;
        k += 1;
    }
    while (k < k_max && d > threshold);
    printf("number of iterations where: %d",k);
}

double one_gauss_seidel(int n, double one_fourth, double delta2, double **u, double **f)
{
    int i, j;
    double temp, change_temp;
    double change_average = 0;
    for (i = 1; i < n - 1; i++)
    {
        for (j = 1; j < n - 1; j++)
        {
            temp = u[i][j];
            u[i][j] = one_fourth * (u[i][j - 1] + u[i][j + 1] + u[i - 1][j] + u[i + 1][j] + delta2 * f[i][j]);

            //calculating the average of the change by first summing up the absolute value of the changes
            change_temp = u[i][j] - temp;
            if (change_temp < 0.0)
            {
                change_temp = 0 - change_temp;
            }
            change_average += change_temp;
        }
    }
    return change_average / ((double) n * (double) n);
}

void gauss_seidel(int n, int k, double delta2, double **u, double **f, double threshold)
{
    double one_fourth = 1.0 / 4.0;
    double d = 0.0;
    int h = 0;
    do
    {
        d = one_gauss_seidel(n, one_fourth, delta2, u, f);
        h += 1;
    }
    while (d > threshold && h < k);
}
double one_omp_jacobi_step(int n, double one_fourth, double delta2, double **u, double **f, double **u_new) {
    double sum = 0;
    int i, j;
    double temp;
    #pragma omp parallel
    {
        #pragma omp for private(i,j,temp) reduction(+:sum)
        for (i = 1; i < n - 1; i++)
        {
            for (j = 1; j < n - 1; j++)
            {
                temp = u_new[i][j];
                u_new[i][j] = (u[i][j - 1] + u[i][j + 1] + u[i - 1][j] + u[i + 1][j] + (delta2 * f[i][j])) * one_fourth;

                temp -= u_new[i][j];
                sum += (temp < 0 ? -temp : temp);
            }
        }
    }
    return sum /= ((double) n) * ((double) n);
}

void omp_jacobi(int n, int k_max, double delta2, double **u, double **f, double **u_out, double threshold) {
    double one_fourth = 1.0 / 4.0, d = 0;
    int k = 0;
    do   //running k times
    {
        d = one_omp_jacobi_step(n, one_fourth, delta2, u, f, u_out);
        //copying back to the initializing array
        double **temp = u;
        u = u_out;
        u_out = temp;
        k += 1;
    }
    while (k < k_max && d > threshold);
    printf("number of iterations where: %d",k);
}


int main(int argc, char *argv[])
{

    int N = 20;
    int iter = 10000;
    double threshold = 0.01; 
    if (argc > 1)
    {
        N = atoi(argv[1]);
    }
    if (argc > 2)
    {
        iter = atoi(argv[2]);
    }
    if (argc > 3)
    {
        threshold = atof(argv[3]);
    }

    double delta = 2.0 / ((double)N - 1.0);
    N += 2;
    double delta2 = delta * delta;
    double (**f) = malloc(sizeof(*f) * N);
    double (**u) = malloc(sizeof(*u) * N);
    double (**u_out) = malloc(sizeof(*u_out) * N);
    int i;
    for (i = 0; i < N; i++)
    {
        f[i] = malloc(sizeof(*f[i]) * N);
        u[i] = malloc(sizeof(*u[i]) * N);
        u_out[i] = malloc(sizeof(*u_out[i]) * N);
    }
    init(N, delta, f, u, u_out);
    //print_matrix(N, f);
    if (argc > 4)
    {  
        if(!strcmp(argv[4], "jac")) {
            jacobi(N, iter, delta2, u, f, u_out, threshold);
            //print_matrix(N, u_out);
        }else if(!strcmp(argv[4], "gs")) {
            gauss_seidel(N, iter, delta2, u, f, threshold);
            //print_matrix(N, u);
        }else if(!strcmp(argv[4], "omp")) {
            omp_jacobi(N, iter, delta2, u, f, u_out, threshold);
            //print_matrix(N, u_out);
        }else {
            printf("Did not give a method as argument, try adding 'jac', 'gs', 'omp'\n");
        }
    }

    free(f);
    free(u);
    free(u_out);
    return 0;
}
