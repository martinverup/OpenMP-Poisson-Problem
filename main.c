#include <stdio.h>
#include <stdlib.h>

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

void one_jacobi_step(int n, double one_fourth, double delta2, double **u, double **f, double **u_new)
{
    int i, j;
    for (i = 1; i < n - 1; i++)
    {
        for (j = 1; j < n - 1; j++)
        {
            u_new[i][j] = (u[i][j - 1] + u[i][j + 1] + u[i - 1][j] + u[i + 1][j] + (delta2 * f[i][j])) * one_fourth;
        }
    }
}

void jacobi(int n, int k, double delta2, double **u, double **f, double **u_out)
{
    
    double one_fourth = 1.0 / 4.0;
    int h, i, j;
    for (h = 0; h < k; h++)   //running k times
    {
        one_jacobi_step(n, one_fourth, delta2, u, f, u_out);
        //copying back to the initializing array
        double **temp = u;
        u = u_out;
        u_out = temp;
    }
}

int main(int argc, char *argv[])
{

    int N = 20;
    int iter = 10000;
    if (argc > 1)
    {
        N = atoi(argv[1]);
    }
    if (argc > 2)
    {
        iter = atoi(argv[2]);
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

    print_matrix(N, f);
    printf("\n");
    jacobi(N, iter, delta2, u, f, u_out);
    print_matrix(N, u_out);
    free(f);
    free(u);
    free(u_out);

    return 0;
}