#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

typedef struct
{
    double a;
    double b;
    double tol;
    double whole;
} Task;

#define TAG_WORK 1
#define TAG_RESULT 2
#define TAG_NEW_TASK 3
#define TAG_STOP 4

typedef double (*func_ptr)(double);

double fx0(const double x) { return sin(x) + 0.5 * (cos(3 * x)); }
double fx1(const double x) { return 1.0 / (1.0 + 100.0 * pow((x - 0.3), 2)); }
double fx2(const double x) { return sin(200.0 * x) * exp(-x); }

double adaptive_simpson(func_ptr f, double a, double b, double tol, double whole)
{
    double m = (a + b) / 2.0;
    double left_simpson = (m - a) / 6.0 * (f(a) + 4.0 * f((a + m) / 2.0) + f(m));
    double right_simpson = (b - m) / 6.0 * (f(m) + 4.0 * f((m + b) / 2.0) + f(b));

    if (fabs(left_simpson + right_simpson - whole) <= 15.0 * tol)
    {
        return left_simpson + right_simpson + (left_simpson + right_simpson - whole) / 15.0;
    }

    return adaptive_simpson(f, a, m, tol / 2.0, left_simpson) +
           adaptive_simpson(f, m, b, tol / 2.0, right_simpson);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 4)
    {
        if (rank == 0)
            printf("Usage: mpirun -np P ./integration {func_id} {mode} {tol}\n");
        MPI_Finalize();
        return 1;
    }

    const int func_id = atoi(argv[1]);
    const int mode = atoi(argv[2]);
    const double tol = atof(argv[3]);

    func_ptr f;

    switch (func_id)
    {
    case 0:
        f = fx0;
        break;
    case 1:
        f = fx1;
        break;
    case 2:
        f = fx2;
        break;
    default:
        if (rank == 0)
            printf("Invalid function ID: %d\n", func_id);
        MPI_Finalize();
        return 1;
    }

    if (mode == 0)
    {
        if (rank == 0)
        {
            double start_time = MPI_Wtime();

            double initial_estimate = (1.0 - 0.0) / 6.0 * (f(0.0) + 4.0 * f(0.5) + f(1.0));
            double result = adaptive_simpson(f, 0.0, 1.0, tol, initial_estimate);

            double end_time = MPI_Wtime();
            printf("Serial Result: %e\n", result);
            printf("Time Taken: %f seconds\n", end_time - start_time);
        }
    }

    if (mode == 1)
    {
        double start_time = MPI_Wtime();

        int K = size;
        double h = 1.0 / K;
        double local_sum = 0.0;

        for (int i = rank; i < K; i += size)
        {
            double a = i * h;
            double b = (i + 1) * h;
            double mid = (a + b) / 2.0;
            double initial_estimate = (h / 6.0) * (f(a) + 4.0 * f(mid) + f(b));

            local_sum += adaptive_simpson(f, a, b, tol / K, initial_estimate);
        }

        double global_sum = 0.0;
        MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        double end_time = MPI_Wtime();

        if (rank == 0)
        {
            printf("Mode 1 Result: %e\n", global_sum);
            printf("Mode 1 Time: %fs\n", end_time - start_time);
            printf("K value used: %d\n", K);
        }
    }

    MPI_Finalize();
    return 0;
}