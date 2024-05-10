#include "simulation.h"
#include <iostream>
std::vector<double> simulation(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> theta, const double R)
{
    int N = x.size();
    std::vector<double> mean_theta(N, 0.0);
    double r_pow2 = R * R;
    for (int b = 0; b < N; ++b)
    {
        double sx = 0.0;
        double sy = 0.0;
        int count = 0;

        for (int i = 0; i < N; ++i)
        {
            if (i != b)
            {
                double distance_squared = pow(x[i] - x[b], 2) + pow(y[i] - y[b], 2);
                if (distance_squared < r_pow2)
                {
                    sx += cos(theta[i]);
                    sy += sin(theta[i]);
                    count++;
                }
            }
        }

        if (count > 0)
        {
            mean_theta[b] = atan2(sy, sx);
        }
    }

    return mean_theta;
}

std::vector<double> simulation_openmp(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> theta, const double R)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get each process rank

    int N = x.size();
    std::vector<double> mean_theta(N, 0.0);
    double r_pow2 = R * R;
#pragma omp parallel for
    for (int b = 0; b < N; ++b)
    {
        double sx = 0.0;
        double sy = 0.0;
        int count = 0;

        for (int i = 0; i < N; ++i)
        {
            if (i != b)
            {
                double distance_squared = pow(x[i] - x[b], 2) + pow(y[i] - y[b], 2);
                if (distance_squared < r_pow2)
                {
                    sx += cos(theta[i]);
                    sy += sin(theta[i]);
                    count++;
                }
            }
        }

        if (count > 0)
        {
            mean_theta[b] = atan2(sy, sx);
        }
    }

    return mean_theta;
}

std::vector<double> simulation_openmp_dy(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> theta, const double R)
{
    int N = x.size();
    std::vector<double> mean_theta(N, 0.0);
    double r_pow2 = R * R;
#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < N; ++b)
    {
        double sx = 0.0;
        double sy = 0.0;
        int count = 0;

        for (int i = 0; i < N; ++i)
        {
            if (i != b)
            {
                double distance_squared = pow(x[i] - x[b], 2) + pow(y[i] - y[b], 2);
                if (distance_squared < r_pow2)
                {
                    sx += cos(theta[i]);
                    sy += sin(theta[i]);
                    count++;
                }
            }
        }

        if (count > 0)
        {
            mean_theta[b] = atan2(sy, sx);
        }
    }

    return mean_theta;
}


std::vector<double> simulation_mpi(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> theta, const double R)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);  //Get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get each process rank




    int N = x.size();
    std::vector<double> mean_theta(N, 0.0);
    double r_pow2 = R * R;

    int chunk_size = N / size;
    int remainder = N % size;

    std::vector<int> recvcounts(size);
    for (int i = 0; i < size; ++i){
        recvcounts[i] = (i < remainder) ? (chunk_size + 1) : chunk_size;
    }

    std::vector<int> displs(size);
    displs[0] = 0;
    for (int i = 1; i < size; ++i){
        displs[i] = displs[i-1] + recvcounts[i-1];
    }

    chunk_size = recvcounts[rank];


    int start = displs[rank];
    int end = start + chunk_size;

    //std::cout << rank << ":" << chunk_size << "|" << start << "-" << end << std::endl;

    std::vector<double> local_mean_theta(chunk_size, 0.0);
    for (int b = start; b < end; ++b)
    {
        double sx = 0.0;
        double sy = 0.0;
        int count = 0;

        for (int i = 0; i < N; ++i)
        {
            if (i != b)
            {
                double distance_squared = pow(x[i] - x[b], 2) + pow(y[i] - y[b], 2);
                if (distance_squared < r_pow2)
                {
                    sx += cos(theta[i]);
                    sy += sin(theta[i]);
                    count++;
                }
            }
        }


        if (count > 0)
        {
            local_mean_theta[b-start] = atan2(sy, sx);
        }
    }

    //gather all results
    MPI_Gatherv(local_mean_theta.data(), chunk_size, MPI_DOUBLE, mean_theta.data(), recvcounts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //broadcast

    MPI_Bcast(mean_theta.data(), N , MPI_DOUBLE, 0 , MPI_COMM_WORLD);


    return mean_theta;    
}
