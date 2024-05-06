#include <iostream>
#include <assert.h>
#include <random>
#include <chrono>

#include "simulation.h"

const int Nt = 500; // number of time steps
const int N = 500;  // number of birds

void benchmark_simulation()
{
    // Finite Volume simulation

    // Simulation parameters
    double v0 = 1.0;  // Initial velocity
    double eta = 0.5; // random fluctuation in angle (in radians)
    double L = 10;    // size of the box
    double R = 1;     // interaction radius
    double dt = 0.2;  // time step

    // Initialise
    std::default_random_engine generator;
    generator.seed(static_cast<unsigned int>(17));
    std::uniform_real_distribution<double> uniform(0, 1);

    // Initialize bird position
    std::vector<double> x(N);
    std::vector<double> y(N);
    for (int i = 0; i < N; ++i)
    {
        x[i] = uniform(generator) * L;
        y[i] = uniform(generator) * L;
    }

    // Initialize bird velocity
    std::vector<double> theta(N);
    std::vector<double> vx(N);
    std::vector<double> vy(N);
    for (int i = 0; i < N; ++i)
    {
        theta[i] = 2 * M_PI * uniform(generator);
        vx[i] = v0 * cos(theta[i]);
        vy[i] = v0 * sin(theta[i]);
    }

    // start time counter
    auto start = std::chrono::high_resolution_clock::now();

    // Simulation
    for (int i = 0; i < Nt; i++)
    {
        for (int i = 0; i < N; ++i)
        {
            // Update position
            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;

            // Apply periodic boundary conditions
            x[i] = fmod(x[i], L);
            y[i] = fmod(y[i], L);
            if (x[i] < 0)
                x[i] += L; // Handle negative values
            if (y[i] < 0)
                y[i] += L; // Handle negative values
        }

        theta = simulation(x, y, theta, R);

        // Update velocities
        for (int i = 0; i < N; ++i)
        {
            // Update theta randomly
            theta[i] += eta * (uniform(generator) - 0.5);
            // Update velocity components
            vx[i] = v0 * cos(theta[i]);
            vy[i] = v0 * sin(theta[i]);
        }
    }
    // stop time counter

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Standard simulation: " << std::endl;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

void benchmark_simulation_openmp()
{ // Finite Volume simulation

    // Simulation parameters
    double v0 = 1.0;  // Initial velocity
    double eta = 0.5; // random fluctuation in angle (in radians)
    double L = 10;    // size of the box
    double R = 1;     // interaction radius
    double dt = 0.2;  // time step

    // Initialise
    std::default_random_engine generator;
    generator.seed(static_cast<unsigned int>(17));
    std::uniform_real_distribution<double> uniform(0, 1);

    // Initialize bird position
    std::vector<double> x(N);
    std::vector<double> y(N);
    for (int i = 0; i < N; ++i)
    {
        x[i] = uniform(generator) * L;
        y[i] = uniform(generator) * L;
    }

    // Initialize bird velocity
    std::vector<double> theta(N);
    std::vector<double> vx(N);
    std::vector<double> vy(N);
    for (int i = 0; i < N; ++i)
    {
        theta[i] = 2 * M_PI * uniform(generator);
        vx[i] = v0 * cos(theta[i]);
        vy[i] = v0 * sin(theta[i]);
    }

    // start time counter
    auto start = std::chrono::high_resolution_clock::now();

    // Simulation
    for (int i = 0; i < Nt; i++)
    {
        for (int i = 0; i < N; ++i)
        {
            // Update position
            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;

            // Apply periodic boundary conditions
            x[i] = fmod(x[i], L);
            y[i] = fmod(y[i], L);
            if (x[i] < 0)
                x[i] += L; // Handle negative values
            if (y[i] < 0)
                y[i] += L; // Handle negative values
        }

        theta = simulation_openmp(x, y, theta, R);

        // Update velocities
        for (int i = 0; i < N; ++i)
        {
            // Update theta randomly
            theta[i] += eta * (uniform(generator) - 0.5);
            // Update velocity components
            vx[i] = v0 * cos(theta[i]);
            vy[i] = v0 * sin(theta[i]);
        }
    }
    // stop time counter

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Openmp optimzation: " << std::endl;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

void benchmark_simulation_openmp_dy()
{ // Finite Volume simulation

    // Simulation parameters
    double v0 = 1.0;  // Initial velocity
    double eta = 0.5; // random fluctuation in angle (in radians)
    double L = 10;    // size of the box
    double R = 1;     // interaction radius
    double dt = 0.2;  // time step

    // Initialise
    std::default_random_engine generator;
    generator.seed(static_cast<unsigned int>(17));
    std::uniform_real_distribution<double> uniform(0, 1);

    // Initialize bird position
    std::vector<double> x(N);
    std::vector<double> y(N);
    for (int i = 0; i < N; ++i)
    {
        x[i] = uniform(generator) * L;
        y[i] = uniform(generator) * L;
    }

    // Initialize bird velocity
    std::vector<double> theta(N);
    std::vector<double> vx(N);
    std::vector<double> vy(N);
    for (int i = 0; i < N; ++i)
    {
        theta[i] = 2 * M_PI * uniform(generator);
        vx[i] = v0 * cos(theta[i]);
        vy[i] = v0 * sin(theta[i]);
    }

    // start time counter
    auto start = std::chrono::high_resolution_clock::now();

    // Simulation
    for (int i = 0; i < Nt; i++)
    {
        for (int i = 0; i < N; ++i)
        {
            // Update position
            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;

            // Apply periodic boundary conditions
            x[i] = fmod(x[i], L);
            y[i] = fmod(y[i], L);
            if (x[i] < 0)
                x[i] += L; // Handle negative values
            if (y[i] < 0)
                y[i] += L; // Handle negative values
        }

        theta = simulation_openmp_dy(x, y, theta, R);

        // Update velocities
        for (int i = 0; i < N; ++i)
        {
            // Update theta randomly
            theta[i] += eta * (uniform(generator) - 0.5);
            // Update velocity components
            vx[i] = v0 * cos(theta[i]);
            vy[i] = v0 * sin(theta[i]);
        }
    }
    // stop time counter

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Openmp dynamic optimzation: " << std::endl;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

void benchmark_simulation_mpi(int argc, char* argv[])
{
    // Finite Volume simulation

    // Simulation parameters
    double v0 = 1.0;  // Initial velocity
    double eta = 0.5; // random fluctuation in angle (in radians)
    double L = 10;    // size of the box
    double R = 1;     // interaction radius
    double dt = 0.2;  // time step

    // Initialise
    std::default_random_engine generator;
    generator.seed(static_cast<unsigned int>(17));
    std::uniform_real_distribution<double> uniform(0, 1);

    // Initialize bird position
    std::vector<double> x(N);
    std::vector<double> y(N);
    for (int i = 0; i < N; ++i)
    {
        x[i] = uniform(generator) * L;
        y[i] = uniform(generator) * L;
    }

    // Initialize bird velocity
    std::vector<double> theta(N);
    std::vector<double> vx(N);
    std::vector<double> vy(N);
    for (int i = 0; i < N; ++i)
    {
        theta[i] = 2 * M_PI * uniform(generator);
        vx[i] = v0 * cos(theta[i]);
        vy[i] = v0 * sin(theta[i]);
    }

    // start time counter
    auto start = std::chrono::high_resolution_clock::now();

    // Simulation
    for (int i = 0; i < Nt; i++)
    {
        for (int i = 0; i < N; ++i)
        {
            // Update position
            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;

            // Apply periodic boundary conditions
            x[i] = fmod(x[i], L);
            y[i] = fmod(y[i], L);
            if (x[i] < 0)
                x[i] += L; // Handle negative values
            if (y[i] < 0)
                y[i] += L; // Handle negative values
        }

        theta = simulation_mpi(x, y, theta, R, argc, argv);

        // Update velocities
        for (int i = 0; i < N; ++i)
        {
            // Update theta randomly
            theta[i] += eta * (uniform(generator) - 0.5);
            // Update velocity components
            vx[i] = v0 * cos(theta[i]);
            vy[i] = v0 * sin(theta[i]);
        }
    }
    // stop time counter

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "mpi simulation: " << std::endl;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

int main(int argc, char* argv[])
{
    std::cout << "Running benchmarks..." << std::endl;
    //benchmark_simulation();
    //benchmark_simulation_openmp();
    //benchmark_simulation_openmp_dy();
    //benchmark_simulation_mpi(argc, argv);

    std::cout
        << "All tests passed!" << std::endl;


    
    int rank, size, provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);  //Get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get each process rank

    std::cout << size;

    int N = x.size();
    std::vector<double> mean_theta(N, 0.0);
    double r_pow2 = R * R;

    int chunk_size = N/ size;

    int start = rank*chunk_size;
    int end = start + chunk_size;

    std::vector<double> local_mean_theta(chunk_size, 0.0);


    //super simple error-handling
    if(end > N){
        exit(1);
    }

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
    MPI_Gather(local_mean_theta.data(), chunk_size, MPI_DOUBLE, mean_theta.data(), chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Finalize();
}