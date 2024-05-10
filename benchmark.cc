#include <iostream>
#include <assert.h>
#include <random>
#include <chrono>
#include <functional>

#include "simulation.h" 

const int Nt = 500; // number of time steps
const int N = 500;  // number of birds



void benchmark_simulation(int rank, std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&, std::vector<double>, double)> sim_function, std::string func_string)
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

        theta = sim_function(x, y, theta, R);

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

    if(rank == 0){
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << func_string << std::endl;
        std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    }
}

int main(int argc, char* argv[])
{

    int provided, rank;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get each process rank

    if(rank==0) //only print from root process
        std::cout << "Running benchmarks..." << std::endl;


    //Run non-MPI benchmarks only on root
    if(rank == 0){
        benchmark_simulation(rank, simulation, "simulation");
        benchmark_simulation(rank, simulation_openmp, "simulation_openmp");
        benchmark_simulation(rank, simulation_openmp_dy, "simulation_openmp_dy");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    benchmark_simulation(rank, simulation_mpi, "simulation_mpi");

    MPI_Finalize();
}