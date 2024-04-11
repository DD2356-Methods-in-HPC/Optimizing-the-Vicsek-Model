#include <iostream>
#include <random>
#include "simulation.h"
#include <vector>

int main()
{
    // Finite Volume simulation

    // Simulation parameters
    double v0 = 1.0;  // Initial velocity
    double eta = 0.5; // random fluctuation in angle (in radians)
    double L = 10;    // size of the box
    double R = 1;     // interaction radius
    double dt = 0.2;  // time step
    int Nt = 200;     // number of time steps
    int N = 500;      // number of birds

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

    // plot
    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (gnuplotPipe)
    {
        // Configurar el gráfico
        fprintf(gnuplotPipe, "set terminal png\n");
        fprintf(gnuplotPipe, "set output 'assets/birds.png'\n");
        fprintf(gnuplotPipe, "set xrange [0:%f]\n", L);
        fprintf(gnuplotPipe, "set yrange [0:%f]\n", L);
        fprintf(gnuplotPipe, "unset border\n");
        fprintf(gnuplotPipe, "unset xtics\n");
        fprintf(gnuplotPipe, "unset ytics\n");
        fprintf(gnuplotPipe, "set size square\n");

        // Dibujar los vectores
        for (int i = 0; i < N; ++i)
        {
            fprintf(gnuplotPipe, "set arrow %d from %f,%f to %f,%f lt 1 lw 0 filled\n", i + 1, x[i], y[i], x[i] + vx[i], y[i] + vy[i]);
        }
        fprintf(gnuplotPipe, "plot NaN title ''\n"); // Esto es necesario para que se muestren las flechas
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
        std::cout << "Gráfico guardado como 'birds.png'" << std::endl;
    }
    else
    {
        std::cerr << "Error: No se pudo abrir la tubería a Gnuplot." << std::endl;
        return 1;
    }
    return 0;
}