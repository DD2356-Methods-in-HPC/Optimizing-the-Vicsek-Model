#include <vector>
#include <cmath>

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
    int N = x.size();
    std::vector<double> mean_theta(N, 0.0);
    double r_pow2 = R * R;
    for (int b = 0; b < N; ++b)
    {
        double sx = 0.0;
        double sy = 0.0;
        int count = 0;

#pragma omp parallel for
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
