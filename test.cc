#include <iostream>
#include <assert.h>

#include "simulation.h"

const double R = 1.0;
const double threshold = 0.0001;
const std::vector<double> CORRECT_RESULT = {1.0, 1.0, 1.0, 1.0};

void init_vectors(std::vector<double> &x, std::vector<double> &y, std::vector<double> &theta)
{
    x = {1.0, 0.0, 1.0, 1.0};
    y = {1.0, 1.0, 0.0, 0.3};
    theta = {1.0, 1.0, 1.0, 1.0};
}

void test_simulation()
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> theta;
    init_vectors(x, y, theta);

    std::cout << "Testing simulation function" << std::endl;
    std::vector<double> res = simulation(x, y, theta, R);

    assert(res.size() == CORRECT_RESULT.size());
    for (int i = 0; i < res.size(); i++)
    {
        // print type of res[i] and CORRECT_RESULT[i]
        double diff = res[i] - CORRECT_RESULT[i];
        if (diff > threshold || diff < -threshold)
        {
            std::cout << "res[i]: " << res[i] << " CORRECT_RESULT[i]: " << CORRECT_RESULT[i] << std::endl;
            std::cout << "diff: " << diff << std::endl;
            assert(res[i] == CORRECT_RESULT[i]);
        }
    }
}

void test_simulation_openmp()
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> theta;
    init_vectors(x, y, theta);

    std::cout << "Testing simulation_openmp function" << std::endl;
    std::vector<double> res = simulation_openmp(x, y, theta, R);

    assert(res.size() == CORRECT_RESULT.size());

    assert(res == CORRECT_RESULT);
}

void test_simulation_mpi()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get each process rank

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> theta;
    init_vectors(x, y, theta);
    if (rank == 0)
        std::cout << "Testing simulation_mpi function" << std::endl;

    std::vector<double> res = simulation_mpi(x, y, theta, R);

    if (rank == 0)
    { // Only assert the result on the root process
        assert(res.size() == CORRECT_RESULT.size());
        assert(res == CORRECT_RESULT);
    }
}

int main(int argc, char *argv[])
{
    int provided, rank;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get each process rank

    if (rank == 0)
    {
        std::cout << "Running tests..." << std::endl;
    }

    if (rank == 0)
    {
        test_simulation();
//        test_simulation_openmp();
    }

  //  test_simulation_mpi();

    if (rank == 0)
        std::cout << "All tests passed!" << std::endl;

    MPI_Finalize();
}
