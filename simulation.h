#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
std::vector<double> simulation(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> theta, const double R);
std::vector<double> simulation_openmp(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> theta, const double R);
#endif
