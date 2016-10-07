#ifndef SVLIB_SIMULATE_H
#define SVLIB_SIMULATE_H

#include <string>
#include <vector>

void parseSimulationOptions(int argc, char** argv);
void runSimulation(int argc, char** argv);
std::string genBreaks();
std::vector<double> parseErrorRates(const std::string& s);
std::string errorRateString(const std::vector<double>& v, const std::string& name);

#endif
