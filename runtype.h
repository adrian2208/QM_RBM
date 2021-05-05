#pragma once

//Runtype
void BruteForce_MC(int argc, char const* argv[]);
void BruteForce_MC_Numerical(int argc, char const* argv[]);
void VMC_w_ImportanceSampling(int argc, char const* argv[]);
void SimpleGaussianGD(int argc, char const* argv[]);
void BruteForce_MC_Correlated(int argc, char const* argv[]);
void VMC_w_ImportanceSampling_Correlated(int argc, char const* argv[]);
void CorrelatedGaussianGD(int argc, char const* argv[]);
void ParallelCorrelated(int argc, char const* argv[]);
void Plot2D_Positions(int argc, char const* argv[]);


//Outputhandling
void printEnergyToFile(std::string WF, std::string H, std::string fileOptString, std::vector<double> energies);
void printAcceptedToFile(std::string WF, std::string H, std::string fileOptString, double accepted);
void printTimeToFile(std::string WF, std::string H, std::string fileOptString, double elapsed);
void printMatrixToFile(std::string WF, std::string H, std::string fileOptString, int** matrix, int NrSamplingLengths);