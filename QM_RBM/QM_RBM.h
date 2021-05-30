#pragma once

//Runtype
void GradientDescent(double LearningRate, int NrHiddenNodes, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions);
void GradientDescent(double LearningRate, int NrHiddenNodes, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions,int seed);
void ParameterSweep(double LR_min, double LR_max, int NH_min, int NH_max, int Nrsamples, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions);
void sigmaSweep();
void PositionSampling(int NrParticles, int NrHiddenNodes, bool interacting, bool importanceSampling, int seed, bool optimize);

//math Tool
std::vector<double> linspace(double start, double stop, int NrVals);
