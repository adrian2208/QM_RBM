#pragma once

//Runtype
void SelectedRuns();

void GradientDescent(double LearningRate, int NrHiddenNodes, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions, int NrMetropolisSteps, double equilibrationFraction, double StepSize, bool Adaptive, int seed);
void ParameterSweep(double LR_min, double LR_max, int NH_min, int NH_max, int Nrsamples, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions, int NrMetropolisSteps, double equilibrationFraction, double StepSize, bool Adaptive, int seed);
void sigmaSweep();
void PositionSampling(double LearningRate, int NrHiddenNodes, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions, int NrMetropolisSteps, double equilibrationFraction, double StepSize, bool Adaptive, int seed);
void PositionSamplingFromParameters(int NrHiddenNodes, bool interacting, bool importanceSampling, int NrParticles, int NrDimensions, int NrMetropolisSteps, double equilibrationFraction, double StepSize, int seed, std::vector<double> WFParameters);
//math Tool
std::vector<double> linspace(double start, double stop, int NrVals);
