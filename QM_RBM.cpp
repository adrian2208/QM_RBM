#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <utility>
#include <omp.h>
#include <numeric>
#include "system.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/qnet.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/InteractingOscillator.h"
#include "Math/random.h"
#include "Benchmark.h"
#include "csvHandler.h"
#include "QM_RBM.h"
#include "Optimizer/Optimizer.h"
#include "Optimizer/Sgd.h"
#include "Optimizer/Adam.h"

using namespace std;

//////////  MAIN BODY  /////////

int main(int argc, char const* argv[]) {
    
    SelectedRuns();

    return 0;
}

//////////  FUNCTION DEFINITIONS /////////

void SelectedRuns() {
    GradientDescent(0.1, 3, true, 40, true, false, 2, 2, pow(2, 21), 1.0, 0.01, true, 1111);
    PositionSampling(0.1, 3, true, 10, true, true, 2, 2, pow(2, 21), 1.0, 0.01, true, 1111);
}

void PositionSampling(double LearningRate, int NrHiddenNodes, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions, int NrMetropolisSteps, double equilibrationFraction, double StepSize, bool Adaptive, int seed) {
    double sigma = 1.0;
    double sigma0 = 0.001;
    double omega = 1.0;
    int posSampLen = 300;

    System* system = new System(seed);
    if (interacting) {
        system->setHamiltonian(new InteractingOscillator(system, omega, NrParticles));
    }
    else {
        system->setHamiltonian(new HarmonicOscillator(system, omega));
    }

    system->setWaveFunction(new qnet(system, NrDimensions, NrParticles, NrHiddenNodes, sigma, sigma0));

    if (Adaptive) {
        system->getWaveFunction()->setOptimizer(new Adam(system->getWaveFunction()->getNrParameters(), LearningRate));
    }
    else {
        system->getWaveFunction()->setOptimizer(new Sgd(system->getWaveFunction()->getNrParameters(), LearningRate));
    }

    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibrationFraction);
    int NrParameters = system->getWaveFunction()->getNrParameters();
    //Outputfile Management
    std::vector<std::pair<std::string, std::vector<double>>> datastruct = std::vector<std::pair<std::string, std::vector<double>>>();
    std::string descriptor = "GD_ls_v_E_LR_" + to_string(LearningRate) + "_NH_" + to_string(NrHiddenNodes);
    std::string fileName;
    std::string csvComment = "#Energy over gradient descent step\n# HyperParameters: LR = " + to_string(LearningRate) + ", NH = " + to_string(NrHiddenNodes) + ", sigma = " + to_string(sigma);

    system->runOptimizationSteps(NrMetropolisSteps, NrLearningSteps, StepSize, LearningRate, datastruct, importanceSampling, OnlyLastEnergies);
    
    system->setSampler(new Sampler(system));
    system->getSampler()->SetupPositionSampling(posSampLen, 3);
    

    if (importanceSampling) {
        system->runMALASteps(NrMetropolisSteps, StepSize);
    }
    else {
        system->runMetropolisSteps(NrMetropolisSteps, StepSize);
    }
    std::vector<std::pair<std::string, std::vector<int>>> ParticleRadii = std::vector<std::pair<std::string, std::vector<int>>>();

    //Create the Output file
    fileName = system->getSampler()->GenerateFileName(descriptor);
    csvHandler outputFile(fileName, csvComment);
    outputFile.WriteMatrixToFile(system->getSampler()->getParticlePos_matrix(), posSampLen);
    ParticleRadii.push_back({ "radii",system->getSampler()->getParticleRadiusVector() });
    csvHandler radiiFile("Particle_Radii", "");
    radiiFile.WriteToCSV(ParticleRadii);
}



void PositionSamplingFromParameters(int NrHiddenNodes, bool interacting,bool importanceSampling, int NrParticles, int NrDimensions, int NrMetropolisSteps, double equilibrationFraction, double StepSize,int seed, std::vector<double> WFParameters) {
    double sigma = 1.0;
    double sigma0 = 0.001;
    double omega = 1.0;
    int posSampLen = 300;

    System* system = new System(seed);
    if (interacting) {
        system->setHamiltonian(new InteractingOscillator(system, omega, NrParticles));
    }
    else {
        system->setHamiltonian(new HarmonicOscillator(system, omega));
    }

    system->setWaveFunction(new qnet(system, NrDimensions, NrParticles, NrHiddenNodes, sigma, sigma0));

    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibrationFraction);
    int NrParameters = system->getWaveFunction()->getNrParameters();

    system->getSampler()->SetupPositionSampling(posSampLen, 3);

    system->getWaveFunction()->setParameters(WFParameters);

    if (importanceSampling) {
        system->runMALASteps(NrMetropolisSteps, StepSize);
    }
    else {
        system->runMetropolisSteps(NrMetropolisSteps, StepSize);
    }
    std::vector<std::pair<std::string, std::vector<int>>> ParticleRadii = std::vector<std::pair<std::string, std::vector<int>>>();

    //Create the Output file
    std::string descriptor = "Position_Sampling";
    std::string fileName;
    std::string csvComment = "#Position Sampling";

    fileName = system->getSampler()->GenerateFileName(descriptor);
    csvHandler outputFile(fileName, csvComment);
    outputFile.WriteMatrixToFile(system->getSampler()->getParticlePos_matrix(), posSampLen);
    ParticleRadii.push_back({ "radii",system->getSampler()->getParticleRadiusVector() });
    csvHandler radiiFile("Particle_Radii", "");
    radiiFile.WriteToCSV(ParticleRadii);
}



void ParameterSweep(double LR_min, double LR_max, int NH_min, int NH_max, int Nrsamples, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions, int NrMetropolisSteps, double equilibrationFraction, double StepSize, bool Adaptive, int seed) {
    int NrThreads = 6;
    std::vector<double> LRs = linspace(LR_min, LR_max, Nrsamples);
    int NrNHS = NH_max - NH_min + 1;
    std::vector<int> NHs(NrNHS);
    std::iota(std::begin(NHs), std::end(NHs), NH_min);
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < Nrsamples; i++) {
        for (int j = 0; j < NrNHS; j++) {
            int id = omp_get_thread_num();
            std::cout << "thread id: " << id << "\n";
            GradientDescent(LRs[i], NHs[j], interacting, NrLearningSteps,importanceSampling,OnlyLastEnergies,NrParticles,NrDimensions,NrMetropolisSteps,equilibrationFraction,StepSize,Adaptive,seed);
        }
    }
}

void GradientDescent(double LearningRate, int NrHiddenNodes, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions, int NrMetropolisSteps, double equilibrationFraction,double StepSize, bool Adaptive, int seed) {
    double sigma = 1.0;
    double sigma0 = 0.001;
    double omega = 1.0;
   

    System* system = new System(seed);
    if (interacting) {
        system->setHamiltonian(new InteractingOscillator(system, omega, NrParticles));
    }
    else {
        system->setHamiltonian(new HarmonicOscillator(system, omega));
    }

    system->setWaveFunction(new qnet(system, NrDimensions, NrParticles, NrHiddenNodes, sigma, sigma0));

    if (Adaptive) {
        system->getWaveFunction()->setOptimizer(new Adam(system->getWaveFunction()->getNrParameters(), LearningRate));
    }
    else {
        system->getWaveFunction()->setOptimizer(new Sgd(system->getWaveFunction()->getNrParameters(), LearningRate));
    }

    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibrationFraction);
    int NrParameters = system->getWaveFunction()->getNrParameters();
    //Outputfile Management
    std::vector<std::pair<std::string, std::vector<double>>> datastruct = std::vector<std::pair<std::string, std::vector<double>>>();
    std::string descriptor = "GD_ls_v_E_LR_" + to_string(LearningRate) + "_NH_" + to_string(NrHiddenNodes) + "_Adaptive_" + to_string(Adaptive);;
    std::string fileName;
    std::string csvComment = "#Energy over gradient descent step\n# HyperParameters: LR = " + to_string(LearningRate) + ", NH = " + to_string(NrHiddenNodes) + ", sigma = " + to_string(sigma);
    
    system->runOptimizationSteps(NrMetropolisSteps, NrLearningSteps, StepSize, LearningRate, datastruct, importanceSampling,OnlyLastEnergies);

    //Create the Output file
    fileName = system->getSampler()->GenerateFileName(descriptor);
    csvHandler outputFile(fileName, csvComment);
    outputFile.WriteToCSV(datastruct);
}

void sigmaSweep() {

    bool interaction = true;
    bool importanceSampling = true;
    int NrDimensions = 2;
    int NrParticles = 2;
    int NrHiddenNodes = 2;

    double equilibrationFraction = 0.5;
    int NrMetropolisSteps = pow(2, 18);

    //Hyper Parameters
    double sigma0 = 0.001;
    double omega = 1.0;
    double stepLength = 0.2;
    double timeStep = 0.3;

    double startSigma = 0.8;
    double stopSigma = 2.8;
    int NrSigmas = 50;
    std::vector<double> sigma = linspace(startSigma,stopSigma,NrSigmas);

    //Outputfile Management
    std::vector<std::pair<std::string, std::vector<double>>> datastruct = std::vector<std::pair<std::string, std::vector<double>>>();
    std::string descriptor = "test";
    std::string fileName;
    std::string csvComment = "#columns of energies for blocking sampled at different sigma\n# HyperParameters: SL = " + to_string(stepLength) + ", TS = " + to_string(timeStep)+", NH = " + to_string(NrHiddenNodes);
    datastruct.push_back({ "Sigma vals",sigma });

    for (int i = 0; i < NrSigmas; i++) {

        //setup
        System* system = new System();
        if (interaction){
            system->setHamiltonian(new InteractingOscillator(system, omega, NrParticles));
        }
        else{
            system->setHamiltonian(new HarmonicOscillator(system, omega));
        }
        system->setWaveFunction(new qnet(system, NrDimensions, NrParticles, NrHiddenNodes, sigma[i], sigma0));
        system->setSampler(new Sampler(system));
        system->setEquilibrationFraction(equilibrationFraction);
        
        //run
        if (importanceSampling){
            system->runMALASteps(NrMetropolisSteps, timeStep);
        }
        else {
            system->runMetropolisSteps(NrMetropolisSteps, stepLength);
        }
        //Data Handling
        datastruct.push_back({ "sigma:" + to_string(sigma[i]),system->getSampler()->getEnergyVector() });
        fileName = system->getSampler()->GenerateFileName(descriptor);
    }
    //Create the Output file
    csvHandler outputFile(fileName,csvComment);
    outputFile.WriteToCSV(datastruct);
}



///////////////////////////////////////////
std::vector<double> linspace(double start, double stop, int NrVals) {
    std::vector<double> output(NrVals);
    double dist = (stop - start) / (NrVals - 1);
    output[0] = start;
    for (int i = 0; i < NrVals-1; i++) {
        output[i + 1] = output[i] + dist;
    }
    return output;
}

