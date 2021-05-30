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

using namespace std;

//////////  MAIN BODY  /////////

int main(int argc, char const* argv[]) {
    
    //GradientDescent(0.1,2,true,80,true);
    //sigmaSearch();
    ParameterSweep(0.01, 0.01, 6, 8,1, true,200,true,true,2,2);
    //GradientDescent(0.05, 2, true, 200, true, false,2,2,2020);
    //GradientDescent(0.05, 2, true, 100, false, false, 2, 2, 2020);

    return 0;
}

//////////  FUNCTION DEFINITIONS /////////
void PositionSampling(int NrParticles,int NrHiddenNodes, bool interacting, bool importanceSampling,  int seed,bool optimize) {
    double sigma = 1.0;
    double sigma0 = 0.001;
    double omega = 1.0;
    double equilibrationFraction = 0.8;
    int NrMetropolisSteps = pow(2, 19);
    int NrDimensions = 2;

    //Hyper Parameters
    double stepLength = 0.1;
    double timeStep = 0.45;

    //Position sampling lengths
    int posSampLen = 300;

    System* system = new System();
    if (interacting) {
        system->setHamiltonian(new InteractingOscillator(system, omega, NrParticles));
    }
    else {
        system->setHamiltonian(new HarmonicOscillator(system, omega));
    }

    system->setWaveFunction(new qnet(system, NrDimensions, NrParticles, NrHiddenNodes, sigma, sigma0));
    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibrationFraction);

    if (optimize){
        std::vector<std::pair<std::string, std::vector<double>>> datastruct = std::vector<std::pair<std::string, std::vector<double>>>();
        system->runOptimizationSteps(NrMetropolisSteps, 200, timeStep, 0.01, datastruct, importanceSampling, true);
        system->setSampler(new Sampler(system));
    }

    system->getSampler()->SetupPositionSampling(posSampLen, 3);
    
    //Outputfile Management
    std::string descriptor = "Position_sampling_P_" + to_string(NrParticles) + "_NH_" + to_string(NrHiddenNodes) + "_I_" + to_string(interacting);
    std::string fileName;
    std::string csvComment = "#Position Sampling\n# HyperParameters:P" + to_string(NrParticles) + "_NH_" + to_string(NrHiddenNodes) + "_I_" + to_string(interacting);
    
    if (importanceSampling) {
        system->runMALASteps(NrMetropolisSteps, timeStep);
    }
    else {
        system->runMetropolisSteps(NrMetropolisSteps, stepLength);
    }

    //Create the Output file
    fileName = system->getSampler()->GenerateFileName(descriptor);
    csvHandler outputFile(fileName, csvComment);
    outputFile.WriteMatrixToFile(system->getSampler()->getParticlePos_matrix(), posSampLen);
}

void ParameterSweep(double LR_min, double LR_max, int NH_min, int NH_max, int Nrsamples, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions) {
    int NrThreads = 6;
    std::vector<double> LRs = linspace(LR_min, LR_max, Nrsamples);
    int NrNHS = NH_max - NH_min + 1;
    std::vector<int> NHs(NrNHS);
    std::iota(std::begin(NHs), std::end(NHs), NH_min);
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < Nrsamples; i++) {
        for (int j = 0; j < NrNHS; j++) {
            int id = omp_get_thread_num();
            std::cout << "thread id: " << id << "\n";// "   Learning Rate: " << LRs[i] << "   #Hidden: " << NHs[j] << "\n";
            GradientDescent(LRs[i], NHs[j], interacting, NrLearningSteps,importanceSampling,OnlyLastEnergies,NrParticles,NrDimensions);
        }
    }
}

void GradientDescent(double LearningRate, int NrHiddenNodes, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions) {
    double sigma = 1.0;
    double sigma0 = 0.001;
    double omega = 1.0;
    double equilibrationFraction = 0.8;
    int NrMetropolisSteps = pow(2, 19);
   
    //Hyper Parameters
    double stepLength = 0.1;
    double timeStep = 0.45;

    System* system = new System();
    if (interacting) {
        system->setHamiltonian(new InteractingOscillator(system, omega, NrParticles));
    }
    else {
        system->setHamiltonian(new HarmonicOscillator(system, omega));
    }

    system->setWaveFunction(new qnet(system, NrDimensions, NrParticles, NrHiddenNodes, sigma, sigma0));
    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibrationFraction);

    //Outputfile Management
    std::vector<std::pair<std::string, std::vector<double>>> datastruct = std::vector<std::pair<std::string, std::vector<double>>>();
    std::string descriptor = "GD_ls_v_E_LR_" + to_string(LearningRate) + "_NH_" + to_string(NrHiddenNodes);
    std::string fileName;
    std::string csvComment = "#Energy over gradient descent step\n# HyperParameters: LR = " + to_string(LearningRate) + ", NH = " + to_string(NrHiddenNodes) + ", sigma = " + to_string(sigma);
    system->runOptimizationSteps(NrMetropolisSteps, NrLearningSteps, timeStep, LearningRate, datastruct, importanceSampling,OnlyLastEnergies);
    //Create the Output file
    fileName = system->getSampler()->GenerateFileName(descriptor);
    csvHandler outputFile(fileName, csvComment);
    outputFile.WriteToCSV(datastruct);
}
void GradientDescent(double LearningRate, int NrHiddenNodes, bool interacting, int NrLearningSteps, bool importanceSampling, bool OnlyLastEnergies, int NrParticles, int NrDimensions, int seed) {
    double sigma = 1.0;
    double sigma0 = 0.001;
    double omega = 1.0;
    double equilibrationFraction = 0.8;
    int NrMetropolisSteps = pow(2, 19);

    //Hyper Parameters
    double stepLength = 0.1;
    double timeStep = 0.45;

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

    //Outputfile Management
    std::vector<std::pair<std::string, std::vector<double>>> datastruct = std::vector<std::pair<std::string, std::vector<double>>>();
    std::string descriptor = "GD_ls_v_E_LR_" + to_string(LearningRate) + "_NH_" + to_string(NrHiddenNodes);
    std::string fileName;
    std::string csvComment = "#Energy over gradient descent step\n# HyperParameters: LR = " + to_string(LearningRate) + ", NH = " + to_string(NrHiddenNodes) + ", sigma = " + to_string(sigma);
    system->runOptimizationSteps(NrMetropolisSteps, NrLearningSteps, timeStep, LearningRate, datastruct, importanceSampling, OnlyLastEnergies);
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

