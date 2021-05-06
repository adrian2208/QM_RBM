#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <omp.h>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/correlatedGaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/ellipticOscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "Benchmark.h"
#include "runtype.h"
#include "sampler_gd.h"

using namespace std;

//////////  MAIN BODY  /////////

int main(int argc, char const* argv[]) {

    BruteForce_MC(argc, argv);

    return 0;
}

//////////  FUNCTION DEFINITIONS /////////

void BruteForce_MC(int argc, char const* argv[]) {
    std::cout << "Running Brute-Force Monte Carlo....\n";

    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles = 3;
    int numberOfSteps = pow(2, 15);
    double omega = 1.0;           // Oscillator frequency.
    double alpha = 0.5;           // Variational parameter.
    double stepLength = 0.1;       // Metropolis step length.
    double equilibration = 0.8;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    double timeStep = 0.01;
    std::string FileOptString = ""; // If empty, no file is printed


    if (argc > 1) {
        try {
            numberOfDimensions = atoi(argv[1]);
            numberOfParticles = atoi(argv[2]);
            numberOfSteps = atoi(argv[3]);
            alpha = atof(argv[4]);        // Variational parameter.
            timeStep = atof(argv[5]);     //Provides an additional identifier at the end of the output file name
            FileOptString = argv[6];
        }
        catch (int c) {
            std::cout << "Error: " << c << ". Enough input arguments?";
        }
    }


    System* system = new System(seed);
    system->setHamiltonian(new HarmonicOscillator(system, omega));
    system->setWaveFunction(new SimpleGaussian(system, alpha));
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setDriftCoefficient(stepLength);
    system->setFileOptString(FileOptString);

    int NrSamplingLengths = 50;
    system->getSampler()->SetupPositionSampling(NrSamplingLengths, 3);

    double elapsed;
    {
        Timer timer("Brute-Force MC", system);
        system->runMetropolisSteps(numberOfSteps);
    }
    int** matrix = system->getSampler()->getParticlePos_matrix();

    printMatrixToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, matrix, NrSamplingLengths);

    elapsed = system->getElapsedTime();
    printEnergyToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, system->getSampler()->getEnergyVector());
    printTimeToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, elapsed);
    printAcceptedToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, system->getSampler()->getNumberAccepted() / numberOfSteps);
}


///////////////////////////////////////////


void printEnergyToFile(std::string WF, std::string H, std::string fileOptString, std::vector<double> energies) {
    //std::string WF = system->getWaveFunction()->getName();
    //std::string H = system->getHamiltonian()->getName();


    //IMPORTANT: when this program is executed through python, the root directory is assumed to be 
    //the directory of the python program, not the executable, meaning the path below 
    // will output to variational-monte-carlo-fys4411/Output when called from python 
    // but to a subdirectory in the x64-Release directory when run from the executable

    std::string path = ".\\Output\\" + WF + "_" + H + "_" + fileOptString + "_Energies" + ".csv";
    std::ofstream OutFile;
    OutFile.open(path, std::fstream::app);
    for (double energy : energies) {
        OutFile << energy << "\n";
    }
    OutFile.close();
}

void printAcceptedToFile(std::string WF, std::string H, std::string fileOptString, double accepted) {
    //IMPORTANT: when this program is executed through python, the root directory is assumed to be 
    //the directory of the python program, not the executable, meaning the path below 
    // will output to variational-monte-carlo-fys4411/Output when called from python 
    // but to a subdirectory in the x64-Release directory when run from the executable

    std::string path = ".\\Output\\" + WF + "_" + H + "_" + fileOptString + "_accepted" + ".csv";
    std::ofstream OutFile;
    OutFile.open(path, std::fstream::app);
    OutFile << accepted << "\n";
    OutFile.close();
}

void printTimeToFile(std::string WF, std::string H, std::string fileOptString, double elapsed) {
    //std::string WF = system->getWaveFunction()->getName();
    //std::string H = system->getHamiltonian()->getName();


    //IMPORTANT: when this program is executed through python, the root directory is assumed to be 
    //the directory of the python program, not the executable, meaning the path below 
    // will output to variational-monte-carlo-fys4411/Output when called from python 
    // but to a subdirectory in the x64-Release directory when run from the executable

    std::string path = ".\\Output\\" + WF + "_" + H + "_" + fileOptString + "_ElapsedTime" + ".csv";
    std::ofstream OutFile;
    OutFile.open(path, std::fstream::app);
    OutFile << elapsed << "\n";

    OutFile.close();
}



void printMatrixToFile(std::string WF, std::string H, std::string fileOptString, int** matrix, int NrSamplingLengths) {
    //std::string WF = system->getWaveFunction()->getName();
    //std::string H = system->getHamiltonian()->getName();


    //IMPORTANT: when this program is executed through python, the root directory is assumed to be 
    //the directory of the python program, not the executable, meaning the path below 
    // will output to variational-monte-carlo-fys4411/Output when called from python 
    // but to a subdirectory in the x64-Release directory when run from the executable

    std::string path = ".\\Output\\" + WF + "_" + H + "_" + fileOptString + "_Matrix" + ".csv";
    std::ofstream OutFile;
    OutFile.open(path, std::fstream::app);
    for (int i = 0; i < NrSamplingLengths; i++) {
        OutFile << matrix[i][0];
        for (int j = 1; j < NrSamplingLengths; j++) {
            OutFile << "," << matrix[i][j];
        }
        OutFile << "\n";
    }

    OutFile.close();
}