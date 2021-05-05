#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <omp.h>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/simplegaussian_numerical.h"
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

int main(int argc, char const *argv[]) {

    //BruteForce_MC(argc, argv);
    //BruteForce_MC_Numerical(argc, argv);
    //BruteForce_MC_Correlated(argc, argv);
    //VMC_w_ImportanceSampling(argc, argv);
    //SimpleGaussianGD(argc, argv);
    //VMC_w_ImportanceSampling_Correlated(argc, argv);
    //CorrelatedGaussianGD(argc, argv);
    //ParallelCorrelated(argc, argv);
    Plot2D_Positions(argc, argv);
    return 0;
}

//////////  FUNCTION DEFINITIONS /////////

void BruteForce_MC(int argc, char const* argv[]) {
    std::cout << "Running Brute-Force Monte Carlo....\n";

    // Seed for the random number generator
    int seed = 2020;
    
    int numberOfDimensions = 3;
    int numberOfParticles = 3;
    int numberOfSteps = pow(2,15);
    double omega = 1.0;           // Oscillator frequency.
    double alpha = 0.5;           // Variational parameter.
    double stepLength = 0.1;       // Metropolis step length.
    double equilibration = 0.8;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    double timeStep = 0.01;
    std::string FileOptString= ""; // If empty, no file is printed
    

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
        Timer timer("Brute-Force MC",system);
        system->runMetropolisSteps(numberOfSteps);
    }
    int** matrix = system->getSampler()->getParticlePos_matrix();
    
    printMatrixToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, matrix, NrSamplingLengths);

    elapsed = system->getElapsedTime();
    printEnergyToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, system->getSampler()->getEnergyVector());
    printTimeToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, elapsed);
    printAcceptedToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, system->getSampler()->getNumberAccepted() / numberOfSteps);
}

void BruteForce_MC_Numerical(int argc, char const* argv[]) {
    std::cout << "Running Numerical Brute-Force Monte Carlo....\n";

    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles = 10;
    int numberOfSteps = pow(2, 18);
    double omega = 1.0;           // Oscillator frequency.
    double alpha = 0.5;           // Variational parameter.
    double stepLength = 0.1;       // Metropolis step length.
    double equilibration = 0.8;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    std::string FileOptString = ""; // If empty, no file is printed


    if (argc > 1) {
        try {
            numberOfDimensions = atoi(argv[1]);
            numberOfParticles = atoi(argv[2]);
            numberOfSteps = atoi(argv[3]);
            alpha = atof(argv[4]);        // Variational parameter.
            FileOptString = argv[5];     //Provides an additional identifier at the end of the output file name
        }
        catch (int c) {
            std::cout << "Error: " << c << ". Enough input arguments?";
        }
    }


    System* system = new System(seed);
    system->setHamiltonian(new HarmonicOscillator(system, omega));
    system->setWaveFunction(new SimpleGaussian_Numerical(system, alpha));
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setDriftCoefficient(stepLength);
    system->setFileOptString(FileOptString);
    double elapsed;
    {
        Timer timer("Brute-Force MC", system);
        system->runMetropolisSteps(numberOfSteps);
    }
    elapsed = system->getElapsedTime();
    printEnergyToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, system->getSampler()->getEnergyVector());
    printTimeToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, elapsed);
}

void BruteForce_MC_Correlated(int argc, char const* argv[]) {
    std::cout << "Running Brute-Force Monte Carlo....\n";

    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles = 10;
    int numberOfSteps = (int)1e6;
    double gamma = 2.82843;
    double hardCoreDiameter = 0.00433;
    double alpha = 0.5;           // Variational parameter.
    double beta = 2.82843;
    double stepLength = 0.1;       // Metropolis step length.
    double equilibration = 0.1;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    std::string FileOptString = ""; // If empty, no file is printed

    System* system = new System(seed);
    system->setHamiltonian(new EllipticOscillator(system, gamma, hardCoreDiameter, numberOfDimensions));
    system->setWaveFunction(new CorrelatedGaussian(system, alpha, beta, hardCoreDiameter));
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setDriftCoefficient(stepLength);
    system->setFileOptString(FileOptString);
    {
        Timer timer("Brute-Force MC",system);
        system->runMetropolisSteps(numberOfSteps);
    }
}

void VMC_w_ImportanceSampling(int argc, char const* argv[]) {
    std::cout << "Running VMC /w Importance Sampling....\n";

    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles = 10;
    int numberOfSteps = (int)1e6;
    double omega = 1.0;           // Oscillator frequency.
    double alpha = 0.5;           // Variational parameter.
    double stepLength = 0.1;       // Metropolis step length.
    double timeStep = 0.01;
    double equilibration = 0.1;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    std::string FileOptString = ""; // If empty, no file is printed

    if (argc > 1) {
        try {
            numberOfDimensions = atoi(argv[1]);
            numberOfParticles = atoi(argv[2]);
            numberOfSteps = atoi(argv[3]);
            alpha = atof(argv[4]);        // Variational parameter.
            timeStep = atof(argv[5]);
            FileOptString = argv[6];     //Provides an additional identifier at the end of the output file name
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
    //system->setStepLength(stepLength);
    system->setTimeStep(timeStep);
    system->setDriftCoefficient(stepLength);
    system->setFileOptString(FileOptString);
    double elapsed;
    {
        Timer timer("VMC_w_ImportanceSampling", system);
        system->runMALASteps(numberOfSteps);
    }
    elapsed = system->getElapsedTime();
    printEnergyToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, system->getSampler()->getEnergyVector());
    printTimeToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, elapsed);
    printAcceptedToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, system->getSampler()->getNumberAccepted()/ numberOfSteps);
}

void VMC_w_ImportanceSampling_Correlated(int argc, char const* argv[]) {
    std::cout << "Running VMC /w Importance Sampling (Correlated wave function)....\n";

    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles = 1;
    int numberOfSteps = (int)1e5;
    double gamma = 2.82843;
    double hardCoreDiameter = 0.00433;
    double alpha = 0.5;           // Variational parameter.
    double beta = 2.82843;
    double stepLength = 0.1;       // Metropolis step length.
    double timeStep = 0.01;
    double equilibration = 0.1;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    std::string FileOptString = ""; // If empty, no file is printed

    if (argc > 1) {
        try {
            numberOfDimensions = atoi(argv[1]);
            numberOfParticles = atoi(argv[2]);
            numberOfSteps = atoi(argv[3]);
            alpha = atof(argv[4]);        // Variational parameter.
            timeStep = atof(argv[5]);
            FileOptString = argv[6];     //Provides an additional identifier at the end of the output file name
        }
        catch (int c) {
            std::cout << "Error: " << c << ". Enough input arguments?";
        }
    }

    System* system = new System(seed);
    system->setHamiltonian(new EllipticOscillator(system, gamma, hardCoreDiameter, numberOfDimensions));
    system->setWaveFunction(new CorrelatedGaussian(system, alpha, beta, hardCoreDiameter));
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStep(timeStep);
    system->setDriftCoefficient(stepLength);
    system->setFileOptString(FileOptString);
    double elapsed;
    {
        Timer timer("VMC /w Importance Sampling (Correlated wave function)",system);
        system->runMALASteps(numberOfSteps);
        
    }
    elapsed = system->getElapsedTime();
    int accepted = system->getSampler()->getNumberAccepted();
    printEnergyToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, system->getSampler()->getEnergyVector());
    printTimeToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, elapsed);
    printAcceptedToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, accepted/ numberOfSteps);
}

void SimpleGaussianGD(int argc, char const* argv[]) {
    std::cout << "Running SimpleGaussian w/ Gradient Descent....\n";

    // Seed for the random number generator
    int seed = 2020;

    // Gradient descent parameters
    double alpha0 = 0.4;       // Initial Guess for the variational parameter
    double learningRate = 1e-2;
    int maxNrSteps = 40;
    double tolerance = 1e-6;

    std::vector<double> alpha;


    // System Parameters
    int numberOfDimensions = 3;
    int numberOfParticles = 10;
    int numberOfSteps = (int)1e6;
    double omega = 1.0;           // Oscillator frequency.
    double stepLength = 0.1;       // Metropolis step length.
    double equilibration = 0.1;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    std::string FileOptString = ""; // If empty, no file is printed

    if (argc > 1) {
        try {
            numberOfDimensions = atoi(argv[1]);
            numberOfParticles = atoi(argv[2]);
            numberOfSteps = atoi(argv[3]);
            omega = atof(argv[4]);        // Oscillator frequency.
            alpha0 = atof(argv[5]);        // Variational parameter.
            stepLength = atof(argv[6]);   // Metropolis step length.
            equilibration = atof(argv[7]);// Amount of the total steps used
            driftCoeff = atof(argv[8]);   //coeff. D used in MALA
            FileOptString = argv[9];     //Provides an additional identifier at the end of the output file name
        }
        catch (int c) {
            std::cout << "Error: " << c << ". Enough input arguments?";
        }
    }
    alpha.push_back(alpha0);

    int i = 0;

    do {
        double alpha_current = alpha[i];

        System* system = new System(seed);
        system->setHamiltonian(new HarmonicOscillator(system, omega));
        system->setWaveFunction(new SimpleGaussian(system, alpha_current));
        system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setSampler(new sampler_gd(system));
        system->setEquilibrationFraction(equilibration);
        system->setStepLength(stepLength);
        system->setDriftCoefficient(stepLength);
        system->setFileOptString(FileOptString);
        system->runMetropolisSteps(numberOfSteps);

        double localEnergyGradient = system->getWaveFunction()->localEnergyDerivative();
        alpha_current -= learningRate * localEnergyGradient;
        std::cout << "local energy gradient = " << localEnergyGradient << "\n";
        alpha.push_back(alpha_current);
        
        i++;

    } while (i < maxNrSteps && tolerance < std::abs(alpha[i] - alpha[i - 1]));

    for (double item : alpha) {
        std::cout << item << "\n";
    }
}

void CorrelatedGaussianGD(int argc, char const* argv[]) {
    std::cout << "Running CorrelatedGaussian w/ Gradient Descent....\n";

    // Seed for the random number generator
    int seed = 2020;

    // Gradient descent parameters
    double alpha0 = 0.40946;       // Initial Guess for the variational parameter
    double learningRate = 5e-5;
    int maxNrSteps = 50;
    double tolerance = 1e-6;

    std::vector<double> alpha;


    // System Parameters
    int numberOfDimensions = 3;
    int numberOfParticles = 30;
    int numberOfSteps = (int)1e5;
    double gamma = 2.82843;
    double hardCoreDiameter = 0.00433;
    double beta = 2.82843;
    double stepLength = 0.1;       // Metropolis step length.
    double equilibration = 0.1;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    std::string FileOptString = ""; // If empty, no file is printed

    if (argc > 1) {
        try {
            numberOfDimensions = atoi(argv[1]);
            numberOfParticles = atoi(argv[2]);
            numberOfSteps = atoi(argv[3]);
            alpha0 = atof(argv[5]);        // Variational parameter.
            stepLength = atof(argv[6]);   // Metropolis step length.
            equilibration = atof(argv[7]);// Amount of the total steps used
            driftCoeff = atof(argv[8]);   //coeff. D used in MALA
            FileOptString = argv[9];     //Provides an additional identifier at the end of the output file name
        }
        catch (int c) {
            std::cout << "Error: " << c << ". Enough input arguments?";
        }
    }

    alpha.push_back(alpha0);

    int i = 0;

    do {
        double alpha_current = alpha[i];

        System* system = new System(seed);
        system->setHamiltonian(new EllipticOscillator(system, gamma, hardCoreDiameter, numberOfDimensions));
        system->setWaveFunction(new CorrelatedGaussian(system, alpha_current, beta, hardCoreDiameter));
        system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setSampler(new sampler_gd(system));
        system->setEquilibrationFraction(equilibration);
        system->setStepLength(stepLength);
        system->setDriftCoefficient(stepLength);
        system->setFileOptString(FileOptString);
        system->runMetropolisSteps(numberOfSteps);

        double localEnergyGradient = system->getWaveFunction()->localEnergyDerivative();
        alpha_current -= learningRate * localEnergyGradient;
        std::cout << "local energy gradient = " << localEnergyGradient << "\n";
        alpha.push_back(alpha_current);

        i++;

    } while (i < maxNrSteps && tolerance < std::abs(alpha[i] - alpha[i - 1]));

    for (double item : alpha) {
        std::cout << item << "\n";
    }
}

void ParallelCorrelated(int argc, char const* argv[]) {
    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles = 30;
    int numberOfSteps = (int)1e6;
    double gamma = 2.82843;
    double hardCoreDiameter = 0.00433;
    double beta = 2.82843;
    double stepLength = 0.1;       // Metropolis step length.
    double equilibration = 0.1;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    std::string FileOptString = ""; // If empty, no file is printed
    
    double alpha_min = 0.40;
    double alpha_max = 0.415;


    int NrProcessors_Available = omp_get_num_procs();
    std::cout << NrProcessors_Available << "\n";

    int    NrAlphas = NrProcessors_Available;

    std::vector<double> alphas(NrAlphas);           // Variational parameter.
    double step = (alpha_max - alpha_min) / (NrAlphas - 1);
    for (int i = 0; i < NrAlphas; i++) {
        alphas[i] = alpha_min + i * step;
    }

    #pragma omp parallel for 
    for (int i = 0; i < NrAlphas; i++) {
            System* system = new System(seed);
            system->setHamiltonian(new EllipticOscillator(system, gamma, hardCoreDiameter, numberOfDimensions));
            system->setWaveFunction(new CorrelatedGaussian(system, alphas[i], beta, hardCoreDiameter));
            system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
            system->setSampler(new Sampler(system));
            system->setEquilibrationFraction(equilibration);
            system->setStepLength(stepLength);
            system->setDriftCoefficient(stepLength);
            system->setFileOptString(FileOptString);
            system->runMetropolisSteps(numberOfSteps);
        }

}

void Plot2D_Positions(int argc, char const* argv[]) {
    std::cout << "Running....\n";

    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;

    double gamma = 2.82843;
    //double hardCoreDiameter = 0.00433;
    double beta = 2.82843;
    double stepLength = 0.1;       // Metropolis step length.
    double equilibration = 0.8;    // Amount of the total steps used
    double driftCoeff = 0.5;       // coeff. D used in MALA
    double timeStep = 0.01;
    std::string FileOptString = ""; // If empty, no file is printed

    std::cout << argv[1] << "\n";
    int numberOfParticles = atoi(argv[1]);
    int numberOfSteps = atoi(argv[2]);
    double alpha = atof(argv[3]);        // Variational parameter.
    double hardCoreDiameter = atof(argv[4]);
    int NrSamplingLengths = atoi(argv[5]);
    int SamplingRadius = atoi(argv[6]);
    FileOptString = argv[7];



    System* system = new System(seed);
    system->setHamiltonian(new EllipticOscillator(system, gamma, hardCoreDiameter, numberOfDimensions));
    system->setWaveFunction(new CorrelatedGaussian(system, alpha, beta, hardCoreDiameter));
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setSampler(new Sampler(system));
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setDriftCoefficient(stepLength);
    system->setFileOptString(FileOptString);

    system->getSampler()->SetupPositionSampling(NrSamplingLengths, SamplingRadius);

    system->runMetropolisSteps(numberOfSteps);

    
    int** matrix = system->getSampler()->getParticlePos_matrix();
    printMatrixToFile(system->getWaveFunction()->getName(), system->getHamiltonian()->getName(), FileOptString, matrix, NrSamplingLengths);

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
            OutFile << "," << matrix[i][j] ;
        }
        OutFile << "\n";
    }

    OutFile.close();
}