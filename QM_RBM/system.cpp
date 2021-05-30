#include "system.h"
#include <cassert>
#include <iostream>
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "Math/random.h"


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

bool System::metropolisStep() {
    int NrParticles = m_waveFunction->getNrParticles();
    int particleIndex = m_random->nextInt(m_waveFunction->getNrParticles() - 1);
    std::vector<double> OldPos = m_waveFunction->getParticlePosition(particleIndex);
    double wfOld = m_waveFunction->evaluate(particleIndex);

    for (int i = 0; i < m_waveFunction->getNrDimensions(); i++) {
        m_waveFunction->adjustPosition(particleIndex*NrParticles+i,m_stepLength * (m_random->nextDouble()-0.5));
    }
    
    double wfNew = m_waveFunction->evaluate(particleIndex);

    if (wfNew * wfNew / (wfOld * wfOld) < m_random->nextDouble()) {
        for (int i = 0; i < m_waveFunction->getNrDimensions(); i++) {
            m_waveFunction->setPosition(particleIndex * NrParticles + i, OldPos[i]);
        }
        return false;
    }
    return true;
}

bool System::MALAStep(double timeStep) {
    double timeStepSqrt = sqrt(timeStep);
    int NrParticles = m_waveFunction->getNrParticles();
    int particleIndex = m_random->nextInt(m_waveFunction->getNrParticles() - 1);
    std::vector<double> OldPos = m_waveFunction->getParticlePosition(particleIndex);
    double wfOld = m_waveFunction->evaluate(particleIndex);
    std::vector<double> driftTerm = m_waveFunction->driftTerm(particleIndex);

    for (int i = 0; i < m_waveFunction->getNrDimensions(); i++) {
        m_waveFunction->adjustPosition(particleIndex * NrParticles + i,0.5 * driftTerm[i] * timeStep
            + m_random->nextGaussian(0, 1) * timeStepSqrt);
    }

    std::vector<double> driftTermNew = m_waveFunction->driftTerm(particleIndex);
    double Greens_exponent = 0;
    for (int i = 0; i < m_waveFunction->getNrDimensions(); i++) {
        Greens_exponent += 0.5 * (driftTerm[i] + driftTermNew[i]) *
            (0.5 * timeStep * 0.5 * (driftTerm[i]
                - driftTermNew[i]) - m_waveFunction->getParticlePosition(particleIndex)[i] + OldPos[i]);
    }

    double wfNew = m_waveFunction->evaluate(particleIndex);

    if (wfNew * wfNew* exp(Greens_exponent) / (wfOld * wfOld+1.0e-10) < m_random->nextDouble()) {
        for (int i = 0; i < m_waveFunction->getNrDimensions(); i++) {
            //std::cout << "particleIndex: " << particleIndex << "        OldPos[i]: " << OldPos[i] << "\n";
            m_waveFunction->setPosition(particleIndex* NrParticles + i, OldPos[i]);
        }
        //std::cout << "\n";
        return false;
    }
    return true;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps, double stepLength) {
    m_stepLength = stepLength;
    m_samplingType = "Metropolis";
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    
    bool acceptedStep;
    int NrEqillibration_steps = (int)(numberOfMetropolisSteps * m_equilibrationFraction);

    for (int i = 0; i < NrEqillibration_steps; i++) {
        acceptedStep = metropolisStep();
    }

    for (int i = 0; i < m_numberOfMetropolisSteps; i++) {
        acceptedStep = metropolisStep();
        m_sampler->sample(acceptedStep);
    }
    
    m_sampler->computeAverages();
    if (!m_fileOptString.empty()) {
        m_sampler->printOutputToFile();
    }
    else {
        m_sampler->printOutputToTerminal();
    }
}

void System::runMALASteps(int numberOfMetropolisSteps, double timeStep) {
    m_samplingType = "Importance";
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    bool acceptedStep;
    int NrEqillibration_steps = (int)(numberOfMetropolisSteps * m_equilibrationFraction);

    for (int i = 0; i < NrEqillibration_steps; i++) {
        acceptedStep = MALAStep(timeStep);
    }
    for (int i = 0; i < m_numberOfMetropolisSteps; i++) {
        acceptedStep = MALAStep(timeStep);
        m_sampler->sample(acceptedStep);
    }

    m_sampler->computeAverages();
    
    if (!m_fileOptString.empty()) {
        m_sampler->printOutputToFile();
    }
    else {
        m_sampler->printOutputToTerminal();
    }
}

void System::runOptimizationSteps(int numberOfMetropolisSteps, int numberOfOptimizationSteps, double timestep,double learningRate, std::vector<std::pair<std::string, std::vector<double>>> &datastruct, bool importanceSampling, bool OnlyLastEnergies){
    m_LearningRate = learningRate;
    m_stepLength = timestep;
    double decayRate = 0.0;// learningRate / numberOfOptimizationSteps;
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    int NrParameters = m_waveFunction->getNrParameters();
    int Nra = m_waveFunction->getNrVisibleNodes();
    int Nrb = m_waveFunction->getNrHiddenNodes();
    int NrH = Nra * Nrb;

    bool acceptedStep;
    int NrEqillibration_steps = (int)(numberOfMetropolisSteps * m_equilibrationFraction);

    for (int optStep = 0; optStep < numberOfOptimizationSteps; optStep++) {
        setSampler(new Sampler(this));
        std::vector<double>  cumulativeEnergyDiffWFProduct(NrParameters, 0);
        std::vector<double>  cumulativeDifferentiatedWF(NrParameters, 0);

        for (int i = 0; i < NrEqillibration_steps; i++) {
            if (importanceSampling){
                acceptedStep = MALAStep(timestep);
            }
            else {
                acceptedStep = metropolisStep();
            }
        }
        for (int i = 0; i < m_numberOfMetropolisSteps; i++) {
            if (importanceSampling) {
                acceptedStep = MALAStep(timestep);
            }
            else {
                acceptedStep = metropolisStep();
            }
            m_sampler->sampleGD(acceptedStep, NrParameters,cumulativeEnergyDiffWFProduct, cumulativeDifferentiatedWF);
        }

        m_sampler->computeAveragesGD(NrParameters, cumulativeEnergyDiffWFProduct, cumulativeDifferentiatedWF);

        for (int i = 0; i < NrParameters; i++) {
            cumulativeEnergyDiffWFProduct[i] = 2 * (cumulativeEnergyDiffWFProduct[i] - cumulativeDifferentiatedWF[i]+1.0e-9);
        }

        m_waveFunction->OptimizeParameters(cumulativeEnergyDiffWFProduct);
        //setLearningRate(m_LearningRate - decayRate);

        //for (int i = 0;i< cumulativeEnergyDiffWFProduct.size();i++)
        //    std::cout << "weight_gradient: " << cumulativeEnergyDiffWFProduct[i] << "\n";

        //m_sampler->printOutputToTerminal();
        //m_sampler->resetSampler();
        if (!OnlyLastEnergies) {
            datastruct.push_back({ "step:" + std::to_string(optStep),m_sampler->getEnergyVector() });
        }
        else {
            if (optStep == numberOfOptimizationSteps - 1) {
                datastruct.push_back({ std::to_string(learningRate),m_sampler->getEnergyVector() });
            }
        }
    }
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setDriftCoefficient(double coeff) {
    assert(coeff >= 0);
    m_driftCoefficient = coeff;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

void System::setFileOptString(std::string fileOptString) {
    m_fileOptString = fileOptString;
}

void System::setSampler(Sampler* sampler) {
    m_sampler = sampler;
}

void System::setElapsedTime(double elapsed) {
    m_elapsed = elapsed;
}

void System::setLearningRate(double LearningRate){
    m_LearningRate = LearningRate;
}

void System::setTimeStep(double timeStep) {
    m_timeStep = timeStep;
    m_timeStepSqrt = sqrt(timeStep);
}