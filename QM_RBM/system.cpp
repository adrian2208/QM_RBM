#include "system.h"
#include <cassert>
#include <iostream>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    int particleIndex = m_random->nextInt(m_numberOfParticles - 1);
    Particle* randParticle = m_particles[particleIndex];
    std::vector<double> OldPos = randParticle->getPosition();
    double wfOld = m_waveFunction->evaluate(m_particles, particleIndex);
    
    for (int i = 0; i < m_numberOfDimensions; i++) {
        randParticle->adjustPosition(m_stepLength * (m_random->nextDouble()-0.5), i);
    }
    
    double wfNew = m_waveFunction->evaluate(m_particles, particleIndex);

    if (wfNew * wfNew / (wfOld * wfOld) < m_random->nextDouble()) {
        randParticle->setPosition(OldPos);
        return false;
    }
    return true;
}

bool System::MALAStep() {
    int particleIndex = m_random->nextInt(m_numberOfParticles - 1);
    Particle* randParticle = m_particles[particleIndex];
    std::vector<double> OldPos = m_particles[particleIndex]->getPosition();
    double wfOld = m_waveFunction->evaluate(m_particles, particleIndex);
    std::vector<double> driftTerm = m_waveFunction->driftTerm(m_particles, particleIndex);
    
    for (int i = 0; i < m_numberOfDimensions; i++) {
        randParticle->adjustPosition(m_driftCoefficient * driftTerm[i] * m_timeStep 
            + m_random->nextGaussian(0, 1) * m_timeStepSqrt, i);
    }

    std::vector<double> driftTermNew = m_waveFunction->driftTerm(m_particles, particleIndex);
    double Greens_exponent = 0;
    for (int i = 0; i < m_numberOfDimensions; i++) {
        Greens_exponent += 0.5 * (driftTerm[i] + driftTermNew[i]) *
            (m_driftCoefficient * m_timeStep * 0.5 * (driftTerm[i]
                - driftTermNew[i]) - randParticle->getPosition()[i] + OldPos[i]);
    }

    double wfNew = m_waveFunction->evaluate(m_particles, particleIndex);

    if (wfNew * wfNew* exp(Greens_exponent) / (wfOld * wfOld) < m_random->nextDouble()) {
        randParticle->setPosition(OldPos);
        return false;
    }
    return true;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    //m_sampler                   = new Sampler(this);
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

void System::runMALASteps(int numberOfMetropolisSteps) {
    m_particles = m_initialState->getParticles();
    //m_sampler = new Sampler(this);
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    bool acceptedStep;
    int NrEqillibration_steps = (int)(numberOfMetropolisSteps * m_equilibrationFraction);

    for (int i = 0; i < NrEqillibration_steps; i++) {
        acceptedStep = MALAStep();
    }
    for (int i = 0; i < m_numberOfMetropolisSteps; i++) {
        acceptedStep = MALAStep();
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

void System::setTimeStep(double timeStep) {
    m_timeStep = timeStep;
    m_timeStepSqrt = sqrt(timeStep);
}