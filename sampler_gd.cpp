#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include "sampler.h"
#include "sampler_gd.h"

sampler_gd::sampler_gd(System* system) : Sampler(system) {};

void sampler_gd::sample(bool acceptedStep) {
    // Making sure the sampling variables are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergyDiffWFProduct = 0;
        m_cumulativeDifferentiatedWF = 0;
        m_NrAcceptedSteps = 0;
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    double localEnergy = m_system->getHamiltonian()->
        computeLocalEnergy(m_system->getParticles());
    double DiffWF = m_system->getWaveFunction()->computeVariationalDerivative(m_system->getParticles());
    m_cumulativeEnergy += localEnergy;
    m_cumulativeDifferentiatedWF += DiffWF;
    m_cumulativeEnergyDiffWFProduct += DiffWF * localEnergy;
    if (acceptedStep) {
        m_NrAcceptedSteps += 1;
    }
    m_stepNumber++;
}

void sampler_gd::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */
    double NrSteps = m_system->getNumberOfMetropolisSteps();
    m_energy = m_cumulativeEnergy / NrSteps;
    m_expectationOfProduct = m_cumulativeEnergyDiffWFProduct / NrSteps;
    m_productOfExpectations = m_cumulativeDifferentiatedWF * m_energy / NrSteps;
}