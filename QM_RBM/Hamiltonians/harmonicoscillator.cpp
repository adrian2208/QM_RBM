#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
           Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_omegaSquared = omega * omega;
}

double HarmonicOscillator::computeLocalEnergy(WaveFunction* wavefunction) {
    std::vector<double> x = wavefunction->getX();
    int NrParticles = wavefunction->getNrParticles();
    int NrDimensions = wavefunction->getNrDimensions();
    int NrVisibleNodes = wavefunction->getNrVisibleNodes();
    double output = 0;

    for (int i = 0; i < NrVisibleNodes; i++) {
        output += pow(x[i], 2);
    }
    output *= 0.5 * m_omegaSquared;
    output += wavefunction->computeDoubleDerivative();


    return output;
}

