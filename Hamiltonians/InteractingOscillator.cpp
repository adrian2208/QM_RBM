#include "InteractingOscillator.h"
#include <cassert>
#include <iostream>
#include <string>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

InteractingOscillator::InteractingOscillator(System* system, double omega, int NrParticles) :
    Hamiltonian(system) {
    assert(omega > 0.0 && NrParticles == 2);
    m_omega = omega;
    m_omegaSquared = omega * omega;
}

double InteractingOscillator::computeLocalEnergy(WaveFunction* wavefunction) {
    std::vector<double> x = wavefunction->getX();
    int NrParticles = wavefunction->getNrParticles();
    int NrDimensions = wavefunction->getNrDimensions();
    int NrVisibleNodes = wavefunction->getNrVisibleNodes();
    double output = 0.0;
    double KE = 0.0;
    double ParticleDistance = 0.0;
    int index;

    for (int i = 0; i < NrVisibleNodes; i++) {
        output += pow(x[i], 2);
    }
    output *= 0.5* m_omegaSquared;

    KE = wavefunction->computeDoubleDerivative();

    output += KE;
    


    double interactionTerm = 0.0;
    double rad = 0.0;
    int idx1, idx2;
    for (int p = 0; p < NrVisibleNodes - NrDimensions; p += NrDimensions) {
        for (int s = (p + NrDimensions); s < NrVisibleNodes; s += NrDimensions) {
            for (int i = 0; i < NrDimensions; i++) {
                idx1 = p + i;
                idx2 = s + i;
                rad += pow(x[idx1] - x[idx2],2.0);
            }
            interactionTerm += 1.0 / sqrt(rad);
            rad = 0.0;
        }
    }

    output += interactionTerm;
    return output;
}

double InteractingOscillator::computeKE(WaveFunction* wavefunction) {
    return wavefunction->computeDoubleDerivative();
}

double InteractingOscillator::computePE(WaveFunction* wavefunction) {
    std::vector<double> x = wavefunction->getX();
    int NrVisibleNodes = wavefunction->getNrVisibleNodes();
    double output = 0.0;

    for (int i = 0; i < NrVisibleNodes; i++) {
        output += pow(x[i], 2.0);
    }
    output *= 0.5 * m_omegaSquared;

    return output;
}

double InteractingOscillator::computeIE(WaveFunction* wavefunction) {
    std::vector<double> x = wavefunction->getX();
    int NrDimensions = wavefunction->getNrDimensions();
    int NrVisibleNodes = wavefunction->getNrVisibleNodes();

    double interactionTerm = 0.0;
    double rad;
    int idx1, idx2;
    for (int p = 0; p < NrVisibleNodes - NrDimensions; p += NrDimensions) {
        for (int s = (p + NrDimensions); s < NrVisibleNodes; s += NrDimensions) {
            rad = 0.0;
            for (int i = 0; i < NrDimensions; i++) {
                idx1 = p + i;
                idx2 = s + i;
                rad += pow(x[idx1] - x[idx2], 2.0);
            }
            interactionTerm += 1.0 / (sqrt(rad)+1.0e-10);
        }
    }
    return interactionTerm;
}