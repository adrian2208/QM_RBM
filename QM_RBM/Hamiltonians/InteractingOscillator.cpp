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
    assert(omega > 0 && NrParticles == 2);
    m_omega = omega;
    m_omegaSquared = omega * omega;
}

double InteractingOscillator::computeLocalEnergy(WaveFunction* wavefunction) {
    std::vector<double> x = wavefunction->getX();
    int NrParticles = wavefunction->getNrParticles();
    int NrDimensions = wavefunction->getNrDimensions();
    int NrVisibleNodes = wavefunction->getNrVisibleNodes();
    double output = 0;
    double KE = 0;
    double ParticleDistance = 0;
    int index;

    for (int i = 0; i < NrVisibleNodes; i++) {
        output += pow(x[i], 2);
    }
    output *= 0.5 * m_omegaSquared;
    //std::cout << "PE:" << std::to_string(output);
    KE = wavefunction->computeDoubleDerivative();
    //std::cout << "    KE:" << std::to_string(KE);
    output += KE;
    

    //for (int i = 0; i < NrDimensions; i++) {
    //    index = i + NrDimensions;
    //    ParticleDistance += pow(x[i] - x[index], 2);
    //}
    
    //ParticleDistance = sqrt(pow(x[0], 2) + pow(x[1], 2)) - sqrt(pow(x[2], 2) + pow(x[3], 2));
    double interactionTerm = 0;
    double rad;
    for (int p = 0; p < NrVisibleNodes - NrDimensions; p += NrDimensions) {
        for (int s = (p + NrDimensions); s < NrVisibleNodes; s += NrDimensions) {
            rad = 0;
            for (int i = 0; i < NrDimensions; i++) {
                rad += pow(x[p + i] - x[s + i],2);
            }
            interactionTerm += 1.0 / sqrt(rad);
        }
    }
    //std::cout << "    Int:" << std::to_string(interactionTerm) << "\n";
    //introducing a small offset to improve numerical stability
    output += interactionTerm;//1.0 / (abs(ParticleDistance)+1.0e-8);
    return output;
}
