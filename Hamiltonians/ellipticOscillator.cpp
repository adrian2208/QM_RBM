#include "ellipticOscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

EllipticOscillator::EllipticOscillator(System* system, double gamma, double hardCoreDiameter, int NrDimensions) :
    Hamiltonian(system) {
    assert(gamma > 0 && hardCoreDiameter > 0 && NrDimensions == 3);
    m_gamma = gamma;
    m_gamma_squared = gamma*gamma;
    m_hardCoreDiameter = hardCoreDiameter;
    m_hardCoreDiameter_squared = hardCoreDiameter* hardCoreDiameter;
}

double EllipticOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    double potentialEnergy = 0;
    double kineticEnergy = -0.5 * m_system->getWaveFunction()->computeDoubleDerivative(particles);
    int numberOfParticles = m_system->getNumberOfParticles();
    std::vector<double> particlePosition(3);
    std::vector<double> particle2Position(3);
    double particleDistanceSquared;
    double localEnergy = 0;
    for (int i = 0; i < numberOfParticles; i++) {
        particlePosition = particles[i]->getPosition();
        potentialEnergy += particlePosition[0] * particlePosition[0] + particlePosition[1] * particlePosition[1] + m_gamma_squared * particlePosition[2] * particlePosition[2];

        for (int j = i + 1; j < numberOfParticles; j++) {
            particle2Position = particles[j]->getPosition();
            particleDistanceSquared = 
              pow(particlePosition[0] - particle2Position[0], 2)
            + pow(particlePosition[1] - particle2Position[1], 2)
            + pow(particlePosition[2] - particle2Position[2], 2);
            if (particleDistanceSquared <= m_hardCoreDiameter_squared) {
                localEnergy += 1000; //Big number- approximating infinity
                std::cout << "Crash!!!!\n";
            }
        }

    }
    localEnergy += kineticEnergy + potentialEnergy * 0.5;
    return localEnergy;
}