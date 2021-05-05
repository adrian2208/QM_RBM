#include "simplegaussian.h"
#include "simplegaussian_numerical.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"


double SimpleGaussian_Numerical::computeDoubleDerivative(std::vector<class Particle*> particles) {
    int NrParticles = m_system->getNumberOfParticles();
    int NrDimensions = m_system->getNumberOfDimensions();
    float stepSize = 0.0001;
    float twice_stepsize = 2 * stepSize;
    float Laplacian = 0;
    float particleLaplacian;
    for (int particleIndex = 0; particleIndex < NrParticles; particleIndex++) {
        particleLaplacian = 0;
        for (int i = 0; i < NrDimensions; i++) {
            particles[particleIndex]->adjustPosition(twice_stepsize, i);
            particleLaplacian += evaluate(particles, particleIndex);

            particles[particleIndex]->adjustPosition(-stepSize, i);
            particleLaplacian -= 2 * evaluate(particles, particleIndex);

            particles[particleIndex]->adjustPosition(-stepSize, i);
        }
        particleLaplacian /= evaluate(particles, particleIndex);
        particleLaplacian += NrDimensions;
        Laplacian += particleLaplacian;
    }
    return (Laplacian / (stepSize * stepSize));
}