#include "simplegaussian.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "../sampler.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
    double alpha = m_parameters[0];
    double sum_position_squared = 0;
    
    for (int i = 0; i < sizeof(particles); i++) {
        sum_position_squared += particles[i]->positionSquared();
    }

    return exp(-alpha * sum_position_squared);
}
/*To avoid wasting cpu cycles, we can evaluate only the factor of the wavefunction that changes 
*during a particle move.*/
double SimpleGaussian::evaluate(std::vector<class Particle*> particles, int particleIndex) {
    return exp(-m_parameters[0] * particles[particleIndex]->positionSquared());
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    double alpha = m_parameters[0];
    int NrParticles = m_system->getNumberOfParticles();
    int NrDimensions = m_system->getNumberOfDimensions();
    double sum_position_squared = 0;

    for (Particle *particle : particles) {
        sum_position_squared += particle->positionSquared();
    }

    double prefactor = -2 * alpha * NrParticles * NrDimensions + 4 * alpha * alpha * sum_position_squared;
   
    return prefactor;
}

std::vector<double> SimpleGaussian::driftTerm(std::vector<class Particle*> particles, int indx) {
    double alpha_times_minusfour = -4 * m_parameters[0];
    std::vector<double> output (m_system->getNumberOfDimensions());
    std::vector<double> position = particles[indx]->getPosition();

    for (double component : position) {
        output.push_back(alpha_times_minusfour * component);
    }
    return output;
}

double SimpleGaussian::computeVariationalDerivative(std::vector<class Particle*> particles) {
    double sum_position_squared = 0;
    for (Particle* particle : particles) {
        sum_position_squared += particle->positionSquared();
    }

    //this output is already normalized with the wavefunction
    return -sum_position_squared;
}

double SimpleGaussian::localEnergyDerivative() {
    double arg1 = m_system->getSampler()->getExpectationOfProduct();
    double arg2 = m_system->getSampler()->getProductOfExpectations();
    //std::cout << "arg1 = " << arg1 << "\n" << "arg2 = " << arg2 << "\n";
    
    return 2 * (arg1 - arg2);
}