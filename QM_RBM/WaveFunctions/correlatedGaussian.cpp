#include "correlatedGaussian.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "../sampler.h"

CorrelatedGaussian::CorrelatedGaussian(class System* system, double alpha, double beta, double hardCoreDiameter) : WaveFunction(system) {
    m_hardCoreDiameter_squared = hardCoreDiameter* hardCoreDiameter;
    m_hardCoreDiameter = hardCoreDiameter;
    m_numberOfParameters = 2;

    assert(alpha >= 0 && beta >=0);
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
}

double CorrelatedGaussian::evaluateGaussian(std::vector<class Particle*> particles) {
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    std::vector<double> particlePostition(3);
    double totalPosSquared = 0;
    for (auto particle : particles) {
        particlePostition = particle->getPosition();
        totalPosSquared +=
            pow(particlePostition[0], 2) +
            pow(particlePostition[1], 2) +
            beta * pow(particlePostition[2], 2);
    }

    return exp(-alpha * totalPosSquared);
}

double CorrelatedGaussian::evaluate(std::vector<class Particle*> particles) {
    double correlation= 1;
    int NrParticles = m_system->getNumberOfParticles();
    for (int i = 0; i < NrParticles; i++) {
        for (int j = i + 1; j < NrParticles; j++) {
            correlation *= correlationFunction(particles[i], particles[j]);
        }
    }

    return evaluateGaussian(particles) * correlation;
}

double CorrelatedGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    double alpha = m_parameters[0];
    double beta = m_parameters[1];
    
    int NrParticles = m_system->getNumberOfParticles();
    double output = 0;
    for (int particleIndex = 0; particleIndex < NrParticles; particleIndex++) {
        std::vector<double> particle1Position = particles[particleIndex]->getPosition();

        double x_squared = pow(particle1Position[0], 2);
        double y_squared = pow(particle1Position[1], 2);
        double z_squared = pow(particle1Position[2], 2);
        output += 4 * alpha * alpha * (x_squared + 
            y_squared + beta * beta * z_squared) -  2 * alpha * (2 + beta);

        double distance_squared;
        double correlationPrefactor;
        std::vector<double> phiGradient(3);
        phiGradient[0] = -4 * m_hardCoreDiameter * x_squared;
        phiGradient[1] = -4 * m_hardCoreDiameter * y_squared;
        phiGradient[2] = -4 * m_hardCoreDiameter * beta * z_squared;
        std::vector<double> particle2Position(3);
        std::vector<double> secondTermSum(3, 0);

        for (int i = 0; i < particleIndex; i++) {
            particle2Position = particles[i]->getPosition();
            distance_squared = P2PDistanceSquared(particles[particleIndex], particles[i]);

            if (distance_squared > m_hardCoreDiameter_squared) {
                correlationPrefactor = -m_hardCoreDiameter / (distance_squared * (m_hardCoreDiameter - sqrt(distance_squared)));
                secondTermSum = vectorAddition(secondTermSum, vectorDifference_scaler(particle1Position, particle2Position, correlationPrefactor));
                output += (m_hardCoreDiameter_squared - 2 * m_hardCoreDiameter * sqrt(distance_squared)) / (distance_squared * pow(sqrt(distance_squared) - m_hardCoreDiameter, 2)) + correlationPrefactor;
            }
        }

        for (int i = particleIndex + 1; i < NrParticles; i++) {
            particle2Position = particles[i]->getPosition();
            distance_squared = P2PDistanceSquared(particles[particleIndex], particles[i]);

            if (distance_squared > m_hardCoreDiameter_squared) {
                correlationPrefactor = -m_hardCoreDiameter / (distance_squared * (m_hardCoreDiameter - sqrt(distance_squared)));
                secondTermSum = vectorAddition(secondTermSum, vectorDifference_scaler(particle1Position, particle2Position, correlationPrefactor));
                output += (m_hardCoreDiameter_squared - 2 * m_hardCoreDiameter * sqrt(distance_squared)) / (distance_squared * pow(sqrt(distance_squared) - m_hardCoreDiameter, 2)) + correlationPrefactor;
            }
        }

        output += scalarProduct(secondTermSum, phiGradient);
        output += scalarProduct(secondTermSum, secondTermSum);
    }
    return output;

}

double CorrelatedGaussian::correlationFunction(Particle* particle1, Particle* particle2) {
    double distance = P2PDistanceSquared(particle1, particle2);
    if (distance > m_hardCoreDiameter_squared) {
            return 1 - m_hardCoreDiameter / sqrt(distance);
    }
    else {
        return 0;
    }
}

double CorrelatedGaussian::localEnergyDerivative() {
    double arg1 = m_system->getSampler()->getExpectationOfProduct();
    double arg2 = m_system->getSampler()->getProductOfExpectations();
    //std::cout << "arg1 = " << arg1 << "\n" << "arg2 = " << arg2 << "\n";

    return 2 * (arg1 - arg2);
}

std::vector<double> CorrelatedGaussian::driftTerm(std::vector<class Particle*> particles, int indx) {
    double alpha_times_minusfour = -4 * m_parameters[0];
    double beta = m_parameters[1];
    int NrParticles = m_system->getNumberOfParticles();

    std::vector<double> position = particles[indx]->getPosition();
    std::vector<double> forceVector(3, alpha_times_minusfour);
    forceVector[0] *= position[0];
    forceVector[1] *= position[1];
    forceVector[2] *= position[2]*beta;
    
    double distance_squared;
    double correlationPrefactor;
    std::vector<double> scaledDifferenceVector(3);

    for (int i = 0; i < indx; i++) {
        distance_squared = P2PDistanceSquared(particles[i], particles[indx]);
        correlationPrefactor = 2*m_hardCoreDiameter / (distance_squared * (m_hardCoreDiameter - sqrt(distance_squared)));
        scaledDifferenceVector = vectorDifference_scaler(particles[i]->getPosition(), position, correlationPrefactor);
        forceVector = vectorAddition(forceVector, scaledDifferenceVector);
    }
    for (int i = indx+1; i < NrParticles; i++) {
        distance_squared = P2PDistanceSquared(particles[i], particles[indx]);
        correlationPrefactor = 2 * m_hardCoreDiameter / (distance_squared * (m_hardCoreDiameter - sqrt(distance_squared)));
        scaledDifferenceVector = vectorDifference_scaler(particles[i]->getPosition(), position, correlationPrefactor);
        forceVector = vectorAddition(forceVector, scaledDifferenceVector);
    }

    return forceVector;
}

double CorrelatedGaussian::computeVariationalDerivative(std::vector<class Particle*> particles) {
    double sum_position_squared = 0;
    for (Particle* particle : particles) {
        sum_position_squared += particle->positionSquared();
    }

    //this output is already normalized with the wavefunction
    return -sum_position_squared;
}

double CorrelatedGaussian::P2PDistanceSquared(Particle* particle1, Particle* particle2) {
    std::vector<double> position1 = particle1->getPosition();
    std::vector<double> position2 = particle2->getPosition();
        
    return 
        pow(position1[0] - position2[0], 2) +
        pow(position1[1] - position2[1], 2) +
        pow(position1[2] - position2[2], 2);
}

std::vector<double> CorrelatedGaussian::vectorDifference_scaler(std::vector<double> vector1, std::vector<double> vector2, double prefactor) {
    std::vector<double> out(3);
    out[0] = prefactor * (vector1[0] - vector2[0]);
    out[1] = prefactor * (vector1[1] - vector2[1]);
    out[2] = prefactor * (vector1[2] - vector2[2]);
    return out;
 }

double CorrelatedGaussian::scalarProduct(std::vector<double> vector1, std::vector<double> vector2) {
    return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
}

std::vector<double> CorrelatedGaussian::vectorAddition(std::vector<double> vector1, std::vector<double> vector2) {
    std::vector<double> out(3);
    out[0] = vector1[0] + vector2[0];
    out[1] = vector1[1] + vector2[1];
    out[2] = vector1[2] + vector2[2];
    return out;
}