#pragma once
#include <string>
#include "wavefunction.h"

class CorrelatedGaussian : public WaveFunction {
public:
    CorrelatedGaussian(class System* system, double alpha, double beta, double hardCoreDiameter);
    double evaluate(std::vector<class Particle*> particles) override;
    //double evaluate(std::vector<class Particle*> particles, int particleIndex) override;
    double evaluateGaussian(std::vector<class Particle*> particles);
    //double evaluateGaussian(std::vector<class Particle*> particles, int particleIndex);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeVariationalDerivative(std::vector<class Particle*> particles) override;
    double localEnergyDerivative() override;
    std::vector<double> driftTerm(std::vector<class Particle*> particles, int indx) override;
    std::string getName() { return m_name; };
private:
    double correlationFunction(Particle* particle1, Particle* particle2);
    double P2PDistanceSquared(Particle* particle1, Particle* particle2);
    double scalarProduct(std::vector<double> vector1, std::vector<double> vector2);
    std::vector<double> vectorDifference_scaler(std::vector<double> vector1, std::vector<double> vector2, double prefactor);
    std::vector<double> vectorAddition(std::vector<double> vector1, std::vector<double> vector2);

    std::string m_name = "CorrelatedGaussian";
    double m_hardCoreDiameter_squared = 0;
    double m_hardCoreDiameter = 0;
};
