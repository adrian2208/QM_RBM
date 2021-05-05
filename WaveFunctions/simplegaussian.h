#pragma once
#include <string>
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles) override;
    double evaluate(std::vector<class Particle*> particles, int particleIndex) override;
    double computeDoubleDerivative(std::vector<class Particle*> particles) override;
    double computeVariationalDerivative(std::vector<class Particle*> particles) override;
    double localEnergyDerivative() override;
    std::vector<double> driftTerm(std::vector<class Particle*> particles, int indx) override;
    std::string getName() { return m_name; };
private:
    std::string m_name = "SimpleGaussian";
};
