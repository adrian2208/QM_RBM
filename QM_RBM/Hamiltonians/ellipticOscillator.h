#pragma once
#include "hamiltonian.h"
#include <vector>
#include <string>

class EllipticOscillator : public Hamiltonian {
public:
    EllipticOscillator(System* system, double gamma, double hardCoreDiameter, int NrDimensions);
    double computeLocalEnergy(std::vector<Particle*> particles) override;
    std::string getName() override { return m_name; };


private:
    double m_hardCoreDiameter = 0;
    double m_gamma = 0;
    double m_gamma_squared = 0;
    double m_hardCoreDiameter_squared = 0;
    std::string m_name = "EllipticOscillator";
};
