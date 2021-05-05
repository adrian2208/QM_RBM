#pragma once
#include "hamiltonian.h"
#include <vector>
#include <string>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles) override;
    std::string getName() override { return m_name; };


private:
    double m_omega = 0;
    std::string m_name = "HarmonicOscillator";
};

