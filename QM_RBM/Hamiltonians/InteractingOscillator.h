#pragma once
#include "hamiltonian.h"
#include <vector>
#include <string>

class InteractingOscillator : public Hamiltonian {
public:
    InteractingOscillator(System* system, double omega,int NrParticles);
    double computeLocalEnergy(WaveFunction* wavefunction) override;
    std::string getName() override { return m_name; };


private:
    double m_omega = 0;
    double m_omegaSquared = 0;
    std::string m_name = "Y";
};

