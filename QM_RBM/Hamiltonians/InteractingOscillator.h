#pragma once
#include "hamiltonian.h"
#include <vector>
#include <string>

class InteractingOscillator : public Hamiltonian {
public:
    InteractingOscillator(System* system, double omega,int NrParticles);
    double computeLocalEnergy(WaveFunction* wavefunction) override;
    std::string getName() override { return m_name; };
    
    double computeKE(WaveFunction* wavefunction);
    double computePE(WaveFunction* wavefunction);
    double computeIE(WaveFunction* wavefunction);

private:
    double m_omega = 0.0;
    double m_omegaSquared = 0.0;
    std::string m_name = "Y";
};

