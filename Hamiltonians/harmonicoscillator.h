#pragma once
#include "hamiltonian.h"
#include <vector>
#include <string>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(WaveFunction* wavefunction) override;
    std::string getName() override { return m_name; };
    double computeKE(WaveFunction* wavefunction) override;
    double computePE(WaveFunction* wavefunction) override;
    double computeIE(WaveFunction* wavefunction) override;

private:
    double m_omega = 0;
    double m_omegaSquared = 0;
    std::string m_name = "N";
};

