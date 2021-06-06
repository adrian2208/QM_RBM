#pragma once
#include <vector>
#include <string>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(class WaveFunction* wavefunction) = 0;
    virtual std::string getName() = 0;
    virtual double computeKE(WaveFunction* wavefunction) = 0;
    virtual double computePE(WaveFunction* wavefunction) = 0;
    virtual double computeIE(WaveFunction* wavefunction) = 0;
protected:
    class System* m_system = nullptr;
};

