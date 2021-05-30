#pragma once
#include <string>
#include "wavefunction.h"

class qnet : public WaveFunction {
public:
    qnet(class System* system, int NrDimensions, int NrParticles, int NrHiddenNodes, double sigma,double sigma0);
    double evaluate() override;
    double evaluate(int particleIndex) override;
    double computeDoubleDerivative() override;
    std::vector<double> computeVariationalDerivative() override;
    std::vector<double> driftTerm(int indx) override;
    virtual std::vector<double> getParticlePosition(int ParticleIndx) override;

    std::string getName() { return m_name; };
    void adjustPosition(int node, double dx) override;
    void setPosition(int node, double x) override;
    void OptimizeParameters(std::vector<double>& ParameterGradientVector) override;

private:
    std::string m_name = "QRBM";
};
