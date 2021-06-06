#pragma once
#include <vector>
#include <string>
#include "system.h"


class WaveFunction {
public:
    WaveFunction(class System* system);

    virtual double evaluate() = 0;
    virtual double evaluate(int ParticleIndx) = 0;

    virtual std::vector<double> computeVariationalDerivative() = 0;
    virtual double computeDoubleDerivative() = 0;
    virtual std::vector<double> driftTerm(int indx) =0;
    virtual std::string getName() = 0;
    virtual void adjustPosition(int node, double dx) = 0;
    virtual void setPosition(int node, double x) = 0;
    virtual void initializePositions() = 0;
    virtual void OptimizeParameters(std::vector<double>& ParameterGradientVector)=0;
    virtual void setParameters(std::vector<double>& ParameterVector) = 0;
    virtual void setOptimizer(class Optimizer* optimizer)=0;
    virtual std::vector<double> getParticlePosition(int ParticleIndx) = 0;
    virtual std::vector<double> getParameters() = 0;
    int getNrVisibleNodes();
    int getNrHiddenNodes();
    int getNrDimensions();
    int getNrParticles();
    int getNrParameters();
    std::vector<double> getX();


    double** m_W = nullptr;
    std::vector<double> m_X = std::vector<double>();
    //std::vector<int> m_H = std::vector<int>();
    std::vector<double> m_a = std::vector<double>();
    std::vector<double> m_b = std::vector<double>();



protected:
    int m_NrVisibleNodes = 0;
    int m_NrHiddenNodes = 0;
    int m_NrParameters = 0;
    int m_NrDimensions = 0;
    int m_NrParticles = 0;
    double m_sigma = 0.0;
    double m_sigmaSquared = 0.0;
    double m_oneOverSigmaSquared = 0.0;

    class System* m_system = nullptr;
};
