#pragma once
#include<string>

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void sampleGD(bool acceptedStep,int NrParameters,std::vector<double>& cumulativeEnergyDiffWFProduct, std::vector<double>& cumulativeDifferentiatedWF);
    
    void computeAverages();
    void computeAveragesGD(int NrParameters, std::vector<double>& cumulativeEnergyDiffWFProduct, std::vector<double>& cumulativeDifferentiatedWF);

    void printOutputToTerminal();
    void printOutputToFile();
    void SetupPositionSampling(int NrSamplingLengths, int intPosSamplingRadius);
    void resetSampler();

    double getEnergy()          { return m_energy; }
    double getNumberAccepted() { return m_NrAcceptedSteps; }
    int**  getParticlePos_matrix() { return m_particlePos_matrix; }
    std::vector<double> getEnergyVector() { return m_energyVector; }

    std::string GenerateFileName(std::string descriptor);

protected:
    // Position Sampling
    void    samplePos();
    int**   m_particlePos_matrix            = nullptr;
    int     m_intPosSamplingRadius          = 0;
    int     m_NrSamplingLengths             = 0;
    double  m_PosSamplingWidth              = 0;
    bool    m_samplingPos                   = false;

    // sampling
    int     m_numberOfMetropolisSteps       = 0;
    int     m_stepNumber                    = 0;
    double  m_energy                        = 0;
    double  m_cumulativeEnergy              = 0;
    int     m_NrAcceptedSteps               = 0;


    std::vector<double> m_energyVector;
    
    class System* m_system = nullptr;
};
