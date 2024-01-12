#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include "sampler.h"
#include "system.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    // Making sure the sampling variables are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0.0;
        m_NrAcceptedSteps = 0;
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    double localEnergy = m_system->getHamiltonian()->computeLocalEnergy(m_system->getWaveFunction());
    if (m_samplingPos) {
        samplePos();
    }
    

    m_cumulativeEnergy  += localEnergy;
    m_NrAcceptedSteps += (int)acceptedStep;
    m_energyVector.push_back(localEnergy);
    m_stepNumber++;
}

void Sampler::sampleGD(bool acceptedStep,int NrParameters,std::vector<double> &cumulativeEnergyDiffWFProduct, std::vector<double> &cumulativeDifferentiatedWF){

    //double localEnergy = m_system->getHamiltonian()->computeLocalEnergy(m_system->getWaveFunction());
    //m_cumulativeEnergy += localEnergy;
    double KE = m_system->getHamiltonian()->computeKE(m_system->getWaveFunction());
    double PE = m_system->getHamiltonian()->computePE(m_system->getWaveFunction());
    double IE = m_system->getHamiltonian()->computeIE(m_system->getWaveFunction());
    double localEnergy = KE + PE + IE;

    m_cumulativeKE += KE;
    m_cumulativePE += PE;
    m_cumulativeIE += IE;
    m_cumulativeEnergy += localEnergy;
    std::vector<double> DiffWF = m_system->getWaveFunction()->computeVariationalDerivative();
    
    m_energyVector.push_back(localEnergy);

    for (int i = 0; i < NrParameters; i++) {
        cumulativeDifferentiatedWF[i] += DiffWF[i];
        cumulativeEnergyDiffWFProduct[i] += DiffWF[i] * localEnergy;
    }
    if (acceptedStep) {
        m_NrAcceptedSteps += 1;
    }
    m_stepNumber++;
}



void Sampler::samplePos() {
    int NrParticles = m_system->getWaveFunction()->getNrParticles();
    int x_index;
    int y_index;
    int r_index;
    for (int i = 0; i < NrParticles;i++) {
        std::vector<double> particlePos = m_system->getWaveFunction()->getParticlePosition(i);
        if (abs(particlePos[0]) < m_intPosSamplingRadius) {
            if (abs(particlePos[1]) < m_intPosSamplingRadius) {
                x_index = floor((particlePos[0] + m_intPosSamplingRadius) / m_PosSamplingWidth);
                y_index = floor((particlePos[1] + m_intPosSamplingRadius) / m_PosSamplingWidth);
                r_index = floor(sqrt(pow(particlePos[0], 2) + pow(particlePos[0], 2)) /(2* m_PosSamplingWidth));
                m_particlePos_matrix[x_index][y_index] ++;
                m_particleRadius[r_index] ++;
                //std::cout << r_index << ",";
            }
        }

    }
        
}

void Sampler::SetupPositionSampling(int NrSamplingLengths, int intPosSamplingRadius) {
    m_samplingPos = true;
    m_particleRadius.resize(NrSamplingLengths);
    m_PosSamplingWidth = 2.0 * intPosSamplingRadius / NrSamplingLengths;
    //std::cout << "Width: " << m_PosSamplingWidth << "input 1 and 2: " << NrSamplingLengths << "   " << intPosSamplingRadius << "\n";
    m_NrSamplingLengths = NrSamplingLengths;
    m_intPosSamplingRadius = intPosSamplingRadius;
    m_particlePos_matrix = new int* [NrSamplingLengths];
    for (int i = 0; i < NrSamplingLengths; i++) {
        m_particlePos_matrix[i] = new int[NrSamplingLengths];
        m_particleRadius[i] = 0;
    }
    for (int i = 0; i < NrSamplingLengths; i++) {
        for (int j = 0; j < NrSamplingLengths; j++) {
            m_particlePos_matrix[i][j] = 0.0;
        }
    }

}
void Sampler::resetSampler(){
    m_stepNumber = 0;
    m_cumulativeEnergy = 0.0;
    m_NrAcceptedSteps = 0;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getWaveFunction()->getNrParticles();
    int     nd = m_system->getWaveFunction()->getNrDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    double  ef = m_system->getEquilibrationFraction();
    //std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms * ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << endl;
}

std::string Sampler::GenerateFileName(std::string descriptor) {
    int     np = m_system->getWaveFunction()->getNrParticles();
    int     nd = m_system->getWaveFunction()->getNrDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    double  ef = m_system->getEquilibrationFraction();
    std::string hamiltonian = m_system->getHamiltonian()->getName();
    std::string samplingType = m_system->getSamplingType();
    std::string out = "D" + std::to_string(nd) + "_P_" + std::to_string(np) + "I_" + hamiltonian + "_" + samplingType +"_S_2pow" + std::to_string((int)floor(std::log2(ms))) + "_eqS_2pow" + std::to_string((int)floor(std::log2(std::round(ms * ef)))) + "_" + descriptor;
    return out;
}

void Sampler::printOutputToFile() {

}



void Sampler::computeAverages() {
    m_energy = m_cumulativeEnergy / m_numberOfMetropolisSteps;
   
}

void Sampler::computeAveragesGD(int NrParameters, std::vector<double>& cumulativeEnergyDiffWFProduct, std::vector<double>& cumulativeDifferentiatedWF){
    double NrSteps = m_system->getNumberOfMetropolisSteps();
    m_energy = m_cumulativeEnergy / NrSteps;
    m_Kenergy = m_cumulativeKE / NrSteps;
    m_Penergy = m_cumulativePE / NrSteps;
    m_Ienergy = m_cumulativeIE / NrSteps;
    m_acceptedRatio =  m_NrAcceptedSteps/ NrSteps;
    for (int i = 0; i < NrParameters; i++) {
        cumulativeEnergyDiffWFProduct[i] = cumulativeEnergyDiffWFProduct[i] / NrSteps;//m_expectationOfProduct
        cumulativeDifferentiatedWF[i] = cumulativeDifferentiatedWF[i] * m_energy / NrSteps;//m_productOfExpectations
    }
}
