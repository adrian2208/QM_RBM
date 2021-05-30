#include "qnet.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <omp.h>
#include "Benchmark.h"
#include "Math/random.h"
#include "wavefunction.h"
#include "../system.h"
#include "../sampler.h"

qnet::qnet(System* system,
            int NrDimensions,
            int NrParticles, 
            int NrHiddenNodes, 
            double sigma,
            double sigma0) :
        WaveFunction(system) {

    assert(NrDimensions > 0 && NrParticles > 0 && NrHiddenNodes > 0);

    m_NrDimensions = NrDimensions;
    m_NrParticles = NrParticles;
    m_NrVisibleNodes = NrDimensions * NrParticles; 
    m_sigma = sigma;
    m_sigmaSquared = sigma * sigma;
    m_oneOverSigmaSquared = 1.0 / m_sigmaSquared;
    m_NrHiddenNodes = NrHiddenNodes;
    m_NrParameters = m_NrVisibleNodes + m_NrHiddenNodes + m_NrVisibleNodes * m_NrHiddenNodes;
    

    m_W = new double* [m_NrVisibleNodes];
    
    for (int i = 0; i < m_NrVisibleNodes; i++) {
        m_W[i] = new double[NrHiddenNodes];
        for (int j = 0; j < NrHiddenNodes; j++) {
            m_W[i][j] = m_system->getRandomEngine()->nextGaussian(0, sigma0);
        }

        m_X.push_back(m_system->getRandomEngine()->nextDouble()-0.5);
        m_a.push_back(m_system->getRandomEngine()->nextGaussian(0, sigma0));
    }
    for (int j = 0; j < NrHiddenNodes; j++ ) {

        m_b.push_back(m_system->getRandomEngine()->nextGaussian(0, sigma0));
        m_H.push_back(m_system->getRandomEngine()->nextInt(1));
    }
}

double qnet::evaluate() {
    double exponent = 0;
    std::vector<double> exponent2_sum(m_NrHiddenNodes);
    double product = 1;
//#pragma omp parallel for reduction(+:exponent)
        for (int i = 0; i < m_NrVisibleNodes; i++) {
            exponent += pow(m_X[i] - m_a[i], 2);
        }

    exponent = -exponent / (2 * m_sigmaSquared);

        int i, j;
        //#pragma omp parallel for private(i)
        for (j = 0; j < m_NrHiddenNodes; j++) {
            for (i = 0; i < m_NrVisibleNodes; i++) {
                exponent2_sum[j] = m_X[i] * m_W[i][j];
            }
        }
//#pragma omp parallel for reduction(*:product)
    for (int j = 0; j < m_NrHiddenNodes; j++) {
        product *= 1 + exp(m_b[j] + m_oneOverSigmaSquared * exponent2_sum[j]);
    }

    return exp(exponent)*product;
}

double qnet::evaluate(int particleIndex){
    int startIndx = particleIndex * m_NrDimensions;
    int endIndx = startIndx + m_NrDimensions;
    double exponent = 0;
    std::vector<double> exponent2_sum(m_NrHiddenNodes);
    double product = 1;
    //#pragma omp parallel for reduction(+:exponent)
    for (int i = startIndx; i < endIndx; i++) {
        exponent += pow(m_X[i] - m_a[i], 2);
    }

    exponent = -exponent / (2 * m_sigmaSquared);

    int i, j;
    //#pragma omp parallel for private(i)
    for (j = 0; j < m_NrHiddenNodes; j++) {
        for (i = 0; i < m_NrVisibleNodes; i++) {
            exponent2_sum[j] = m_X[i] * m_W[i][j];
        }
    }
    //#pragma omp parallel for reduction(*:product)
    for (int j = 0; j < m_NrHiddenNodes; j++) {
        product *= 1 + exp(m_b[j] + m_oneOverSigmaSquared * exponent2_sum[j]);
    }

    return exp(exponent) * product;
}

void qnet::adjustPosition(int node, double dx) {
    m_X[node] += dx;
}

void qnet::setPosition(int node, double x){
    m_X[node] = x;
}

void qnet::OptimizeParameters(std::vector<double>& ParameterGradientVector){
    double Lr = m_system->getLearningRate();
    int idx = 0;
    for (int i = 0; i < m_NrVisibleNodes; i++) {
        m_a[i] = m_a[i]-Lr*ParameterGradientVector[idx];
        //std::cout << "a = " << m_a[i] << "\n";
        idx++;
     }
    for (int i = 0; i < m_NrHiddenNodes; i++) {
        m_b[i] = m_b[i]-Lr* ParameterGradientVector[idx];
        //std::cout << "b = " << m_b[i] << "\n";
        idx++;
    }

    for (int i = 0; i < m_NrVisibleNodes; i++) {
        for (int j = 0; j < m_NrHiddenNodes; j++) {
            m_W[i][j] = m_W[i][j] - Lr* ParameterGradientVector[idx];
            idx++;
        }
    }
}



double qnet::computeDoubleDerivative() {
    /// <summary>
    /// Computes Equation XX.XX 
    /// </summary>
    /// <returns>kinetic energy</returns>
    std::vector<double> Q(m_NrVisibleNodes);
    double xWSum = 0;//needed for calculating k
    double sigmaWSum = 0;
    double sigmaWSquaredSum = 0;
    double output = 0;
    std::vector<double> sigma(m_NrHiddenNodes);

    // calculating sigma
    for (int j = 0; j < m_NrHiddenNodes; j++) {
        for (int i = 0; i < m_NrVisibleNodes; i++) {
            xWSum += m_X[i] * m_W[i][j];
        }
        sigma[j] = 1.0 / (1 + exp(-m_b[j] - m_oneOverSigmaSquared * xWSum));
        xWSum = 0;
    }

    // calculating sums and output
    for (int i = 0; i < m_NrVisibleNodes; i++) {
        for (int j = 0; j < m_NrHiddenNodes; j++) {
            sigmaWSum += m_W[i][j] * sigma[j];
            sigmaWSquaredSum += sigmaWSum * m_W[i][j] * (1 - sigma[j]);
        }
        output += pow(m_X[i] - m_a[i]-sigmaWSum,2)+ sigmaWSquaredSum;
    }
    output *= -0.5 * m_oneOverSigmaSquared * m_oneOverSigmaSquared;
    output += 0.5 * m_NrVisibleNodes * m_oneOverSigmaSquared;

    return output;
}


std::vector<double> qnet::computeVariationalDerivative() {
    int middleIndx = m_NrVisibleNodes + m_NrHiddenNodes;
    int endIndx = m_NrParameters;
    std::vector<double> gradientVector(endIndx);
    std::vector<double> exponentSum(m_NrHiddenNodes);

    double xWSum = 0;
    std::vector<double> sigma(m_NrHiddenNodes);

    // calculating sigmoid
    for (int j = 0; j < m_NrHiddenNodes; j++) {
        for (int i = 0; i < m_NrVisibleNodes; i++) {
            xWSum += m_X[i] * m_W[i][j];
        }
        sigma[j] = 1.0 / (1 + exp(-m_b[j] - m_oneOverSigmaSquared * xWSum));
        xWSum = 0;
    }
    
    int i = 0;
    for (; i < m_NrVisibleNodes; i++) {
        gradientVector[i] = (m_X[i] - m_a[i]) * m_oneOverSigmaSquared;
    }

    int j = 0;
    for (; i < middleIndx; i++) {
        gradientVector[i] = sigma[j];
        j++;
    }

    for (int k = 0; k < m_NrVisibleNodes; k++) {
        for (int l = 0; l < m_NrHiddenNodes; l++) {
            gradientVector[i] = m_X[k] * m_oneOverSigmaSquared * sigma[l];
            i++;
        }
    }

    return gradientVector;
}



std::vector<double> qnet::driftTerm(int indx){

    double xWSum = 0;//needed for calculating k
    std::vector<double> sigma(m_NrHiddenNodes);

    // calculating sigmoid
    for (int j = 0; j < m_NrHiddenNodes; j++) {
        for (int i = 0; i < m_NrVisibleNodes; i++) {
            xWSum += m_X[i] * m_W[i][j];
        }
        sigma[j] = 1.0 / (1 + exp(-m_b[j] - m_oneOverSigmaSquared * xWSum));
        xWSum = 0;
    }

    std::vector<double> driftVector(m_NrDimensions,0);
    int vector_index = 0;
    for (int i = indx;i< indx+m_NrDimensions;i++){
        for (int j = 0; j < m_NrHiddenNodes; j++) {
            driftVector[vector_index] += sigma[j] * m_W[i][j];
        }
        driftVector[vector_index] -= m_X[i] - m_a[i];
        driftVector[vector_index] *= 2 * m_oneOverSigmaSquared;
        vector_index++;

    }

    return driftVector;
}

std::vector<double> qnet::getParticlePosition(int ParticleIndx){
    int startingIndx = ParticleIndx * m_NrDimensions;
    std::vector<double> pos(m_NrDimensions);
    for (int i = 0; i < m_NrDimensions; i++) {
        pos[i] = m_X[startingIndx + i];
    }

    return pos;
}


