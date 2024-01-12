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
#include "Optimizer/Optimizer.h"

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
            m_W[i][j] = m_system->getRandomEngine()->nextGaussian(0.0, sigma0);
        }

        m_X.push_back(m_system->getRandomEngine()->nextDouble()-0.5);
        m_a.push_back(m_system->getRandomEngine()->nextGaussian(0.0, sigma0));
    }
    for (int j = 0; j < NrHiddenNodes; j++ ) {

        m_b.push_back(m_system->getRandomEngine()->nextGaussian(0.0, sigma0));
        //m_H.push_back(m_system->getRandomEngine()->nextInt(1));
    }
}

void qnet::initializePositions() {
    for (int i = 0; i < m_NrVisibleNodes; i++) {
        m_X[i] = m_system->getRandomEngine()->nextDouble() - 0.5;
    }
}

double qnet::evaluate() {
    double exponent = 0.0;
    std::vector<double> exponent2_sum(m_NrHiddenNodes,0);
    double product = 1.0;

        for (int i = 0; i < m_NrVisibleNodes; i++) {
            exponent += pow(m_X[i] - m_a[i], 2.0);
        }

    exponent = -exponent / (2.0 * m_sigmaSquared);

        for (int j = 0; j < m_NrHiddenNodes; j++) {
            for (int i = 0; i < m_NrVisibleNodes; i++) {
                exponent2_sum[j] += m_X[i] * m_W[i][j];
            }
        }

    for (int j = 0; j < m_NrHiddenNodes; j++) {
        product *= 1.0 + exp(m_b[j] + m_oneOverSigmaSquared * exponent2_sum[j]);
    }

    return exp(exponent)*product;
}

double qnet::evaluate(int particleIndex){
    int startIndx = particleIndex * m_NrDimensions;
    int endIndx = startIndx + m_NrDimensions;
    double exponent = 0.0;
    std::vector<double> exponent2_sum(m_NrHiddenNodes,0);
    double product = 1.0;

    for (int i = startIndx; i < endIndx; i++) {
        exponent += pow(m_X[i] - m_a[i], 2);
    }

    exponent = -exponent / (2.0 * m_sigmaSquared);


    for (int j = 0; j < m_NrHiddenNodes; j++) {
        for (int i = 0; i < m_NrVisibleNodes; i++) {
            exponent2_sum[j] += m_X[i] * m_W[i][j];
        }
    }

    for (int j = 0; j < m_NrHiddenNodes; j++) {
        product *= (1.0 + exp(m_b[j] + m_oneOverSigmaSquared * exponent2_sum[j]));
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
    m_optimizer->optimize(m_a, m_b, m_W, ParameterGradientVector);
}

void qnet::setParameters(std::vector<double>& ParameterVector) {
    int idx = 0;
    for (int i = 0; i < m_NrVisibleNodes; i++) {
        m_a[i] = ParameterVector[idx];
        idx++;
    }
    for (int i = 0; i < m_NrHiddenNodes; i++) {
        m_b[i] = ParameterVector[idx];
        idx++;
    }

    for (int i = 0; i < m_NrVisibleNodes; i++) {
        for (int j = 0; j < m_NrHiddenNodes; j++) {
            m_W[i][j] = ParameterVector[idx];
            idx++;
        }
    }
}

std::vector<double> qnet::getParameters() {
    std::vector<double> output(m_NrParameters);
    int idx = 0;
    for (int i = 0; i < m_NrVisibleNodes; i++) {
        output[idx]= m_a[i];
        idx++;
    }
    for (int i = 0; i < m_NrHiddenNodes; i++) {
        output[idx] = m_b[i];
        idx++;
    }

    for (int i = 0; i < m_NrVisibleNodes; i++) {
        for (int j = 0; j < m_NrHiddenNodes; j++) {
            output[idx] = m_W[i][j];
            idx++;
        }
    }
    return output;
}

void qnet::setOptimizer(Optimizer* optimizer){
    m_optimizer = optimizer;
}





double qnet::computeDoubleDerivative() {
    double xWSum = 0.0;
    double sigmoidWSum = 0.0;
    double sigmoidWSquaredSum = 0.0;
    double output = 0.0;
    std::vector<double> sigmoid(m_NrHiddenNodes);

    // calculating sigmoid
    for (int j = 0; j < m_NrHiddenNodes; j++) {
        for (int k = 0; k < m_NrVisibleNodes; k++) {
            xWSum += m_X[k] * m_W[k][j];
        }
        sigmoid[j] = 1.0 / (1.0 + exp(-m_b[j] - m_oneOverSigmaSquared * xWSum));
        xWSum = 0.0;
    }
    double temp = 0.0;
    for (int i = 0; i < m_NrVisibleNodes; i++) {
        for (int j = 0; j < m_NrHiddenNodes; j++) {
            temp = m_W[i][j] * sigmoid[j];
            sigmoidWSum += temp;
            sigmoidWSquaredSum += temp * m_W[i][j] * (1.0 - sigmoid[j]);
        }
        output += pow(m_oneOverSigmaSquared * (m_a[i] - m_X[i] + sigmoidWSum), 2.0);
        output += -m_oneOverSigmaSquared + pow(m_oneOverSigmaSquared, 2.0) * sigmoidWSquaredSum;
        sigmoidWSum = 0.0;
        sigmoidWSquaredSum = 0.0;
    }

    output *= -0.5;


    return output;
}


std::vector<double> qnet::computeVariationalDerivative() {
    int middleIndx = m_NrVisibleNodes + m_NrHiddenNodes;
    int endIndx = m_NrParameters;
    std::vector<double> gradientVector(endIndx);

    double xWSum = 0.0;
    std::vector<double> sigma(m_NrHiddenNodes);

    // calculating sigmoid
    for (int j = 0; j < m_NrHiddenNodes; j++) {
        for (int i = 0; i < m_NrVisibleNodes; i++) {
            xWSum += m_X[i] * m_W[i][j];
        }
        sigma[j] = 1.0 / (1.0 + exp(-m_b[j] - m_oneOverSigmaSquared * xWSum));
        //std::cout << "sigmoid: " << sigma[j] << "\n";
        xWSum = 0.0;
    }
    
    int i = 0;
    for (; i < m_NrVisibleNodes; i++) {
        gradientVector[i] = (m_X[i] - m_a[i]) * m_oneOverSigmaSquared;
    }
    //std::cout << "d_da done: " << "index: " << i << "\n";
    int j = 0;
    for (; i < middleIndx; i++) {
        gradientVector[i] = sigma[j];
        j++;
    }
    //std::cout << "d_db done: " << "index: " << i << "\n";
    for (int k = 0; k < m_NrVisibleNodes; k++) {
        for (int l = 0; l < m_NrHiddenNodes; l++) {
            gradientVector[i] = m_X[k] * m_oneOverSigmaSquared * sigma[l];
            i++;
        }
    }
    //std::cout << "d_dW done: " << "index: " << i << "\n";
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
        sigma[j] = 1.0 / (1.0 + exp(-m_b[j] - m_oneOverSigmaSquared * xWSum));
        xWSum = 0.0;
    }

    std::vector<double> driftVector(m_NrDimensions,0);
    int vector_index = 0;
    for (int i = indx;i< indx+m_NrDimensions;i++){
        for (int j = 0; j < m_NrHiddenNodes; j++) {
            driftVector[vector_index] += sigma[j] * m_W[i][j];
        }
        driftVector[vector_index] -= m_X[i] - m_a[i];
        driftVector[vector_index] *= 2.0 * m_oneOverSigmaSquared;
        vector_index++;

    }

    return driftVector;
}

std::vector<double> qnet::getParticlePosition(int ParticleIndx){
    int startingIndx = ParticleIndx * m_NrDimensions;
    int idx;
    std::vector<double> pos(m_NrDimensions);
    for (int i = 0; i < m_NrDimensions; i++) {
        idx = startingIndx + i;
        pos[i] = m_X[idx];
    }

    return pos;
}


