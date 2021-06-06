#include "Adam.h"

Adam::Adam(int NrParameters, double LR) :
    Optimizer(NrParameters) {
    m_LR = LR;
}

void Adam::optimize(std::vector<double> &a, std::vector<double>& b, double** &W, std::vector<double>& ParameterGradientVector){
    int NrVisibleNodes = a.size();
    int NrHiddenNodes = b.size();

    for (int i = 0; i < m_NrParameters; i++) {
        m_Mw[i] = m_b1 * m_Mw[i] + (1 - m_b1) * ParameterGradientVector[i];
        m_Vw[i] = m_b2 * m_Vw[i] + (1 - m_b2) * pow(ParameterGradientVector[i], 2.0);
        
        m_MHatw[i] = m_Mw[i]/(1-pow(m_b1,m_OptimizationCycle+1));
        m_VHatw[i] = m_Vw[i] / (1 - pow(m_b2, m_OptimizationCycle+1));
    }

    int idx = 0;
    for (int i = 0; i < NrVisibleNodes; i++) {
        a[i] = a[i] - m_LR * m_MHatw[idx]/(sqrt(m_VHatw[idx])+1.0e-8);
        idx++;
    }
    for (int i = 0; i < NrHiddenNodes; i++) {
        b[i] = b[i] - m_LR * m_MHatw[idx] / (sqrt(m_VHatw[idx]) + 1.0e-8);
        idx++;
    }

    for (int i = 0; i < NrVisibleNodes; i++) {
        for (int j = 0; j < NrHiddenNodes; j++) {
            W[i][j] = W[i][j] - m_LR * m_MHatw[idx] / (sqrt(m_VHatw[idx]) + 1.0e-8);
            idx++;
        }
    }

    m_OptimizationCycle++;
}

