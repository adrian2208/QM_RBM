#include "Sgd.h"

Sgd::Sgd(int NrParameters, double LR) :
    Optimizer(NrParameters) {
    m_LR = LR;
}

void Sgd::optimize(std::vector<double>& a, std::vector<double>& b, double**& W, std::vector<double>& ParameterGradientVector) {
    int NrVisibleNodes = a.size();
    int NrHiddenNodes = b.size();

    int idx = 0;
    for (int i = 0; i < NrVisibleNodes; i++) {
        a[i] = a[i] - m_LR * ParameterGradientVector[idx];
        idx++;
    }
    for (int i = 0; i < NrHiddenNodes; i++) {
        b[i] = b[i] - m_LR * ParameterGradientVector[idx];
        idx++;
    }

    for (int i = 0; i < NrVisibleNodes; i++) {
        for (int j = 0; j < NrHiddenNodes; j++) {
            W[i][j] = W[i][j] - m_LR * ParameterGradientVector[idx];
            idx++;
        }
    }

    m_OptimizationCycle++;
}
