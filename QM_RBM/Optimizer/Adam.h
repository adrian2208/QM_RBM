#pragma once
#include "Optimizer.h"

class Adam :public Optimizer {
public:
	Adam(int NrParameters, double LR);
	void optimize(std::vector<double>& a, std::vector<double>& b, double**& W, std::vector<double>& ParameterGradientVector) override;
private:
	std::vector<double> m_Mw = std::vector<double>(m_NrParameters,0);
	std::vector<double> m_Vw = std::vector<double>(m_NrParameters,0);
	std::vector<double> m_MHatw = std::vector<double>(m_NrParameters, 0);
	std::vector<double> m_VHatw = std::vector<double>(m_NrParameters, 0);
	double m_b1 = 0.9;
	double m_b2 = 0.999;
	double m_LR = 0.0;
	int m_OptimizationCycle = 0;
};