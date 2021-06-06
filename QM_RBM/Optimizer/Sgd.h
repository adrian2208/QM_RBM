#pragma once
#include "Optimizer.h"

class Sgd :public Optimizer {
public:
	Sgd(int NrParameters, double LR);
	void optimize(std::vector<double>& a, std::vector<double>& b, double**& W, std::vector<double>& ParameterGradientVector) override;
private:
	double m_LR = 0;
	int m_OptimizationCycle = 0;
};