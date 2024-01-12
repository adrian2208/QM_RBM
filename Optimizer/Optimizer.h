#pragma once
#include "wavefunction.h"
class Optimizer {
public:
	Optimizer(int NrParameters);
	virtual void optimize(std::vector<double>& a, std::vector<double>& b, double**& W, std::vector<double>& ParameterGradientVector) = 0;
protected:
	int m_NrParameters = 0;
};