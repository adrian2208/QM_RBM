#pragma once
#include <string>
#include "simplegaussian.h"


class SimpleGaussian_Numerical : public SimpleGaussian {
public:
	SimpleGaussian_Numerical(class System* system, double alpha) : SimpleGaussian(system, alpha) {}
	double computeDoubleDerivative(std::vector<class Particle*> particles) override;
	std::string getName() { return m_name; };
private:
	std::string m_name = "SimpleGaussianNumerical";
};