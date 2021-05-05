#pragma once
#include "sampler.h"

class sampler_gd : public Sampler {
public:
	sampler_gd(class System* system);
	void sample(bool acceptedStep) override;
	void computeAverages() override;
};