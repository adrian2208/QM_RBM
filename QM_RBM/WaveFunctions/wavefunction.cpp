#include "wavefunction.h"


WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

int WaveFunction::getNrVisibleNodes()
{
    return m_NrVisibleNodes;
}
int WaveFunction::getNrHiddenNodes()
{
    return m_NrHiddenNodes;
}
int WaveFunction::getNrDimensions()
{
    return m_NrDimensions;
}

int WaveFunction::getNrParticles()
{
    return m_NrParticles;
}

int WaveFunction::getNrParameters()
{
    return m_NrParameters;
}

std::vector<double> WaveFunction::getX()
{
    return m_X;
}