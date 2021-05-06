#pragma once
#include <vector>
#include <string>
#include <Math/random.h>

class System {
public:
    System();
    System(int seed);
    bool metropolisStep             ();
    bool MALAStep                   ();//Metropolis- Adjusted Langevin Algorithm
    void setDriftCoefficient        (double stepLength);
    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void runMALASteps               (int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setTimeStep                (double timeStep);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setFileOptString           (std::string fileOptString);
    void setSampler                 (class Sampler* sampler);
    void setElapsedTime             (double elapsed);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    class Random*                   getRandomEngine()   { return m_random; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getStepLength()              { return m_stepLength; }
    double getElapsedTime()             { return m_elapsed; }
private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_timeStep = 0;
    double                          m_timeStepSqrt = 0;
    double                          m_rootTimeStep = sqrt(m_stepLength);
    double                          m_driftCoefficient = 0.5;
    double                          m_elapsed = 0;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    class Random*                   m_random = nullptr;
    std::string                     m_fileOptString = "";
};

