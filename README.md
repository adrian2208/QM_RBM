# Restricted Boltzmann Machines for trapped electron systems



image

![Screenshot](./Selected_Runs/Particle_Radii.pdf)

image

This program utilizes RBMs and variational Monte- Carlo methods to approximate the ground states of trapped electrons.

The Paper can be found in the Report Folder as main.pdf.
## Program structure

The main program body is implemented in c++ with separate classes for the physical system

```c++
class Hamiltonian; 
class WaveFunction;
class System;
class Sampler;
```
with various support classes
```c++
class Optimizer;// - Gradient descent of network parameters
class Timer;//"Benchmark.h" - for timing 
class csvHandler;//"csvHandler.h" - handles output files
class Random;// - generates all randomly distributed values
```
The main scope is contained in QM_RBM.cpp. All files are output to the Ouput folder, then fetched and processed using 
DataHandler.py with the library analysis.py.

## Usage
### See the Selected_Runs folder for an example.
In QM_RBM.cpp the basic structure of a simulation looks like this
```c++
int main(int argc, char const* argv[]) {
    //Setup
    class System* system = new System();
    system->setHamiltonian(new Hamiltonian());
    
    class Wavefunction* wavefunction = new Wavefunction();
    wavefunction->setOptimizer(new Optimizer());
    
    system->setWaveFunction(wavefunction);
    
    system->setSampler(new Sampler());
    system->setEquilibrationFraction();
    //Initialize datastruct for storing data for blocking - It needs to have this form
    std::vector<std::pair<std::string, std::vector<double>>> datastruct = std::vector<std::pair<std::string, std::vector<double>>>();
    
    //Run
    system->runOptimizationSteps(datastruct); //for gradient descent optimization or
    system->runMetropolisSteps(); //for brute force sampling or
    system->runMALASteps() // for Importance Sampling
    //Output handling
    csvHandler outputFile(fileName,csvComment);
    outputFile.WriteToCSV(datastruct);
    
    return 0;
}
```

