# Restricted Boltzmann Machines for trapped electron systems

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
See the Selected_Runs folder for an example.
