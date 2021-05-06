import os
import glob
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import dataModule as data
#############################################################
###################     IMPORTATNT    #######################
### COMPILE C++ PROGRAM AS X64- RELEASE #####################
### PYTHON ONLY ACESSES THE X64- RELEASE FOLDER #############
#############################################################

###############################################################################
##############			OUTPUT FILE COLUMN NAME BY INDEX	    ###############
#PARTICLES[0], #DIMENSIONS[1], #METROPOLIS STEPS[2], EQUILIBRATION FRACTION[3]#
#ALPHA[4], STEP LENGTH[5], ENERGY[6]	#######################################

##############          Alpha dependence - local energy       #################

def plot_Energy_v_alpha_for_DifferentNrParticles(numberOfDimensions, numberOfParticles_list, alphas):
    #numberOfDimensions  = 1
    #numberOfParticles_list  = [1,2,3,4,5]
    numberOfSteps       = 2**16
    omega				= 1
    #alphas				= np.linspace(0.2,0.8,4)
    stepLength			= 0.1
    x_list = []
    y_list = []
    error_list = []
    legend_list = []

    for numberOfParticles in numberOfParticles_list:
        x,y,errors, NrParticles, NrDims, NrSteps, wf,Hamiltonian, elapsedTime = data.RunVMC_w_blocking(numberOfDimensions,numberOfParticles,numberOfSteps,alphas,4,6)
        x_list.append(x)
        y_list.append(y)
        error_list.append(errors)
        legend_list.append("#Particles = " + str(NrParticles))
    
    data.createFig(x_list,y_list,error_list, Hamiltonian+" "+wf, "alpha", "Mean Energy", legend_list, errorbars = True, )
    
    return 0

def plot_PnrNormalizedEnergy_v_alpha_for_DifferentNrParticles(numberOfDimensions, numberOfParticles_list, alphas):
    #numberOfDimensions  = 1
    #numberOfParticles_list  = [1,2,3,4,5]
    exponent = 20
    numberOfSteps       = 2**exponent
    omega				= 1
    #alphas				= np.linspace(0.2,0.8,4)
    stepLength			= 0.1
    x_list = []
    y_list = []
    error_list = []
    legend_list = []

    for numberOfParticles in numberOfParticles_list:
        x,y,errors, NrParticles, NrDims, NrSteps, wf,Hamiltonian, elapsedTime = data.RunVMC_w_blocking(numberOfDimensions,numberOfParticles,numberOfSteps,alphas,4,6)
        x_list.append(x)
        y_list.append(y/NrParticles)
        error_list.append(errors/NrParticles)
        legend_list.append("#Particles = " + str(NrParticles))
    
    data.createFig(x_list,y_list,error_list,"{}_{}_D{}_P{}_2pow{}_ImpSamp".format(Hamiltonian, wf,numberOfDimensions,numberOfParticles_list,exponent), r"$\alpha$", r"$\frac{\langle E\rangle}{N}$", legend_list, errorbars = True)
    
    return 0

def task_b_plot_analytical(numberOfDimensions, numberOfParticles_list, alphas):
    #numberOfDimensions  = 1
    #numberOfParticles_list  = [1,2,3,4,5]
    exponent = 22
    numberOfSteps       = 2**exponent
    omega				= 1
    #alphas				= np.linspace(0.2,0.8,4)
    stepLength			= 0.1
    x_list = []
    y_list = []
    error_list = []
    legend_list = []
    elapsed_time_list = []

    for numberOfParticles in numberOfParticles_list:
        x,y,errors, NrParticles, NrDims, NrSteps, wf,Hamiltonian, elapsedTime = data.RunVMC_w_blocking(numberOfDimensions,numberOfParticles,numberOfSteps,alphas,4,6)
        x_list.append(x)
        y_list.append(y/NrParticles)
        error_list.append(errors/NrParticles)
        legend_list.append("#P = {} ({:.2f} ms)".format(NrParticles,elapsedTime))
    data.createFig(x_list,y_list,error_list, "{}_{}_D{}_P{}_2pow{}_analytic".format(Hamiltonian, wf,numberOfDimensions,numberOfParticles_list,exponent), r"$\alpha$", r"$\frac{\langle E\rangle}{N}$", legend_list, errorbars = True)
    
    return 0

def task_b_plot_numerical(numberOfDimensions, numberOfParticles_list, alphas):
    #numberOfDimensions  = 1
    #numberOfParticles_list  = [1,2,3,4,5]
    exponent = 22
    numberOfSteps       = 2**exponent
    omega				= 1
    #alphas				= np.linspace(0.2,0.8,4)
    stepLength			= 0.1
    x_list = []
    y_list = []
    error_list = []
    legend_list = []
    elapsed_time_list = []

    for numberOfParticles in numberOfParticles_list:
        x,y,errors, NrParticles, NrDims, NrSteps, wf,Hamiltonian, elapsedTime = data.RunVMC_w_blocking(numberOfDimensions,numberOfParticles,numberOfSteps,alphas,4,6)
        x_list.append(x)
        y_list.append(y/NrParticles)
        error_list.append(errors/NrParticles)
        legend_list.append("#P = {} ({:.2f} ms)".format(NrParticles,elapsedTime))
    data.createFig(x_list,y_list,error_list, "{}_{}_D{}_P{}_2pow{}_numerical".format(Hamiltonian, wf,numberOfDimensions,numberOfParticles_list,exponent), r"$\alpha$", r"$\frac{\langle E\rangle}{N}$", legend_list, errorbars = True)
    
    return 0

def task_c_plot(numberOfDimensions, numberOfParticles_list, alphas, timeStep):
    #numberOfDimensions  = 1
    #numberOfParticles_list  = [1,2,3,4,5]
    exponent = 22
    numberOfSteps       = 2**exponent
    omega				= 1
    #alphas				= np.linspace(0.2,0.8,4)
    stepLength			= 0.1
    x_list = []
    y_list = []
    error_list = []
    legend_list = []
    elapsed_time_list = []
    timeStep_list = []

    for numberOfParticles in numberOfParticles_list:
        x,y,errors, NrParticles, NrDims, NrSteps, wf,Hamiltonian, elapsedTime = data.RunVMC_w_blocking(numberOfDimensions,numberOfParticles,numberOfSteps,alphas,4,6)
        x_list.append(x)
        y_list.append(y/NrParticles)
        error_list.append(errors/NrParticles)
        legend_list.append("#P = {} ({:.2f} ms)".format(NrParticles,elapsedTime))
    data.createFig(x_list,y_list,error_list, "{}_{}_D{}_P{}_2pow{}_task_C_1".format(Hamiltonian, wf,numberOfDimensions,numberOfParticles_list,exponent), r"$\alpha$", r"$\frac{\langle E\rangle}{N}$", legend_list, errorbars = True)
    
    return 0

def task_c_plot_timeStep(numberOfDimensions, numberOfParticles_list, alpha, timeSteps):
    #numberOfDimensions  = 1
    #numberOfParticles_list  = [1,2,3,4,5]
    exponent = 18
    numberOfSteps       = 2**exponent
    omega				= 1
    #alphas				= np.linspace(0.2,0.8,4)
    stepLength			= 0.1
    y_list = []
    timeSteps_list = []
    legend_list = []
    elapsed_time_list = []

    for numberOfParticles in numberOfParticles_list:
        y,NrParticles, NrDims,Nrsteps, wf, Hamiltonian  = data.RunVMC_w_blocking_varyTimestep_constantAlpha(numberOfDimensions,numberOfParticles,numberOfSteps,alpha,timeSteps)
        y_list.append(y)
        timeSteps_list.append(timeSteps)
        legend_list.append("#P = {}".format(numberOfParticles))
    data.createFig(timeSteps_list,y_list, [] ,"{}_{}_D{}_P{}_2pow{}_task_C_2".format(Hamiltonian, wf,numberOfDimensions,numberOfParticles_list,exponent), "timestep","accepted Proportion" , legend_list)
    
    return 0

def Plot_2D_particlePositions():
    NrParticles = 150
    NrSteps = 2**20
    alpha = 0.5
    NrSamplingLengths = 5
    SamplingRadius = 2
    data.clearOutputDirectory()

    hardCoreDiameter = 0
    FileOptString = "HC_diameter_zero"
    data.execute([NrParticles, NrSteps, alpha, hardCoreDiameter, NrSamplingLengths, SamplingRadius, FileOptString])
    
    hardCoreDiameter = 0.00433
    FileOptString = "HC_diameter_regular"
    data.execute([NrParticles, NrSteps, alpha, hardCoreDiameter, NrSamplingLengths, SamplingRadius, FileOptString])

    #hardCoreDiameter = 0.433
    #FileOptString = "HC_diameter_bigger"
    #data.execute([NrParticles, NrSteps, alpha, hardCoreDiameter, NrSamplingLengths, SamplingRadius, FileOptString])
    
    NoJastrow_matrix = data.getMatrixFromFile(".\\Output\\CorrelatedGaussian_EllipticOscillator_HC_diameter_zero_Matrix.csv")
    Jastrow_matrix = data.getMatrixFromFile(".\\Output\\CorrelatedGaussian_EllipticOscillator_HC_diameter_regular_Matrix.csv")
    #BigJastrow_matrix = data.getMatrixFromFile(".\\Output\\CorrelatedGaussian_EllipticOscillator_HC_diameter_bigger_Matrix.csv")
    
    data.createMatrixPlot(NoJastrow_matrix)
    data.createMatrixPlot(Jastrow_matrix)
    #data.createMatrixPlot(BigJastrow_matrix)
    data.createMatrixPlot(Jastrow_matrix-NoJastrow_matrix)


#plot_PnrNormalizedEnergy_v_alpha_for_DifferentNrParticles(3, [5], np.linspace(0.5,0.7,1))
#plot_PnrNormalizedEnergy_v_alpha_for_DifferentNrParticles(2, [1,10,100,500], np.linspace(0.3,0.8,9))
#plot_PnrNormalizedEnergy_v_alpha_for_DifferentNrParticles(3, [1,10,100,500], np.linspace(0.3,0.8,9))
#task_c_plot_timeStep(3, [2], 0.5, [0.1,0.1])
#task_c_plot(3, [1,10,100,500], np.linspace(0.3,0.8,10), 0.01)

#Plot_2D_particlePositions()


