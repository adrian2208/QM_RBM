import os
import glob
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import analysis

def plot_Energy_v_sigma(filename):
    # Data importing and conditioning
    data_array = analysis.csvToArray(filename)
    energies = data_array[1:]
    sigmas = data_array[0]
    #Blocking analysis for errorbars
    exp_vals, std_devs = analysis.runBlocking(energies)
    #Plotting
    analysis.createPlot(sigmas,exp_vals,std_devs,filename,r"$\sigma$", r"$\left\langle E \right\rangle$" )

def plot_Energy_v_LearningStep(filename):
    # Data importing and conditioning
    data_array = analysis.csvToArray(filename)

    energies = data_array
    learningSteps = np.arange(len(data_array))
    #Blocking analysis for errorbars
    exp_vals, std_devs = analysis.runBlocking(energies)
    #Plotting
    analysis.createPlot(learningSteps,exp_vals,std_devs,filename,"Learning Step", r"$\left\langle E \right\rangle$" )

def GetBlockingResults(filename):
    # Data importing and conditioning
    data_array = analysis.csvToArray(filename)
    energies = data_array
    #Blocking analysis for errorbars
    exp_vals, std_devs = analysis.runBlocking(energies)
    return exp_vals, std_devs

def plot_Radial_Distribution(filename):
    data_array = analysis.csvToArray(filename)
    NrSampled = data_array[0]
    boxes= np.arange(len(NrSampled))
    plt.style.use("ggplot")
    plt.figure()
    plt.xlim(0,80)
    plt.scatter(boxes,NrSampled/(sum(NrSampled)),label = "doublet")
    
    plt.ylabel("Normalized Distribution")
    plt.xlabel(r"$1/600(a_0)$")
    
    plt.legend()
    figPath = ".\\Results\\{}.pdf".format(filename)
    plt.savefig(figPath)
    plt.close()

def plot_Parameter_sweep():
    # Data importing and conditioning

    fileInfo = analysis.GetParameterSweepFileNames_and_params()
    for info in fileInfo:
        filename = info[0]
        data_array = analysis.csvToArray(filename)

        energies = data_array
        learningSteps = np.arange(len(data_array))

        #Blocking analysis for errorbars
        exp_vals, std_devs = analysis.runBlocking(energies)
        #Plotting
        analysis.createPlot(learningSteps,exp_vals,std_devs,filename,"Learning Step", r"$\left\langle E \right\rangle$" )

def plot_seaborn(filename_expVal,filename_StdDev):
    df = pd.read_csv("Output\{}.csv".format(filename_expVal), usecols = ['h','0.01','0.013103','0.016207','0.01931','0.022414','0.028621','0.031724','0.037931','0.041034','0.047241','0.050345','0.056552','0.059655','0.065862'],index_col ='h')
    mask = df < 0

    sns.heatmap(df,vmax = 3.5,vmin = 2.4,annot = True,mask = mask,fmt = '.2f', square = True,annot_kws={"fontsize":8})
    figPath = ".\\Results\\{}.pdf".format(filename_expVal)
    plt.savefig(figPath)
    plt.close()
    df = pd.read_csv("Output\{}.csv".format(filename_StdDev), usecols = ['h','0.01','0.013103','0.016207','0.01931','0.022414','0.028621','0.031724','0.037931','0.041034','0.047241','0.050345','0.056552','0.059655','0.065862'],index_col ='h')
    mask = df == 0

    sns.heatmap(df,vmax = 1.0,vmin = 0.0,annot = True,mask = mask,fmt = '.2f', square = True,annot_kws={"fontsize":8})
    figPath = ".\\Results\\{}.pdf".format(filename_StdDev)
    plt.savefig(figPath,bbox_inches="tight")

def plot_DataFrame():
    exp_val, std_dev = analysis.process_seaborn_plot_data()
    mask = exp_val < 0
    sns.heatmap(exp_val,vmax = 2.00001,vmin = 1.99992,annot = True,mask = mask,fmt = '.2f', square = True,annot_kws={"fontsize":8})
    figPath = ".\\Results\\{}.pdf".format("exp_val")
    plt.savefig(figPath)
    plt.close()
    mask = std_dev == 0
    sns.heatmap(std_dev,vmax = 0.0001,vmin = 0.0000,annot = True,mask = mask,fmt = '.2f', square = True,annot_kws={"fontsize":8})
    figPath = ".\\Results\\{}.pdf".format("std_dev")
    plt.savefig(figPath,bbox_inches="tight")


def Plot_2D_particlePositions(filename):
    matrix = analysis.getMatrixFromFile(".\\Output\\{}.csv".format(filename))
    plt.style.use("ggplot")

    sns.heatmap(matrix,vmax = 450,vmin = 10)
    figPath = ".\\Results\\{}.pdf".format(filename)
    plt.savefig(figPath,bbox_inches="tight")
    plt.close()

def plot_Selected_Runs():
    plot_Radial_Distribution("Particle_Radii")
    Plot_2D_particlePositions("D2_P_2I_Y_Importance_S_2pow21_eqS_2pow21_GD_ls_v_E_LR_0.100000_NH_3")
    plot_Energy_v_LearningStep("D2_P_2I_Y__S_2pow21_eqS_2pow21_GD_ls_v_E_LR_0.100000_NH_3_Adaptive_1")


plot_Selected_Runs()

