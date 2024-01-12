import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os.path
from itertools import chain
from collections import Counter

def block(x):
	from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt
	from numpy.linalg import inv
	"""
	This code is taken verbatim from Marius Jonsson's 
	blocking code
	"""
	# preliminaries

	n = len(x)
	d = int(log2(n))
	s, gamma = zeros(d), zeros(d)
	mu = mean(x)
	# estimate the auto-covariance and variances
	# for each blocking transformation
	for i in arange(0,d):
		n = len(x)
		# estimate autocovariance of x
		gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
		# estimate variance of x
		s[i] = var(x)
		# perform blocking transformation
		x = 0.5*(x[0::2] + x[1::2])

	# generate the test observator M_k from the theorem
	M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

	# we need a list of magic numbers
	q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

	# use magic to determine when we should have stopped blocking
	for k in arange(0,d):
		if(M[k] < q[k]):
			break
	if (k >= d-1):
		print("Warning: Use more data")
	return mu, s[k]/2**(d-k)

def runBlocking(array):
	out = np.zeros((2,len(array)))
	for i in range(len(array)):
		print(array[i])
		out[:,i] = block(array[i])
	mean = out[0,:]
	std_dev = np.sqrt(out[1,:])
	return mean,std_dev

def createPlot(x,y,error=[], filename="", x_label="", y_label="", savefig = True, showfig = False):

	plt.style.use("ggplot")

	plt.figure()
	plt.xlabel(x_label)
	plt.ylabel(y_label)

	if not (error == []) :
		plt.errorbar(x,y,error, fmt = ".", capsize = 2.0, ecolor = "grey", color = "orange")

	else:
		plt.plot(x,y)
	#plt.ylim(0,10)
	if (showfig):
		plt.show()
	if (savefig):
		figPath = ".\\Results\\{}.pdf".format(filename)
		plt.savefig(figPath)

def csvToArray(filename):
    df = pd.read_csv('Output\{}.csv'.format(filename), comment='#')
    data_array = df.to_numpy()
    output_array = []
    for i in range(len(data_array[0,:])):
        col = data_array[:,i] 
        col  = col[~np.isnan(col)]
        output_array.append(col)
    return output_array

def PlotReadMe(filename):
	comments = getCSVcomment(filename)
	file = open("Results\PlotReadMe.txt","a")
	file.write("{} : \n".format(filename))
	for comment in comments:
		file.write("{}".format(comment))
	file.write("--------------------------------------------------------------------------------\n")
	file.close()


def getCSVcomment(filename):
	from itertools import takewhile
	with open('Output\{}.csv'.format(filename), 'r') as file:
		lines = takewhile(lambda s: s.startswith('#'), file)
		comments = list(lines)
	return comments

def clearOutputDirectory():
	""" Thank you! StackOverflow users Nick Stinemates and Mark Amery """
	directory = ".\Output"
	for file in os.listdir(directory):
		path = os.path.join(directory,file)
		try:
			if os.path.isfile(path) or os.path.islink(path):
				os.unlink(path)
			elif os.path.isdir(path):
				shutil.rmtree(path)
		except exception as e:
			print('Failed to delete %s. Reason: %s' % (file_path, e))
	return 0

def GetParameterSweepFileNames_and_params():
	directory = ".\Output"
	output = []
	NH_vals = []
	LR_vals = []

	for file in os.listdir(directory):
		name = os.path.splitext(file)[0]
		if (name.split("_")[-2] == "NH" and name.split("_")[-4] == "LR" ):
			NH_val = int(name.split("_")[-1])
			LR_val = float(name.split("_")[-3])
			if NH_val not in NH_vals:
				NH_vals.append(NH_val)
			if LR_val not in LR_vals:
				LR_vals.append(LR_val)
			output.append([name,NH_val,LR_val])
	y_axis = sorted(NH_vals)
	x_axis = sorted(LR_vals)
	exp_val_matrix = np.zeros((len(NH_vals),len(LR_vals)))
	std_dev_matrix = np.zeros((len(NH_vals),len(LR_vals)))
	print(NH_vals)
	print(LR_vals)
	output.sort(key=lambda x: int(x[1]))
	cols = np.array_split(output, len(NH_vals))
	i = 0
	j = 0

	for col in cols:
		for arr in col:
			#print(csvToArray(arr[0])[-1])     arr[0]
			#print(arr)#D2_P_2I_Y__S_2pow18_eqS_2pow17_GD_ls_v_E_LR_0.060000_NH_3
			dataframe = pd.read_csv('Output\{}.csv'.format(arr[0]), comment='#')
			dataframe = dataframe.apply(pd.to_numeric, errors = 'coerce')
			energies = dataframe.to_numpy()[:,-1]
			if np.isnan(energies).any():
				exp_val_matrix[i][j] = np.nan
				std_dev_matrix[i][j] = np.nan
			else:
				exp_val, variance = block(energies)
				std_dev = np.sqrt(variance)
				exp_val_matrix[i][j] = exp_val
				std_dev_matrix[i][j] = std_dev
			print(j)
			j+=1
		i+=1
		j=0
	return exp_val_matrix, std_dev_matrix, x_axis, y_axis

#def sortParameterSweep(ParameterSweepFileNames_and_params):

def process_seaborn_plot_data():
	directory = ".\Output"
	output = []
	actual_output = []
	Hs  = []
	LRs = []
	for file in os.listdir(directory):
		name = os.path.splitext(file)[0]
		if (name.split("_")[-2] == "NH" and name.split("_")[-4] == "LR" ):
			NH_val = int(name.split("_")[-1])
			LR_val = float(name.split("_")[-3])
			output.append([name,NH_val,LR_val])
	output.sort(key=lambda x: int(x[1]))
	for array in output:
		dataframe = pd.read_csv('Output\{}.csv'.format(array[0]), comment='#')
		dataframe = dataframe.apply(pd.to_numeric, errors = 'coerce')
		energies = dataframe.iloc[:, 0].to_numpy()
		if np.isnan(energies).any():
			mean = np.nan
			Variance = np.nan
		else:
			mean, Variance = block(energies)
		actual_output.append([mean, np.sqrt(Variance),array[1],array[2]])
		if array[1] not in Hs:
			Hs.append(array[1])
		if array[2] not in LRs:
			LRs.append(array[2])
	
	exp_val_matrix = np.zeros((len(Hs),len(LRs)))
	std_dev_matrix = np.zeros((len(Hs),len(LRs)))
	
	for item in actual_output:
		row_idx = Hs.index(item[2])
		col_index= LRs.index(item[3])
		exp_val_matrix[row_idx][col_index] = item[0]
		std_dev_matrix[row_idx][col_index] = item[1]
	exp_val_dataframe = pd.DataFrame(exp_val_matrix, columns = LRs, index = Hs)
	std_dev_dataframe = pd.DataFrame(std_dev_matrix, columns = LRs, index = Hs)
	return exp_val_dataframe,std_dev_dataframe

def getMatrixFromFile(FilePath):
	file = pd.read_csv(FilePath, header = None)
	return file.values
