import os
import shutil
import glob
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if not os.path.exists(".\Output"): #Checks is path = "./directory" exists
	os.makedirs(".\Output")

if not os.path.exists(".\Results"): #Checks is path = "./directory" exists
	os.makedirs(".\Results")


def Get_newest_file_path():
	All_files_in_Output = glob.glob('.\\Output\\*')
	NewestFile = max(All_files_in_Output, key = os.path.getctime)
	return NewestFile


def execute(input_args):
	path = ".\\out\\build\\x64-Release\\vmc.exe"
	for arg in input_args:
		path+= " " + str(arg)
		p = subprocess.Popen(path, shell = True)
		p.wait()
	return 0

def variationalCommandLineArgs(NrDimensions, NrParticles, NrStepsPow2, alphas, timeStep = 0.01):
	n = len(alphas)
	out = []
	for i in range(n):
		out.append([NrDimensions, NrParticles, NrStepsPow2, alphas[i],timeStep,i])
	return out

def varyTimeStepCommandLineArgs(NrDimensions, NrParticles, NrStepsPow2,  timeSteps, alpha = 0.5):
	n = len(timeSteps)
	out = []
	for i in range(n):
		out.append([NrDimensions, NrParticles, NrStepsPow2, alpha,timeSteps[i],i])
	return out


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

def RunVMC_w_blocking(NrDimensions, NrParticles, NrStepsPow2, alphas,x_col, y_col, timeStep = 0.01):
	clearOutputDirectory()
	commandLineArgs = variationalCommandLineArgs(NrDimensions, NrParticles, NrStepsPow2, alphas)
	for args in commandLineArgs:
		execute(args)
	#run blocking on the data that was just printed to .csv files
	blockingResult, Main_csv_file, time_file = runBlocking()
	#Append those results of the existing dataframe that will be used for plotting
	x,y = getXYData(Main_csv_file,x_col,y_col)
	errors = [np.sqrt(e) for e in blockingResult[:,1]]
	elapsed = getTimeData(time_file)
	NrParticles, NrDims, NrSteps, wf,Hamiltonian = getRunType(Main_csv_file)
	return x,y,errors, NrParticles, NrDims, NrSteps, wf,Hamiltonian, elapsed

def RunVMC_w_blocking_varyTimestep_constantAlpha(NrDimensions, NrParticles, NrStepsPow2, alpha, timeSteps):
	clearOutputDirectory()
	commandLineArgs = varyTimeStepCommandLineArgs(NrDimensions, NrParticles, NrStepsPow2, timeSteps,alpha)
	for args in commandLineArgs:
		execute(args)

	directory = ".\Output"
	filelist = []
	for file in os.listdir(directory):
		name = os.path.splitext(file)[0]
		if name.split("_")[-1] == "accepted":
			filelist.append(file)
		elif name.split("_")[-1] == "Energies":
			energies_file = file
		elif name.split("_")[-1] == "ElapsedTime":
			elapsedTime_file = file
		elif name.split("_")[-1] == "Matrix":
			matrix_file = file
		else:
			Main_csv_file = file
	print(Main_csv_file)
	n = len(filelist)
	sortedlist = [None]*n
	for filename in filelist:
		name = os.path.splitext(filename)[0]
		sortedlist[int(filename.split("_")[-2])] = filename

	accepted = []
	for file in sortedlist:
		accepted.append(getTimeData(".\\Output\\" + file))
	NrParticles, NrDims, NrSteps, wf,Hamiltonian = getRunType(".\\Output\\"+ Main_csv_file)
	return accepted, NrParticles, NrDims, NrSteps, wf,Hamiltonian


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

def runBlocking():
	# Find all .csv files in Output directory
	# containing list of energies (identified by "Energies" in their filename)
	directory = ".\Output"
	filelist = []
	for file in os.listdir(directory):
		name = os.path.splitext(file)[0]
		if name.split("_")[-1] == "Energies":
			filelist.append(file)
		elif name.split("_")[-1] == "ElapsedTime":
			time_file = file
		elif name.split("_")[-1] == "accepted":
			accepted_fil = file
		elif name.split("_")[-1] == "Matrix":
			matrix_file = file
		else:
			Main_csv_file = file
	#Then make sure those file's filenames are sorted by their identifying number
	n = len(filelist)
	sortedlist = [None]*n
	for filename in filelist:
		name = os.path.splitext(filename)[0]
		sortedlist[int(filename.split("_")[-2])] = filename
	# Finally, we perform the blocking function on each file 
	# and append the result to the output
	output = np.zeros((n,2))
	path = ".\\Output\\"
	for i in range(n):
		y = pd.read_csv(path + sortedlist[i], header = None)[0].to_numpy()
		output[i,:] = block(y)
	return output, path + Main_csv_file, path + time_file

def getXYData(FilePath,x_col_indx, y_col_indx):
	file = pd.read_csv(FilePath, header = None)

	y =  file[y_col_indx][:].to_numpy()
	x =  file[x_col_indx][:].to_numpy()
	return x, y

def getMatrixFromFile(FilePath):
	file = pd.read_csv(FilePath, header = None)
	return file.values

def getTimeData(FilePath):
	file = pd.read_csv(FilePath, header = None)
	time =  file[0][0]
	return time

def getRunType(FilePath):
	file = pd.read_csv(FilePath, header = None)

	NrParticles =  file[0][0]
	NrDims =  file[1][0]
	NrSteps =  file[2][0]
	name = os.path.splitext(os.path.basename(FilePath))[0].split("_")
	wf = name[0]
	Hamiltonian = name[1]
	return NrParticles, NrDims, NrSteps, wf,Hamiltonian



def createFig(x,y,error, title, x_label, y_label, legend, errorbars = False, savefig = True, showfig = False, elapsedTime = []):
	""" 
	enter x, y and error each as an array of arrays containing the plotting data for each
    plot line 
    """
	plt.style.use("ggplot")

	plt.figure()
	plt.xlabel(x_label)
	plt.ylabel(y_label)

	if errorbars:
		for i in range(len(x)):
			plt.errorbar(x[i],y[i],error[i],label = legend[i])
			if elapsedTime:
				for j in range(len(elapsedTime)):
					plt.annotate("{:.2f} ms".format(elapsedTime[i]), (x[i][j],y[i][j]))
	else:
		for i in range(len(x)):
			plt.plot(x[i],y[i],label = legend[i])
			if elapsedTime:
				for j in range(len(elapsedTime)):
					plt.annotate("{:.2f} ms".format(elapsedTime[i]), (x[i][j],y[i][j]))


	plt.legend()
	if (showfig):
		plt.show()
	if (savefig):
		figPath = ".\\Results\\{}_NoTitle.pdf".format(title)
		plt.savefig(figPath)
		plt.title(title)
		figPath = ".\\Results\\{}_.pdf".format(title)
		plt.savefig(figPath)

def createMatrixPlot(matrix):
	""" 
	enter x, y and error each as an array of arrays containing the plotting data for each
    plot line 
    """
	plt.style.use("ggplot")

	plt.imshow(matrix, interpolation = "none")
	plt.colorbar()
	plt.show()
	#figPath = ".\\Results\\{}.pdf".format(title)
	#plt.savefig(figPath)


###################ignore us, we're old and stupid ###################################
"""
def Read_file_and_plot(variational_parameter,y_cols, FilePath, Xname, Yname):
	file = pd.read_csv(FilePath, header = None)

	y =  file[y_cols[0]][:]
	x = variational_parameter
	plt.plot(x,y)
	plt.xlabel(Xname)
	plt.ylabel(Yname)
	fig = plt.gcf()
	plt.show()
	User_input = input("Save results= Y/N?").lower()
	if (User_input == "y"):
		figName = ".\\Results\\{}_x_{}_y_{}.pdf".format(os.path.splitext(os.path.basename(FilePath))[0], Xname, Yname)
		plt.draw()
		fig.savefig(figName)
		os.rename(FilePath,".\\Results\\{}".format(os.path.basename(FilePath)))
	else:
		os.remove(Get_newest_file_path())

def execute_cpp_script(input_args, Variational_param_Index):
	#os.remove(Get_newest_file_path())
	var_args = input_args[Variational_param_Index]
	del input_args[Variational_param_Index]
	for var_arg in var_args:
		args  = input_args.copy()
		args.insert(Variational_param_Index, var_arg)
		path = ".\\out\\build\\x64-Release\\vmc.exe"
		for arg in args:
			path+= " " + str(arg)
		p = subprocess.Popen(path, shell = True)
		p.wait()
	return 0
"""