This folder contains selected results from the program   ##
and how to reproduce them								 ##
###########################################################

figure : D2_P_2I_Y__S_2pow19_eqS_2pow18_GD_ls_v_E_LR_0.050000_NH_2.pdf
description: the Local energy of two interacting particles in two dimensions plotted over the learning step in SGD with learning rate 0.05 and two hidden units
	In QM_RBM.cpp : 
		- Insert GradientDescent(0.05, 2, true, 100, true, false,2,2,2020); into the main{} scope
		- Comment out any other functions in the main{} scope
		- Delete any potential files in the Output directory
		- compile from CMakeLists.txt and run QM_RBM.exe
	In DataHandler.py :
		- At the mottom of the file, comment out any existing function calls
		- insert plot_Energy_v_LearningStep("D2_P_2I_Y__S_2pow18_eqS_2pow14_GD_ls_v_E_LR_0.050000_NH_2")
		- run the python program
		- the resulting graph should be saved to the Results folder, but not in the Report_Results folder

---------------------------------------------------------------------------------------------------------------------------------------------------------------------

figure : D2_P_2I_Y_Importance_S_2pow19_eqS_2pow18_Position_sampling_P_2_NH_2_I_1.pdf
description: Position sampling of two interacting particles in two dimensions
	In QM_RBM.cpp : 
		- Insert PositionSampling(2, 2, true, true, 2020, true); into the main{} scope
		- Comment out any other functions in the main{} scope
		- Delete any potential files in the Output directory
		- compile from CMakeLists.txt and run QM_RBM.exe
	In DataHandler.py :
		- At the mottom of the file, comment out any existing function calls
		- insert Plot_2D_particlePositions("D2_P_2I_Y_Importance_S_2pow19_eqS_2pow18_Position_sampling_P_2_NH_2_I_1")
		- run the python program
		- the resulting graph should be saved to the Results folder, but not in the Report_Results folder