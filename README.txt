dQCuts-pipeline is a full package for the analysis of proteins in python.  It begins by inputting trajectories generated from a Molecular Dynamics simulator such as AMBER, and continues by computing the relvant coordinate or angular data, reducing the dimensionality of the data, generating a similarity graph, and finally clustering the datapoints, each of which correspond to a single conformation the protein forms.

Dependencies:
	I have tried my best to keep the dependencies down, but there are a few:
		- MDAnalysis	- tool to manipulate the trajectories in python
		- Numpy		- ubiquitous numerical analysis package
		- Scipy		- provides sparse matrices
		- Matplotlib	- ubiquitous graphing package
		- SciKit-Learn	- provides knn-functionality

Setup:
	As there are two git-repositories housed under this repository, I wanted to alleviate the hassle of setup.  Thus, I created a small script to make all the necessary directories and pull all the necessary git repos.  Just run:
	$	sh setup.sh		(Internet connection needed)

	After setting up the directory, it is imperative that you edit the file 'generateConfig.py'.  This file is where you tell the entire program how you want it to run.  Be sure to examine each closely as they do affect performance and accuracy.

Use:
	Pass 'main.py' the flag '-h' or '--help' to get a feel for any flags you may want to pass when running the code.  Look through the code to get an understanding of its actions; all files should be commented well.  The most complex part (hopefully) will be telling 'dQAA' or 'cQAA' where to find your trajectories.  You need to structure the filenames, and also edit 'dQAA.py' or 'cQAA.py' to fit those formats, ideally in a way where you can loop through several trajectories.  The file 'wqaa.ipynb' has a basic example on how I have used it so far.  On a standard desktop machine (Quad-Core, 4GB RAM, 2GHz), the code should complete in a little over 30 minutes, but less than 45 minutes.  It has a timing function built in to show how long everything takes (approximately).
	
