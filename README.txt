dQCuts-pipeline is a package to analyze MD simulations.  It begins by inputting trajectories generated from a Molecular Dynamics simulator such as AMBER, and continues by computing the relvant coordinate or angular data, reducing the dimensionality of the data, generating a similarity graph, and finally clustering the datapoints, each of which correspond to a single conformation the protein forms.

Dependencies:
		- MDAnalysis
		- Numpy
		- Scipy
		- Matplotlib
		- SciKit-Learn
		- Misc. Python packages ( argparse, collections, etc... )

Setup:
	$ sh setup.sh  ( makes some directories and downloads other git directories as necessary )

	After running setup, edit the file 'generateConfig.py'.  It allows you to configure the code, but it does have default values coded in.  

Use:
	Pass 'main.py' the flag '-h' or '--help' to get a feel for any flags you may want to pass when running the code.  Look through the code to get an understanding of its actions; all files should be commented well.  To run any data other than the default data, edit 'cQAA.py' or 'dQAA.py'.

Note: The code will run fine on a Linux PC, but on a mac the clustering code will run into issues.  I have a precompiled shared object running some C-code, so to run it on a mac you'll need to recompile that source.  The source for the shared object is located at www.github.com/gvacaliuc/dncuts_eigensolver/ under packages/c_extensions.  (I believe you should simply be able to run setup_c.py to create a new so, but if you run into errors just send me an email at gabe.vacaliuc@gmail.com and I'll get a dynamic library compiled and post a tutorial for OSX.
