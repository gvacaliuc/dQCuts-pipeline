import numpy as np

def genConf_dqaa():
	config = {};

	#	Trajectory to start on -- 0 to start on first
	startTraj = 1;
	config['startTraj'] = startTraj;

	#	Number of trajectories to include -- EDIT
	numberOfTraj = 25;
	config['numOfTraj'] = numberOfTraj;

	#	Residue to start on (0 for second -- read ipynb on wqaa on why we exclude 1st res) -- EDIT
	startRes = 0;
	config['startRes'] = startRes;

	#	Residue to iterate through -- EDIT
	numRes = 5;
	config['numRes'] = numRes;

	#	Number of dimensions JADE should consider -- (base off point in plot of cumulative variance with greatest curvature)
	icadim = 40;
	config['icadim'] = icadim;

	#	Number of neighbors to query for when assembling affinity matrix
	n_neighbors = 31;
	config['n_neighbors'] = n_neighbors;

	#	Name of protein -- appended to beginning of all save files -- EDIT
	pname = 'penta_dihedral';
	config['pname'] = pname;

	#	Number of clusters dncuts should find
	numClusters = 16;
	config['numClusters'] = numClusters;

	affdim = 10;
	config['affdim'] = affdim;
	
	#	Only uses every x number of coordinates -- useful for HUGE trajectories.
	slice_val = 10;
	config['slice_val'] = slice_val;

	kurtosis_window = 100;
	config['kurtosis_window'] = kurtosis_window;

	numBlock = 10;
	config['numBlock'] = numBlock;

	spmult_blocksize = 5;
	config['spmult_blocksize'] = spmult_blocksize;

	np.savez('config.npz', config=config);

def genConf_cqaa():
	config = {};

	#	Trajectory to start on -- 0 to start on first
	startTraj = 0;
	config['startTraj'] = startTraj;

	#	Number of trajectories to include -- EDIT
	numberOfTraj = 1;
	config['numOfTraj'] = numberOfTraj;

	#	Residue to start on (inclusive for cQAA, exclusive for dQAA, first is 0) -- EDIT
	startRes = 0;
	config['startRes'] = startRes;

	#	Number of Residues to include (inclusive for cQAA, exclusive for dQAA) -- EDIT
	numRes = 69;
	config['numRes'] = numRes;

	#	Number of dimensions JADE should consider -- (base off point in plot of cumulative variance with greatest curvature)
	icadim = 60;
	config['icadim'] = icadim;

	#	Number of neighbors to query for when assembling affinity matrix
	n_neighbors = 31;
	config['n_neighbors'] = n_neighbors;

	#	Name of protein -- appended to beginning of all save files -- EDIT
	pname = 'UBQ';
	config['pname'] = pname;

	#	Number of clusters dncuts should find
	numClusters = 16;
	config['numClusters'] = numClusters;

	affdim = 10;
	config['affdim'] = affdim;

	#	Only uses every x number of coordinates -- useful for HUGE trajectories.
	slice_val = 10;
	config['slice_val'] = slice_val;

	kurtosis_window = 100;
	config['kurtosis_window'] = kurtosis_window;

	numBlock = 10;
	config['numBlock'] = numBlock;

	spmult_blocksize = 5;
	config['spmult_blocksize'] = spmult_blocksize;

	np.savez('config.npz', config=config);

def generateConfig(i):
	if i == 0: genConf_cqaa();
	else: genConf_dqaa();

if __name__ == '__main__':
	generateConfig(0);
	
