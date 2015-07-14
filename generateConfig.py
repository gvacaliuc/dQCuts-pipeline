import numpy as np

def genConf_dqaa():
	config = {};

	#	Trajectory to start on -- 0 to start on first
	startTraj = 1;
	config['startTraj'] = startTraj;

	#	Number of trajectories to include -- EDIT
	numberOfTraj = 1;
	config['numOfTraj'] = numberOfTraj;

	#	Residue to start on (0 for second -- read ipynb on wqaa on why we exclude 1st res) -- EDIT
	startRes = 0;
	config['startRes'] = startRes;

	#	Residue to iterate through (includes ith residue if you give positive number, -2 for all but last -- again, read wqaa.ipynb) -- EDIT
	endRes = -2;
	config['endRes'] = endRes;

	#	Number of dimensions JADE should consider -- (base off point in plot of cumulative variance with greatest curvature)
	icadim = 40;
	config['icadim'] = icadim;

	#	Number of neighbors to query for when assembling affinity matrix
	n_neighbors = 31;
	config['n_neighbors'] = n_neighbors;

	#	Name of protein -- appended to beginning of all save files -- EDIT
	pname = 'kbh';
	config['pname'] = pname;

	#	Number of clusters dncuts should find
	numClusters = 16;
	config['numClusters'] = numClusters;

	affdim = 10;
	config['affdim'] = affdim;
	
	np.savez('config.npz', config=config);

def genConf_cqaa():
	config = {};

	#	Trajectory to start on -- 0 to start on first
	startTraj = 0;
	config['startTraj'] = startTraj;

	#	Number of trajectories to include -- EDIT
	numberOfTraj = 1;
	config['numOfTraj'] = numberOfTraj;

	#	Residue to start on (0 for first) -- EDIT
	startRes = 0;
	config['startRes'] = startRes;

	#	Residue to iterate through (includes ith residue if you give positive number, -1 for all) -- EDIT
	endRes = -1;
	config['endRes'] = endRes;

	#	Number of dimensions JADE should consider -- (base off point in plot of cumulative variance with greatest curvature)
	icadim = 40;
	config['icadim'] = icadim;

	#	Number of neighbors to query for when assembling affinity matrix
	n_neighbors = 31;
	config['n_neighbors'] = n_neighbors;

	#	Name of protein -- appended to beginning of all save files -- EDIT
	pname = 'kbh';
	config['pname'] = pname;

	#	Number of clusters dncuts should find
	numClusters = 16;
	config['numClusters'] = numClusters;

	affdim = 10;
	config['affdim'] = affdim;

	np.savez('config.npz', config=config);

if __name__ == '__main__':
	generateConfig();
	