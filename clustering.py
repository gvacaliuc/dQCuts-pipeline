import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from MDAnalysis import *
import argparse
import timing
import time

def getClusters(config, val):

	eigvec = np.load('savefiles/%s_eigv_%ic.npy' %(config['pname'], config['numClusters']));
	clusters = [];
	
	timing.log('Data file accepted, computing cluster indices...');

	#	Finds all indices outside of 2 stddevs of eigenvector mean
	mean = np.mean(eigvec, axis = 0);
	std = np.std(eigvec, axis = 0);
	inClust = np.greater(np.abs(np.subtract(eigvec, mean)),2*std);
	cluster_size = np.sum(inClust, axis = 0);
	indices = np.argsort(-inClust, axis = 0);
	
	for i in range(config['numClusters']):
		clusters.append(indices[:cluster_size[i],i]);

	if val.save: np.save('savefiles/%s_clustind_%in_%ic.npy' %(config['pname'], config['n_neighbors'], config['numClusters']), clusters);

	timing.log('Clusters saved!');

def getMean(config, val):

	coords = np.load('savefiles/%s_coords.npy' %(config['pname']));
	indices = np.load('savefiles/%s_clustind_%in_%ic.npy' %(config['pname'], config['n_neighbors'], config['numClusters']));
	meancoords = [];
	
	for i in range(config['numClusters']):
		meancoords.append( np.mean(coords[:, indices[i]], axis = 1) );

	meancoords = np.array(meancoords);	
	
	if val.graph:
		
		for i in range(config['numClusters']):
			fig = plt.figure()
			ax = fig.gca(projection='3d');
			ax.plot(meancoords[i][::3], meancoords[i][1::3], meancoords[i][2::3]);
			plt.show();
	
	if val.save: np.save('savefiles/%s_meancoords_%in_%ic.npy' %(config['pname'], config['n_neighbors'], config['numClusters']), meancoords);

	
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', action='store_true', dest='graph', default=False, help='Shows graphs.')
	parser.add_argument('-v', action='store_true', dest='verbose', default=False, help='Runs program verbosely.')
	parser.add_argument('-s', '--save', action='store_true', dest='save', default=False, help='Saves clusters.')
	parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Prints debugging help.')

	values = parser.parse_args();
	global val;
	val = values;

	config = np.load('config.npz')['config'];
	config = config.reshape(1)[0];

	getClusters(config, val);
	getMean(config, val);
	
