import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from MDAnalysis import *
import argparse
import timing
import time
from wqaa.KabschAlign import *

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
		print 'Cluster-Size: %i' %(cluster_size[i]);
		avg_time = .9777642*np.median(clusters[i]);
		beg_time = .9777642*np.mean(np.sort(clusters[i])[0]);
		end_time = .9777642*np.mean(np.sort(clusters[i])[-1]);

		print 'Cluster %i median Time: %.2fns' %(i+1, avg_time);
		print 'Cluster %i Time Range: %.2fns - %.2fns' %(i+1, beg_time, end_time);
		

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

def getRMSD(config, val):
	"""
	Pull in mean structure, align coords from cluster indices (need coords then) with Kabsch
	"""
	kalign = KabschAlign();
	meancoords = np.load('savefiles/%s_meancoords_%in_%ic.npy' %(config['pname'], config['n_neighbors'], config['numClusters']));
	clustind = np.load('savefiles/%s_clustind_%in_%ic.npy' %(config['pname'], config['n_neighbors'], config['numClusters']));
	coords = np.load('savefiles/coords_%s.npy' %(config['pname']));
	
	print meancoords.shape
	print coords.shape

	ermsd = [];
	for i in range(len(clustind)):
		temp = [];
		for j in clustind[i]:
			a, b, c, d = kalign.kabsch(meancoords[i].reshape((3,-1)),coords[:,j].reshape((3,-1)));
			temp.append(c);
		ermsd.append(temp);

	if val.save: np.save('savefiles/clusterRMSD_%s.npy' %(config['pname']), ermsd);	

	
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', action='store_true', dest='graph', default=False, help='Shows graphs.')
	parser.add_argument('-v', action='store_true', dest='verbose', default=False, help='Runs program verbosely.')
	parser.add_argument('-s', '--save', action='store_true', dest='save', default=False, help='Saves clusters.')
	parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Prints debugging help.')
	parser.add_argument('-c', '--clusters', action='store_true', dest='clust', default=False, help='Gets cluster indices.');
	parser.add_argument('-m', '--mean', action='store_true', dest='mean', default=False, help='Gets mean cluster structures.');
	parser.add_argument('-r', '--rmsd', action='store_true', dest='rmsd', default=False, help='Computes cluster RMSD values.');
	values = parser.parse_args();
	global val;
	val = values;

	config = np.load('config.npz')['config'];
	config = config.reshape(1)[0];

	if val.clust: getClusters(config, val);
	if val.mean: getMean(config, val);
	if val.rmsd: getRMSD(config, val);
	
