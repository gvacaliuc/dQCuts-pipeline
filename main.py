import numpy as np
from numpy import *
import argparse
import wqaa.dQAA as d
import wqaa.cQAA as c
from sklearn.neighbors import NearestNeighbors as nb
from scipy.sparse import *
from dncuts_eigensolver.packages.dncuts import *
import collections
import timing
from generateConfig import *

parser = argparse.ArgumentParser()
parser.add_argument('-g', action='store_true', dest='graph', default=False, help='Shows graph of Affinity Matrix and eigenv\'s, depending on flags.')
parser.add_argument('-v', action='store_true', dest='verbose', default=False, help='Runs program verbosely.')
parser.add_argument('-s', '--save', action='store_true', dest='save', default=False, help='Saves affinity matrix and eigenv\'s.')
parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Prints debugging help.')
parser.add_argument('-c', '--coord', action='store_true', dest='coord', default=False, help='Runs cQAA instead of dQAA.')
parser.add_argument('--setup', action='store_true', dest='setup', default=False, help='Runs setup calculations: Cum. Sum. of cov. spectrum\nand unit radius neighbor search.')

values = parser.parse_args()
global val;
val = values;

if val.debug: val.verbose = True;

def save_sparse_csc(filename,array):
	np.savez(filename, data = array.data, indices=array.indices, indptr=array.indptr, shape=array.shape )

#	Runs generateConfig again, to make sure all config settings are updated
if val.coord: genConf_cqaa();
else: genConf_dqaa();

#Unpacks config settings
config = np.load('config.npz')['config'];
config = config.reshape(1)[0];

#Refuses to run if save file's prefix is not set up
if (config['pname'] == 'PNAME'):
	raise ValueError('Set \'pname\' to your protein name in \'generateConfig.py\'.');
	exit()

#	cQAA	-----------------------------------------------------------------------------
if val.coord:
	if val.verbose: timing.log('Beginning cQAA for %s with %i trajectories @ %i dimensions...' %(config['pname'], config['numOfTraj'], config['icadim']));
	icamat = c.qaa(config, val);
	if val.save: np.save('savefiles/coord_icajade%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icamat['icajade']);
	if val.save: np.save('savefiles/coord_icacoffs%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icamat['icacoffs']);
	icacoffs = icamat['icacoffs'];
	numSamples = icacoffs.shape[1];
	if val.verbose: timing.log('cQAA complete...\n\nBeginning search for %i nearest neighbors...' %(config['n_neighbors']));
else:
#	dQAA	-----------------------------------------------------------------------------
	if val.verbose: timing.log('Beginning dQAA for %s with %i trajectories @ %i dimensions...' %(config['pname'], config['numOfTraj'], config['icadim']));
	icamat = d.qaa(config, val);
	if val.save: np.save('savefiles/dih_icajade%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icamat['icajade']);
	if val.save: np.save('savefiles/dih_icacoffs%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icamat['icacoffs']);
	icacoffs = icamat['icacoffs'];
	numSamples = icacoffs.shape[1];
	if val.verbose: timing.log('dQAA complete...\n\nBeginning tree construction...');

#	KNN		-----------------------------------------------------------------------------	
if val.setup:
	pts = [];
	timing.log('setup to determine adequate affinity subspace...');
	for i in range(2,config['icadim']):
		print('Finding neighbors in %i dimensions...' %(i));
		#	Need to slice every 10 due to computation difficulties with the whole dataset
		nbrs = nb(algorithm='kd_tree').fit(icacoffs[:i,::10].T);
		ind = nbrs.radius_neighbors(icacoffs[:i,::10].T, return_distance=False);
		sz = [];
		for n in ind:
			sz.append(len(n));
		pts.append(np.mean(sz));
		if np.mean(sz) < 10:
			print 'Average number of neighbors has fallen below 10...'
			break;
	plt.plot( np.arange(2,len(pts)+2) , np.array(pts) );
	plt.yscale('log');
	plt.xscale('log');
	plt.xlabel('Number of Dimensions Considered')
	plt.ylabel('Average Number of Neighbors')
	plt.title('Average Neighbors within a Unit Radius of a Point') 
	plt.show();
	a = input('Enter desired affinity dimension (enter -1 for default): ');
	if (a > 0):
		config['affdim'] = a;

	#	Have to transpose icacoffs b/c sklearn wants Samples x Sensors
nbrs = nb(n_neighbors=config['n_neighbors'], algorithm='kd_tree').fit(icacoffs[ :config['affdim'] ].T);

if val.verbose: timing.log('Tree constructed...\n\nBeginning search for %i nearest neighbors using only %i dimensions...' %(config['n_neighbors'], config['affdim']));

distances, indices = nbrs.kneighbors(icacoffs[ :config['affdim'] ].T);

if val.debug and val.save:
	np.save('savefiles/distances_%s_%in.npy' %(config['pname'], config['n_neighbors']), distances);
	np.save('savefiles/indices_%s_%in.npy' %(config['pname'], config['n_neighbors']), indices);

if val.verbose: timing.log('Search complete...\n\nBeginning affinity matrix generation...');

#	Affgen	-----------------------------------------------------------------------------

indptr = np.arange(0,config['n_neighbors']*numSamples+1, config['n_neighbors']);
indices = indices.flatten();
data = np.exp(-(distances.flatten() ** 2.0));

affMat = csc_matrix( (data, indices, indptr), shape = (numSamples,numSamples) );

if val.graph:
	plt.spy(affMat, precision = .1, markersize = 1);
	plt.show();

if val.save:
	save_sparse_csc('savefiles/%s_aff_%id_%in.npz' %(config['pname'], config['icadim'], config['n_neighbors']), affMat);

if val.verbose: timing.log('Generation complete...\n\nBeginning clustering for %i clusters...' %(config['numClusters']));

#	DNCuts	-----------------------------------------------------------------------------

	#	eig is a collection of two matrices, eig.ev = eigenvectors, eig.evl = eigenvalues
eig = dncuts(affMat, val, nvec=config['numClusters'], n_downsample=2, decimate=2);
if val.save:
	np.save('savefiles/%s_eigv_%ic.npy' %(config['pname'], config['numClusters']), eig.ev);
	np.save('savefiles/%s_eigval_%ic.npy' %(config['pname'], config['numClusters']), eig.evl);

if val.verbose: timing.log('Clustering complete...');
#	K-means? 

