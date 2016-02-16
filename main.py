from __future__ import division
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
parser.add_argument('--setup', action='store_true', dest='setup', default=False, help='Runs setup calculations: Cum. Sum. of cov. spectrum\nand unit radius neighbor search.');
parser.add_argument('--single', action='store_true', dest='single', default=False, help='Runs jade w/ single precision. NOT recommended.')
parser.add_argument('--smart', action='store_true', dest='smart', default=False, help='Runs jade using an alternative diagonalization setup. Refer to Cardoso\'s code for more details.')
parser.add_argument('--ica', type=str, dest='icafile', default='null', help='Skips QAA and allows user to input ICA.');

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

if val.icafile == 'null':
#	cQAA	-----------------------------------------------------------------------------
	if val.coord:
		if val.verbose: timing.log('Beginning cQAA for %s with %i trajectories @ %i dimensions...' %(config['pname'], config['numOfTraj'], config['icadim']));
		icamat = c.qaa(config, val);
		if val.save: np.save('savefiles/coord_icajade%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icamat['icajade']);
		if val.save: np.save('savefiles/coord_icacoffs%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icamat['icacoffs']);
		icacoffs = icamat['icacoffs'];
		if val.verbose: timing.log('cQAA complete...\n\nBeginning search for %i nearest neighbors...' %(config['n_neighbors']));
	else:
#	dQAA	-----------------------------------------------------------------------------
		if val.verbose: timing.log('Beginning dQAA for %s with %i trajectories @ %i dimensions...' %(config['pname'], config['numOfTraj'], config['icadim']));
		icamat = d.qaa(config, val);
		if val.save: np.save('savefiles/dih_icajade%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icamat['icajade']);
		if val.save: np.save('savefiles/dih_icacoffs%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icamat['icacoffs']);
		icacoffs = icamat['icacoffs'];
		if val.verbose: timing.log('dQAA complete...\n\nBeginning tree construction...');

else: 
	icacoffs = np.load(val.icafile);
	timing.log('Skipping QAA, using \'./%s\' as ICA.' %(val.icafile));

numSamples = icacoffs.shape[1];
if val.debug: print numSamples;
numSamples = numSamples - numSamples % 10;
if val.debug: print "new numSamples: ", numSamples;

#	KNN		-----------------------------------------------------------------------------

#	This step determines how many dimensions we should consider when computing the k-neighbors	
if val.setup:
	pts = [];
	timing.log('Setup to determine adequate affinity subspace...');
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

subset = icacoffs[ :config['affdim'], :numSamples ];
nbrs = nb(n_neighbors=config['n_neighbors'], algorithm='kd_tree').fit(subset[ :, ::10 ].T);

if val.verbose: timing.log('Tree constructed...\n\nBeginning search for %i nearest neighbors using only %i dimensions...' %(config['n_neighbors'], config['affdim']));

ind10 = np.zeros((numSamples, config['n_neighbors']*10));
junk, ind10[::10,::10] = nbrs.kneighbors(subset[ :, ::10 ].T);
ind10 = (ind10*10).astype(int);
ind10[:,0] = np.arange(numSamples);
add = np.tile(np.arange(10), numSamples).reshape((numSamples, -1), order = 'C');

#for i in range(numSamples / 10):
#	for j in range(1,config['n_neighbors']):
#		ind10[ 10*i:10*(i+1), 10*j] = np.arange(ind10[ 10*i, 10*j], ind10[ 10*i, 10*j] + 10);
"""
sub_ind *= 10;
above = sub_ind.flatten() > numSamples-10
argsort = np.argsort(above);
for i in argsort:
	if i < 1:
		num_above = i;
		break;
above_ind = map( lambda x,y: (x,y), (argsort[:num_above]//config['n_neighbors'])*10, (argsort[:num_above]%config['n_neighbors'])*10);
for i in above_ind:
	ind10[i] -= 9;"""


for i in range(9):
	ind10[(i+1)::10, (config['n_neighbors']-1)::config['n_neighbors']] = ind10[::10,(config['n_neighbors']-1)::config['n_neighbors']];

for i in range(config['n_neighbors']):
	ind10[:,10*i:10*(i+1)] = np.array( map( lambda x,y: x+y, np.tile(ind10[:,10*i], 10).reshape((numSamples, -1), order='F'), add) );

ind10[numSamples-10:,:10] -= 2*add[numSamples-10:];

""""data_vectors = subset[:, ind10];
tiled_subset = np.tile(subset, data_vectors.shape[2]).reshape((config['affdim'], numSamples, -1), order = 'F');
diff = data_vectors - tiled_subset;
dist = np.sum(diff ** 2.0, axis = 0);
data = np.exp( - dist.flatten() );"""

#	Basic loop implementation as vectorized version takes too much memory now
data = np.zeros((numSamples,config['n_neighbors']*10))
for i in range(numSamples):
	for j in range(1,config['n_neighbors']*10):
		data[i,j] = np.exp(-np.sum((subset[:,ind10[i,j]]-subset[:,i])**2.0));
	if i%(numSamples//100) == 0: timing.log( 'Affinity: %i%% complete!' %(i//(numSamples//100)) );
		

if val.debug and val.save:
	np.save('savefiles/indices_%s_%in.npy' %(config['pname'], config['n_neighbors']), ind10);

if val.verbose: timing.log('Search complete...\n\nBeginning affinity matrix generation...');

#	Affgen	-----------------------------------------------------------------------------

indptr = np.arange(0,10*config['n_neighbors']*numSamples+1, config['n_neighbors']*10);
indices = ind10.flatten();
#data = np.exp( -(dist.flatten()) );
data = data.flatten();
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

