import numpy as np
from numpy import *
import argparse
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
parser.add_argument('--num', type=int, dest='num_neigh', default=0, help='Number of neighbors to find.')

#parser.add_argument('--setup', action='store_true', dest='setup', default=False, help='Runs setup calculations: Cum. Sum. of cov. spectrum\nand unit radius neighbor search.')

values = parser.parse_args()
global val;
val = values;

if val.debug: val.verbose = True;

def save_sparse_csc(filename,array):
	np.savez(filename, data = array.data, indices=array.indices, indptr=array.indptr, shape=array.shape )


config = np.load('config.npz')['config'];
config = config.reshape(1)[0];
if val.num_neigh != 0: config['n_neighbors'] = val.num_neigh

icacoffs = np.load('savefiles/coord_icacoffs%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj'] + config['numOfTraj']))
numSamples = icacoffs.shape[1]
nbrs = nb(n_neighbors=config['n_neighbors'], algorithm='kd_tree').fit(icacoffs[ :config['affdim'] ].T);

if val.verbose: timing.log('Tree constructed...\n\nBeginning search for %i nearest neighbors using only %i dimensions...' %(config['n_neighbors'], config['affdim']));

distances, indices = nbrs.kneighbors(icacoffs[ :config['affdim'] ].T);


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
	np.save('savefiles/%s_eigv_%ic_%i.npy' %(config['pname'], config['numClusters'], config['n_neighbors']), eig.ev);
	np.save('savefiles/%s_eigval_%ic_%i.npy' %(config['pname'], config['numClusters'], config['n_neighbors']), eig.evl);
