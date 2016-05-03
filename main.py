from __future__ import division
import numpy as np
from numpy import *
import argparse
#import wqaa.dQAA as d
#import wqaa.cQAA as c
import cQAA as c
import cQAA as d
from sklearn.neighbors import NearestNeighbors as nb
from scipy.sparse import *
from dncuts_eigensolver.packages.dncuts import *
import collections
import timing
from generateConfig import *
import sys

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

#================================================================================
#	Supplementary Code
def save_sparse_csc(filename,array):
	np.savez(filename, data = array.data, indices=array.indices, indptr=array.indptr, shape=array.shape )

def update_progress(progress):
    barLength = 40 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
#=================================================================================

#	Runs generateConfig again, to make sure all config settings are updated based off saved file
if val.coord: genConf_cqaa();
else: genConf_dqaa();

#	Unpacks config settings
config = np.load('config.npz')['config'];
config = config.reshape(1)[0];

#	Refuses to run if save file's prefix is not set up
if (config['pname'] == 'PNAME'):
	raise ValueError('Set \'pname\' to your protein name in \'generateConfig.py\'.');
	exit()

if val.icafile == 'null':
#	cQAA	-----------------------------------------------------------------------------
	if val.coord:
		if val.verbose: timing.log('Beginning cQAA for %s with %i trajectories @ %i dimensions...' %(config['pname'], config['numOfTraj'], config['icadim']));
		icajade, filename, mapshape = c.qaa(config, val);
		icacoffs = np.memmap(filename, dtype='float64', shape=mapshape);
		if val.save: np.save('savefiles/coord_icajade%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icajade);
		if val.save: np.save('savefiles/coord_icacoffs%s_%id_%i-%it.npy' %(config['pname'], config['icadim'], config['startTraj']+1, config['startTraj']+config['numOfTraj']), icacoffs[:,:]);
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
	if val.icafile[-3:] == 'npy': icacoffs = np.load(val.icafile);
	elif val.icafile[-5:] == 'array': icacoffs = np.memmap(val.icafile, dtype='float64').reshape(config['icadim'],-1);
	else: raise ValueError('Please enter valid array or memmapped array for ICA.');
	timing.log('Skipping QAA, using \'./%s\' as ICA.' %(val.icafile));

numSamples = icacoffs.shape[1];
if val.debug: print numSamples;
numSamples = numSamples - numSamples % config['numBlock'];	#	Just so it fits in nicely (Losing a few samples won't matter.)
if val.debug: print "new numSamples: ", numSamples;

#	KNN		-----------------------------------------------------------------------------

#	This step determines how many dimensions we should consider when computing the k-neighbors	
if val.setup:
	pts = [];
	timing.log('Setup to determine adequate affinity subspace...');
	for i in range(3,config['icadim']):
		print('Finding neighbors in %i dimensions...' %(i));
		#	Slice every 10 to save time
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

subset = np.memmap('.memmapped/subset.array', dtype='float64', mode='w+', shape=(config['affdim'], numSamples));
subset[:,:] = icacoffs[ :config['affdim'], :numSamples ];
blockwidth = numSamples//config['numBlock'];

tmp_ind = np.memmap('.memmapped/tmpind.array', mode='w+', dtype='int64', shape=(blockwidth,config['numBlock']*config['n_neighbors']));
tmp_dist = np.memmap('.memmapped/tmpdist.array', mode='w+', dtype='float64', shape=(blockwidth,config['numBlock']*config['n_neighbors']));

distances = np.memmap('.memmapped/distances.array', mode='w+', dtype='float64', shape=(numSamples,config['n_neighbors']) );
indices = np.memmap('.memmapped/indices.array', mode='w+', dtype='int64', shape=(numSamples,config['n_neighbors']) );

if val.verbose: timing.log('Beginning search for %i nearest neighbors using only %i dimensions...' %(config['n_neighbors'], config['affdim']));

for i in range(config['numBlock']):

	for j in range(config['numBlock']):
		#	Would like to use symmetry to our advantage at some point, but hasn't been implemented yet.
		#if (i <= j):
			nbrs = nb(n_neighbors=config['n_neighbors'], algorithm='kd_tree').fit(subset[:,j*blockwidth:(j+1)*blockwidth].T);
			dist, ind = nbrs.kneighbors(subset[:,i*blockwidth:(i+1)*blockwidth].T);
			ind += blockwidth*j;
			tmp_ind[:, j*config['n_neighbors']:(j+1)*config['n_neighbors']] = ind;
			tmp_dist[:, j*config['n_neighbors']:(j+1)*config['n_neighbors']] = dist;
			del ind, dist, nbrs;
			update_progress( (j+i*config['numBlock']+1) / float(config['numBlock']**2) ); # +1 gets us to 100%
	indsort = np.argsort(tmp_dist,axis=1);
	
	for j in range(blockwidth):
		tmp_ind[j,:] = tmp_ind[j,indsort[j]];
		tmp_dist[j,:] = tmp_dist[j,indsort[j]];

	indices[i*blockwidth:(i+1)*blockwidth,:] = tmp_ind[:,:config['n_neighbors']];
	distances[i*blockwidth:(i+1)*blockwidth,:] = tmp_dist[:,:config['n_neighbors']];
	indices.flush(); distances.flush();
			
		#"""else:
		#	indices[i*blockwidth:(i+1)*blockwidth, j*config['n_neighbors']:(j+1)*config['n_neighbors']] = \
		#		indices[i*blockwidth:(i+1)*blockwidth, j*config['n_neighbors']:(j+1)*config['n_neighbors']]
		#"""

distances = np.memmap('.memmapped/distances.array', dtype='float64');
indices = np.memmap('.memmapped/indices.array', dtype='int64');

affdata = np.memmap('.memmapped/data.array', mode='w+', dtype='float64', shape=(numSamples*config['n_neighbors']));
affdata[:] = np.exp(-distances);
		

if val.debug and val.save:
	np.save('savefiles/indices_%s_%in.npy' %(config['pname'], config['n_neighbors']), indices[:]);

if val.verbose: timing.log('Search complete...\n\nBeginning affinity matrix generation...');

#	Affgen	-----------------------------------------------------------------------------

indptr = np.arange(0,config['n_neighbors']*numSamples+1, config['n_neighbors']);
affMat = csc_matrix( (affdata[:], indices[:], indptr), shape = (numSamples,numSamples) );

if val.save:
	save_sparse_csc('savefiles/%s_aff_%id_%in.npz' %(config['pname'], config['icadim'], config['n_neighbors']), affMat);

#if val.graph:
#	plt.spy(affMat, precision = .1, markersize = 1);
#	plt.show();

if val.verbose: timing.log('Generation complete...\n\nBeginning clustering for %i clusters...' %(config['numClusters']));

#	DNCuts	-----------------------------------------------------------------------------

	#	eig is a collection of two matrices, eig.ev = eigenvectors, eig.evl = eigenvalues
eig = dncuts(affMat, val, config, nvec=config['numClusters'], n_downsample=2, decimate=2);
if val.save:
	np.save('savefiles/%s_eigv_%ic.npy' %(config['pname'], config['numClusters']), eig.ev);
	np.save('savefiles/%s_eigval_%ic.npy' %(config['pname'], config['numClusters']), eig.evl);

if val.verbose: timing.log('Clustering complete...');

