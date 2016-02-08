from __future__ import division
import numpy as np
from MDAnalysis import *
from MDAnalysis.core.Timeseries import *
from MDAnalysis.core import Timeseries

names = ['ALA', 'ARG','ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',\
 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'];

atom = ['N-CA-CB-CG', 'N-CA-CB-SG', 'N-CA-CB-CG1', 'N-CA-CB-OG', \
'N-CA-CB-OG1', 'CA-CB-CG-CD', 'CA-CB-CG-OD1', 'CA-CB-CG-ND1', 'CA-CB-CG1-CD1', \
'CA-CB-CG-CD1', 'CA-CB-CG-SD', 'CB-CG-CD-NE', 'CB-CG-CD-OE1', 'CB-CG-CD-CE', \
'CB-CG-SD-CE', 'CG-CD-NE-CZ', 'CG-CD-CE-NZ', 'CD-NE-CZ-NH1']

chi1 = np.zeros(20);
chi1[0] = -1;
chi1[7] = -1;
chi1[4] = 1;
chi1[9] = 2;
chi1[15] = 3;
chi1[16] = 4;
chi1[-1] = 2;

chi2 = np.array([-1,5,6,6,-1,5,5,-1,7,8,9,5,10,9,5,-1,-1,9,9,-1]);

chi3 = -np.ones(20)
chi3[1] = 11;
chi3[5] = 12;
chi3[6] = 12;
chi3[11] = 13;
chi3[-8] = 14;

chi4 = -np.ones(20);
chi4[1] = 15;
chi4[11] = 16;

chi5 = -np.ones(20);
chi5[1] = 17;

global acids
acids = names;

global atomsel
atomsel = atom;

global chi
chi = np.zeros((5,20));
chi[0] = chi1;
chi[1] = chi2;
chi[2] = chi3;
chi[3] = chi4;
chi[4] = chi5;

def getAtomName(strid, chid, ind):
	
	return atomsel[chi[chid, strid].astype(int)].split('-')[ind];
	
def getStrID(resname):
	for i in range(20):
		if acids[i] == resname:
			return i;

def getChi(u, resid, ind, sliceval):
	atom = u.selectAtoms(' resid %i ' %(resid) );
	atom.set_segid('A');
	angles = []
	resName = atom.resnames()[0];
	#print resName;
	if resName == 'HIE':
		resName = 'HIS';
	strid = getStrID(resName);

	for i in range(5):
		if chi[i,strid] >= 0:
			tmp = [];
			atom = u.selectAtoms( ' atom A %i %s ' %(resid, getAtomName(strid, i, 0)), \
							  	  ' atom A %i %s ' %(resid, getAtomName(strid, i, 1)), \
							  	  ' atom A %i %s ' %(resid, getAtomName(strid, i, 2)), \
							  	  ' atom A %i %s ' %(resid, getAtomName(strid, i, 3)) 	);

			u.trajectory[int(ind*sliceval)];
			tmp.append( atom.dihedral() );

			angles.append(tmp);
		else: break;

	return angles;

if __name__ == '__main__':
	angles = [];
	ind = np.load('savefiles/hivp_clustind_31n_16c.npy');
	u = [];
	sliceval = 10;

	for i in range(5):
		u.append( Universe('hivp/hivp.pdb', 'hivp/hivp_%i.dcd' %(i+1), permissive = False) );

	trajlen = len(u[0].trajectory);

	for i in range(len(ind)):
		print 'On cluster: ', i+1;
		clust = [];
		for k in ind[i]:
			print 'On index: ', k;
			tmp = [];
			for j in range( 198 ):
				a = np.array(getChi( u[int(k // (trajlen/sliceval))] , j+1 , k % (trajlen/sliceval), sliceval ));
				tmp.append(a);
			clust.append(tmp);
		angles.append(clust);
		print len(angles);
		print len(angles[0]);
		print len(angles[0][0]);
		print angles[0][0];

	np.save('savefiles/angles.npy', angles);

	
