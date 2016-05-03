from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time


def histgen( angles ):

	numclust = len(angles);
	hist = [];
	tmphist = [];
	anglehist = [];
	
	for i in range(numclust):

		tmphist = [];
		for k in range( len(angles[i][0]) ):
		
			for l in range( len(angles[i][0][k]) ):
			
				anglehist = [];
				for j in range( len(angles[i]) ):
				
					anglehist.append( angles[i][j][k][l] );
				
				histogram = np.histogram( anglehist, bins = 10, range=(-180,180) );
				tmphist.append( histogram[0]/( np.sum(histogram[0]) ) );
				print histogram[0], np.sum(histogram[0]), histogram[0]/( np.sum(histogram[0]) );
		hist.append( tmphist );

	np.save('savefiles/hist.npy', hist);

def histgraph( hist ):

	plt.ion();
	fig = plt.figure();
	ax = fig.add_subplot(111);
	
	for i in range(370):
		
		print "Chi #%i" %(i+1);
		for j in range(16):
			ax.plot( hist[j, i] );
		plt.draw();
		plt.show();
		time.sleep(1);
		plt.cla()

if __name__ == '__main__':
	#histgen( np.load('savefiles/angles.npy') );
	histgraph( np.load('hist.npy') );
