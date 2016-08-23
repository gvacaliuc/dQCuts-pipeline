from __future__ import division
import numpy as np
import argparse
import sys
from wqaa.main import main as wqaa_main
from qaa_simgen.main import main as simgen_main
from dncuts_eigensolver.main import main as dncuts_main
from scipy.sparse import *

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

def main(config):

    #   Init
    
    #   QAA
    config['numblock'] = config['aff_numblock'];
    wqaa_main(config);

    #   Similarity Matrix Generation
    config['icafile'] = os.path.join(config['saveDir'],'{0}_icacoffs_{1}dim.array'.format(config['pname'],config['icaDim']));
    simgen_main(config);

    #   Clustering
    config['numblock'] = config['mult_numblock'];
    config['aff'] = os.path.join(saveDir, 
                                '{0}_aff_{1}d_{2}n.npz'.format(config['pname'], 
                                                               icacoffs.shape[0],
                                                               config['n_neighbors'],
                                                              ),
                                );
    dncuts_main(config);

def validate(config):

    required = ['numblock', 'logfile', 'saveDir', 'figDir','pname', 
                'n_neighbors', 'affdim', 'analysis', 'pdbfile', 'startRes',
                'endRes',];
    for field in required:
        if not field in config:
            raise ValueError('You didn\'t provide a value for the field: {0}'.format(field));

    directories = ['saveDir', 'figDir',];
    for directory in directories:
        if not os.path.isdir(config[directory]):
            os.makedirs(config[directory]);

    if not os.path.isdir(os.path.join(config['saveDir'], '.memmapped')):
        os.makedirs(os.path.join(config['saveDir'], '.memmapped'));
    if not os.path.isdir(os.path.join(config['saveDir'], '.debug',)):
        os.makedirs(os.path.join(config['saveDir'], '.debug',));
    if not os.path.isdir(os.path.join(config['saveDir'], '.memmapped', 'spmultiply')):
        os.makedirs(os.path.join(config['saveDir'], '.memmapped', 'spmultiply'));
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', action='store_true', dest='graph', default=False, help='Shows graph of Affinity Matrix and eigenv\'s, depending on flags.')
    parser.add_argument('-v', action='store_true', dest='verbose', default=False, help='Runs program verbosely.')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Prints debugging help.')
    parser.add_argument('--setup', action='store_true', dest='setup', default=False, help='Runs setup calculations: Cum. Sum. of cov. spectrum\nand unit radius neighbor search.');
    parser.add_argument('--config', type=str, dest='configpath', default='config.yaml',
                        help='Input other configuration file.');
    values = parser.parse_args()

    #   Get config from file
    with open(values.configpath) as f:
        conf_file = f.read();
        config = yaml.load(conf_file);
    if not 'config' in locals(): raise IOError(
    'Issue opening and reading configuration file: {0}'.format(os.path.abspath(values.configpath)) );

    validate(config);

    #   Update config with CLARGS
    level = 30;
    if values.verbose: level = 20;
    elif values.debug: level = 10;
    config['graph'] = values.graph;
    config['setup'] = values.setup;

    #   Setup stream logger
    ch = logging.StreamHandler(sys.stdout);
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s');
    ch.setLevel(level);
    ch.setFormatter(formatter);

    log.addHandler(ch);

    log.debug('Configuration File:\n'+conf_file);
    log.info('Using Configuration File: {0}'.format(os.path.abspath(values.configpath)));

    main(config);
