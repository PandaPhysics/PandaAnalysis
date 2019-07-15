#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector 
import PandaCore.Tools.Functions
import numpy as np
import json 
import os

args = parse('--out', '--name', '--json')

try:
    os.makedirs(args.out)
except OSError:
    pass

with open(args.json) as jsonfile:
    payload = json.load(jsonfile)
    weight = payload['weight']
    basedir = payload['base']
    features = payload['features']
    cut = payload['cut']
    for i,s in enumerate(payload['samples']):
        if s['name'] == args.name:
            samples = s['samples']
            y = i
            break
    else:
        logger.error(sys.argv[0], 'Could not identify process '+args.name)
        sys.exit(1)

s = Selector()
chain = root.TChain('events')
for sample in samples:
    chain.AddFile(basedir + '/' + sample + '.root')

logger.info(sys.argv[0], 'Reading files for process '+args.name)
s.read_tree(chain, branches=(features+[weight]), cut=cut)

X = np.vstack([s[f] for f in features]).T 
W = s[weight]
#W *= 1000 / W.sum()
Y = y * np.ones(shape=W.shape)

def save(arr, label, additional_label=None):
    if additional_label == None:
        fout = args.out+'/'+args.name+'_'+label+'.npy'
    else:
        fout = args.out+'/'+args.name+'_'+additional_label+'_'+label+'.npy'
    np.save(fout, arr)
    logger.info(sys.argv[0], 'Saved to '+fout)

save(X, 'x')
save(Y, 'y')
save(W, 'w')

#save(X, 'x', 'with_fj')
#save(Y, 'y', 'with_fj')
#save(W, 'w', 'with_fj')
