#!/usr/bin/env python

from sys import argv, exit
import numpy as np 
from os import getenv, system
from PandaCore.Utils.logging import logger
import PandaAnalysis.Deep.job_deep_utilities as deep_utils
from glob import glob 

NORM = True
singletons = deep_utils.singleton_branches
truth = ['jotGenPt[hbbjtidx[0]]/jotPt[hbbjtidx[0]]']
events = ['eventNumber']
fractions = {'train':0.7, 'test':0.15}
fcfg = open(argv[1])
name = argv[2]
outdir = getenv('SUBMIT_NPY')
datadir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/'
me = argv[0].split('/')[-1]
argv = []


data = {}
for fpath in fcfg.readlines():
    d = np.load(fpath.strip())
    mask = (d['jotGenPt[hbbjtidx[0]]/jotPt[hbbjtidx[0]]'] > 0.5)
    for k,v in d.iteritems():
        if v.shape[0]:
            if k not in data:
                data[k] = []
            data[k].append(v[mask])

if not len(data):
    logger.info(me, 'This was an empty config!')
    exit(0)

for k,v in data.iteritems():
    data[k] = np.concatenate(v)

if not data['pt'].shape[0]:
    logger.info(me, 'Nothing passed the mask')
    exit(0)

if NORM:
    deep_utils.normalize_arrays(data, 'pf')
    deep_utils.normalize_arrays(data, 'sv')


def reweight(x_pt):
#    x_pt = 400 + (600 * x_pt)
    return h_pt.GetBinContent(h_pt.FindBin(x_pt))
reweight = np.vectorize(reweight)

def reweight_s(x_pt):
#    x_pt = 400 + (600 * x_pt)
    return h_pt_scaled.GetBinContent(h_pt_scaled.FindBin(x_pt))
reweight_s = np.vectorize(reweight_s)


data['ptweight'] = reweight(data['rawpt'])
data['ptweight_scaled'] = reweight_s(data['rawpt'])


def dump(idx, partition):
    outpath = 'tmp/' + partition + '/' + name + '_%s.npy'

    # singletons
    d = np.vstack([data[x][idx] for x in singletons]).T 
    np.save(outpath%'singletons', d)

    # events
    d = np.vstack([data[x][idx] for x in events]).T 
    np.save(outpath%'events', d)

    # pf
    d = data['pf'][idx, :, :]
    np.save(outpath%'pf', d)

    # sv
    d = data['sv'][idx, :, :]
    np.save(outpath%'sv', d)

    # truth
    d = np.vstack([data[x][idx] for x in truth]).T 
    np.save(outpath%'truth', d)

    # pt weights
    d = data['ptweight'][idx]
    np.save(outpath%'ptweight', d)

    d = data['ptweight_scaled'][idx]
    np.save(outpath%'ptweight_scaled', d)
    


ptratio = data['jotGenPt[hbbjtidx[0]]/jotPt[hbbjtidx[0]]']
mask = np.logical_and(ptratio > 0.5, ptratio < 2)

indices = np.array(range(data['eventNumber'].shape[0]))
indices = indices[mask] # only within pT window
np.random.shuffle(indices)

N = {k:int(len(indices) * v) for k,v in fractions.iteritems()}

for d in ['train', 'test', 'validate']:
    system('mkdir -p tmp/%s'%d)

dump(indices[:N['train']], 'train')
dump(indices[N['train']:N['train']+N['test']], 'test')
dump(indices[N['train']+N['test']:], 'validate')

ret = None
for d in ['train', 'test', 'validate']:
    for ftmp in glob('tmp/'+d+'/*npy'):
        cmd = 'cp -v %s %s/%s'%(ftmp,outdir,ftmp.replace('tmp/',''))
        # cmd = 'gfal-copy -f file://$PWD/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=%s/%s'%(ftmp,outdir,ftmp.replace('tmp/',''))
        logger.info(me, cmd)
        ret = max(ret, system(cmd))

system('rm -rf tmp')
logger.debug(me, 'exit code %i'%ret)
exit(ret)
