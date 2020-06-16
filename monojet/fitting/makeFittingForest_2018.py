#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv
from glob import glob
import argparse
from pprint import pprint

parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--basedir',metavar='region',type=str,default=None)
args = parser.parse_args()
out_region = args.region
region = out_region.split('_')[0]
if region=='test':
    is_test = True 
    region = 'signal'
else:
    is_test = False
sname = argv[0]

argv=[]
import PandaAnalysis.Flat.fitting_forest as forest 
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions # kinematics
import PandaAnalysis.monojet.PandaSelection_2018 as sel
import ROOT as root

basedir = getenv('PANDA_FLATDIR')+'/'
basedir = args.basedir
outdir = getenv('PWD')+'/'
lumi = 59700.

def f(x):
    return basedir + x + '.root'

# variables to import
vmap = {}
mc_vmap = {'genBosonPt':'genBosonPt'}

nontest=True
if region in ['signal','test']:
    u,uphi = ('pfmet','pfmetphi')
elif 'pho' in region:
    mc_vmap = {'genBosonPt':'genPhotonPt'}
    u,uphi = ('pfUAmag','pfUAphi')

elif 'single' in region:
    u,uphi = ('pfUWmag','pfUWphi')
elif 'di' in region:
    u,uphi = ('pfUZmag','pfUZphi')
else:
    u,uphi = ('pfmet','pfmetphi')
    nontest=False
vmap['met'] = u 
#print 'region=', region
#print 'weights=',sel.weights[region]
weights = {'nominal' : sel.weights[region]%lumi}


# build the factory
factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts[region],
                               variables = vmap, 
                               mc_variables = mc_vmap, 
                               mc_weights = weights)

'''
# create some TChains for the V+jets
tAllW = root.TChain('events')
for f_ in ['WJets','WJets_EWK']:
    tAllW.AddFile(f(f_))
tAllZvv = root.TChain('events')
for f_ in ['ZtoNuNu','ZtoNuNu_EWK']:
    tAllZvv.AddFile(f(f_))
tAllZll = root.TChain('events')
for f_ in ['ZJets','ZJets_EWK']:
    tAllZll.AddFile(f(f_))
'''

if region == 'pho':
    factory.add_process(f('WJetsToLNu_HT'),'WJ',is_data=False,extra_weights = {'nominal':sel.weights['wjk_p']%lumi},extra_cut=sel.triggers['pho'])
    factory.add_process(f('GJets_DR-0p4_HT'),'GJ',is_data=False,extra_weights={'nominal':sel.weights['gjk']%lumi},extra_cut=sel.triggers['pho'])
# data
    factory.add_process(f('EGamma'),'Data',is_data=True,extra_cut=sel.triggers['pho'])
    vmap['met']='sqrt(pfmet*pfmet+FakePho1Pt*FakePho1Pt+2*cos(pfmetphi-FakePho1Phi)*pfmet*FakePho1Pt)'
    factory.add_process(f('EGamma'),'QCD',is_data=True, extra_weights = {'nominal':sel.weights['fkp']%lumi},replace_cut=sel.cuts['fkp'])
    vmap['met']= u

if region == 'signal':
    factory.add_process(f('WJetsToLNu_HT'),'WJ',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig']%lumi},extra_cut=sel.triggers['met'])
    factory.add_process(f('ZJetsToNuNu_HT'),'Zvv',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig']%lumi},extra_cut=sel.triggers['met'])

    factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.triggers['unblind'])
    factory.add_process(f('ZJetsToNuNu_HT-100To200'),'zh',is_data=False,extra_cut=sel.triggers['met'])
    factory.add_process(f('ZJetsToNuNu_HT-100To200'),'ggh',is_data=False,extra_cut=sel.triggers['met'])
    factory.add_process(f('ZJetsToNuNu_HT-100To200'),'vbf',is_data=False,extra_cut=sel.triggers['met'])
    factory.add_process(f('ZJetsToNuNu_HT-100To200'),'wh',is_data=False,extra_cut=sel.triggers['met'])
    factory.add_process(f('ZJetsToNuNu_HT-100To200'),'Diboson',is_data=False,extra_weights={'nominal':'0.000000001'})
    factory.add_process(f('ST'),'ST',is_data=False,extra_cut=sel.triggers['met'])
#    factory.add_process(f('JetHT'),'QCD',is_data=True,extra_weights={'nominal':'0.2*20'},replace_cut=sel.cuts['hemsim'])
    factory.add_process(f('ZJetsToNuNu_HT-100To200'),'QCD',is_data=False,extra_weights={'nominal':'0.000000001'})
if region == 'singlemuon':
    factory.add_process(f('WJetsToLNu_HT'),'WJ',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu']%lumi},extra_cut=sel.triggers['met'])

    factory.add_process(f('ST'),'ST',is_data=False,extra_cut=sel.triggers['met'])
    factory.add_process(f('DYJetsToLL_M-50_HT'),'Zll',is_data=False,extra_weights = {'nominal':sel.weights['zjk_1mu']%lumi},extra_cut=sel.triggers['met'])
    factory.add_process(f('TTJets'),'TT',is_data=False,extra_cut=sel.triggers['met'])
    factory.add_process(f('Diboson'),'Diboson',is_data=False,extra_cut=sel.triggers['met'])    
    factory.add_process(f('ST'),'QCD',is_data=False,extra_weights={'nominal':'0.000000001'})
# data
    factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.triggers['met'])
   

if region == 'singleelectron':
    factory.add_process(f('WJetsToLNu_HT'),'WJ',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e']%lumi},extra_cut=sel.triggers['ele'])
    factory.add_process(f('ST'),'ST',is_data=False,extra_cut=sel.triggers['ele'])
    factory.add_process(f('DYJetsToLL_M-50_HT'),'Zll',is_data=False,extra_weights = {'nominal':sel.weights['zjk_1e']%lumi},extra_cut=sel.triggers['ele'])
    factory.add_process(f('TTJets'),'TT',is_data=False,extra_cut=sel.triggers['ele'])
    factory.add_process(f('Diboson'),'Diboson',is_data=False,extra_cut=sel.triggers['ele'])

    factory.add_process(f('ST'),'QCD',is_data=False,extra_weights={'nominal':'0.000000001'})
    factory.add_process(f('GJets_DR-0p4_HT'),'GJ',is_data=False,extra_cut=sel.triggers['ele'])
# data
    factory.add_process(f('EGamma'),'Data',is_data=True,extra_cut=sel.triggers['ele'])

if region == 'dielectron':
    factory.add_process(f('WJetsToLNu_HT'),'WJ',is_data=False,extra_weights = {'nominal':sel.weights['wjk_2e']%lumi},extra_cut=sel.triggers['ele'])
    factory.add_process(f('ST'),'ST',is_data=False,extra_cut=sel.triggers['ele'])
    factory.add_process(f('DYJetsToLL_M-50_HT'),'Zll',is_data=False,extra_weights = {'nominal':sel.weights['zjk_2e']%lumi},extra_cut=sel.triggers['ele'])
    factory.add_process(f('TTJets'),'TT',is_data=False,extra_cut=sel.triggers['ele'])
    factory.add_process(f('Diboson'),'Diboson',is_data=False,extra_cut=sel.triggers['ele'])

# data
    factory.add_process(f('EGamma'),'Data',is_data=True,extra_cut=sel.triggers['ele'])


if region == 'dimuon':
    factory.add_process(f('WJetsToLNu_HT'),'WJ',is_data=False,extra_weights = {'nominal':sel.weights['wjk_2mu']%lumi},extra_cut=sel.triggers['met'])
    factory.add_process(f('ST'),'ST',is_data=False,extra_cut=sel.triggers['met'])
    factory.add_process(f('DYJetsToLL_M-50_HT'),'Zll',is_data=False,extra_weights = {'nominal':sel.weights['zjk_2mu']%lumi},extra_cut=sel.triggers['met'])
    factory.add_process(f('TTJets'),'TT',is_data=False,extra_cut=sel.triggers['met'])
    factory.add_process(f('Diboson'),'Diboson',is_data=False,extra_cut=sel.triggers['met'])

# data
    factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.triggers['met'])

factory.run(outdir+'/fittingForest_%s.root'%out_region)
