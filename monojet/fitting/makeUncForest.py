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
import PandaAnalysis.monojet.PandaSelection as sel
import ROOT as root

basedir = getenv('PANDA_FLATDIR')+'/'
basedir = args.basedir
outdir = getenv('PWD')+'/'
lumi = 41500.

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
weights = {'nominal' : sel.weights['nom']%lumi}


# build the factory
factory1 = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts['signal'],
                               variables = vmap, 
                               mc_variables = mc_vmap, 
                               mc_weights = weights)

factory2 = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts['singleelectron'],
                               variables = vmap, 
                               mc_variables = mc_vmap, 
                               mc_weights = weights)
factory3 = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts['singlemuon'],
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

if region == 'unctau':

    vmap['met']='pfUWmag'
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e']%lumi},extra_cut=sel.ogcut['eletau'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_tauCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_taucen']%lumi},extra_cut=sel.triggers['ele'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_tauUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_tauup']%lumi},extra_cut=sel.triggers['ele'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_tauDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_taudow']%lumi},extra_cut=sel.triggers['ele'])

    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu']%lumi},extra_cut=sel.ogcut['mettau'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_tauCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_taucen']%lumi},extra_cut=sel.triggers['met'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_tauUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_tauup']%lumi},extra_cut=sel.triggers['met'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_tauDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_taudow']%lumi},extra_cut=sel.triggers['met'])

    vmap['met']='pfmet'
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig']%lumi},extra_cut=sel.ogcut['mettau'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_tauCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_taucen']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_tauUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_tauup']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_tauDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_taudow']%lumi},extra_cut=sel.triggers['met'])

    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig']%lumi},extra_cut=sel.ogcut['mettau'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_tauCen',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_taucen']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_tauUp',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_tauup']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_tauDow',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_taudow']%lumi},extra_cut=sel.triggers['met'])

if region == 'uncele':

    vmap['met']='pfUWmag'
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e']%lumi},extra_cut=sel.triggers['ele'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_eleCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_elecen']%lumi},extra_cut=sel.triggers['ele'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_eleUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_eleup']%lumi},extra_cut=sel.triggers['ele'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_eleDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_eledow']%lumi},extra_cut=sel.triggers['ele'])

    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu']%lumi},extra_cut=sel.ogcut['metele'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_eleCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_elecen']%lumi},extra_cut=sel.triggers['met'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_eleUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_eleup']%lumi},extra_cut=sel.triggers['met'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_eleDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_eledow']%lumi},extra_cut=sel.triggers['met'])

    vmap['met']='pfmet'
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig']%lumi},extra_cut=sel.ogcut['metele'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_eleCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_elecen']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_eleUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_eleup']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_eleDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_eledow']%lumi},extra_cut=sel.triggers['met'])

    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig']%lumi},extra_cut=sel.ogcut['metele'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_eleCen',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_elecen']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_eleUp',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_eleup']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_eleDow',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_eledow']%lumi},extra_cut=sel.triggers['met'])

if region == 'uncbtag':
    vmap['met']='pfUWmag'
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e']%lumi},extra_cut=sel.ogcut['elebtag'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_bCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_bcen']%lumi},extra_cut=sel.triggers['ele'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_bUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_bup']%lumi},extra_cut=sel.triggers['ele'])
    factory2.add_process(f('WJetsToLNu_HT'),'WJ_1e_bDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1e_bdow']%lumi},extra_cut=sel.triggers['ele'])

    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu']%lumi},extra_cut=sel.ogcut['metbtag'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_bCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_bcen']%lumi},extra_cut=sel.triggers['met'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_bUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_bup']%lumi},extra_cut=sel.triggers['met'])
    factory3.add_process(f('WJetsToLNu_HT'),'WJ_1m_bDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_1mu_bdow']%lumi},extra_cut=sel.triggers['met'])

    vmap['met']='pfmet'
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig']%lumi},extra_cut=sel.ogcut['metbtag'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_bCen',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_bcen']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_bUp',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_bup']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('WJetsToLNu_HT'),'WJ_sig_bDow',is_data=False,extra_weights = {'nominal':sel.weights['wjk_sig_bdow']%lumi},extra_cut=sel.triggers['met'])
    
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig']%lumi},extra_cut=sel.ogcut['metbtag'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_bCen',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_bcen']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_bUp',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_bup']%lumi},extra_cut=sel.triggers['met'])
    factory1.add_process(f('ZJetsToNuNu_HT'),'Zvv_sig_bDow',is_data=False,extra_weights = {'nominal':sel.weights['zjk_sig_bdow']%lumi},extra_cut=sel.triggers['met'])


factory1.run(outdir+'/uncForest1_%s.root'%out_region)
factory2.run(outdir+'/uncForest2_%s.root'%out_region)
factory3.run(outdir+'/uncForest3_%s.root'%out_region)


