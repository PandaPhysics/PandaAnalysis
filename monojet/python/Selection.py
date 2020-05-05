from PandaCore.Tools.Misc import *

cuts = {}
weights = {}
triggers = {}

eventsel = ''

baseline = 'metFilter==1'

#vbf cuts
cuts['baseline'] = baseline

#regions
#cuts['pho'] = tAND(baseline,'nLooseElectron<1 && nLooseMuon<1 && nTightPhoton==1 && nLoosePhoton==1 && nTau<1 && Alt$(deltaPhi(pfmetphi,jetPhi[0]),99)>0.5 && Alt$(deltaPhi(jetPhi[1],pfmetphi),99)>0.5 && Alt$(deltaPhi(jetPhi[2],pfmetphi),99)>0.5 && Alt$(deltaPhi(jetPhi[3],pfmetphi),99)>0.5 && ((trigger&8)==8)')
cuts['pho'] = tAND(baseline,'nLooseElectron<1 && nLooseMuon<1 && nTightPhoton==1 && nLoosePhoton==1 && nTau<1 && ((trigger&8)==8)')
weights['pho'] = '%f*normalizedWeight*sf_pu'
#*photrigsf(loosePho1Pt)*GJNLO2(genPhotonPt)
triggers['pho'] = "((trigger&8)==8)"
