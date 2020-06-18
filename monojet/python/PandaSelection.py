from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&4)==4',
    'ele':'((trigger&1)==1 || (trigger&8)==8)',
    'pho':'(trigger&8)==8',
    'unblind':'(trigger&4)==4 && (eventNumber%5)==0',
}

ogcut = {
    'mettau':'(trigger&4)==4 && nTau==0',
    'eleele':'((trigger&1)==1 || (trigger&8)==8) && nLooseElectron==0',
    'metele':'(trigger&4)==4 && nLooseElectron==0',
    'eletau':'((trigger&1)==1 || (trigger&8)==8) && nTau==0',
    'metbtag':'(trigger&4)==4 && jetNMBtags==0',
    'elebtag':'((trigger&1)==1 || (trigger&8)==8) && jetNMBtags==0',
}

metFilter='comFilter==1'

UAcut = 'pfUAmag>250 && fabs(pfmet-calomet)/pfUAmag<0.5 && Alt$(DeltaPhi(pfUAphi,jetPhi[0]),99)>0.5 && Alt$(DeltaPhi(pfUAphi,jetPhi[1]),99)>0.5 && Alt$(DeltaPhi(pfUAphi,jetPhi[2]),99)>0.5 && Alt$(DeltaPhi(pfUAphi,jetPhi[3]),99)>0.5'
UZcut = 'pfUZmag>250 && fabs(calomet-pfmet)/pfUZmag<0.5 && Alt$(DeltaPhi(pfUZphi,jetPhi[0]),99)>0.5 && Alt$(DeltaPhi(pfUZphi,jetPhi[1]),99)>0.5 && Alt$(DeltaPhi(pfUZphi,jetPhi[2]),99)>0.5 && Alt$(DeltaPhi(pfUZphi,jetPhi[3]),99)>0.5'
UWcut = 'pfUWmag>250 && fabs(calomet-pfmet)/pfUWmag<0.5 && Alt$(DeltaPhi(pfUWphi,jetPhi[0]),99)>0.5 && Alt$(DeltaPhi(pfUWphi,jetPhi[1]),99)>0.5 && Alt$(DeltaPhi(pfUWphi,jetPhi[2]),99)>0.5 && Alt$(DeltaPhi(pfUWphi,jetPhi[3]),99)>0.5'
Ucut = 'pfmet>250  && fabs(calomet-pfmet)/pfmet<0.5 && Alt$(DeltaPhi(pfmetphi,jetPhi[0]),99)>0.5 && Alt$(DeltaPhi(pfmetphi,jetPhi[1]),99)>0.5 && Alt$(DeltaPhi(pfmetphi,jetPhi[2]),99)>0.5 && Alt$(DeltaPhi(pfmetphi,jetPhi[3]),99)>0.5'
presel = 'nTau==0 && jetNMBtags==0 && GoodLeadJet>2 && LeadJetPt > 100'
#presel = 'jetNMBtags==0 && GoodLeadJet>2 && LeadJetPt > 100'

cuts = {
    'signal'             : tAND(metFilter,tAND(presel,tAND(Ucut,'nLooseMuon==0 && nLooseElectron==0 && nLoosePhoton==0'))), 
    'singlemuon'         : tAND(metFilter,tAND(presel,tAND(UWcut,'nLoosePhoton==0 && nLooseElectron==0 && nLooseMuon==1 && nTightMuon==1 && MT(muonPt[0],muonPhi[0],pfmet,pfmetphi)<160'))),
    'singleelectron'     : tAND(metFilter,tAND(presel,tAND(UWcut,'nLoosePhoton==0 && nLooseElectron==1 && nTightElectron==1 && nLooseMuon==0 && MT(electronPt[0],electronPhi[0],pfmet,pfmetphi)<160 && pfmet > 50'))),
    'dimuon'             : tAND(metFilter,tAND(presel,tAND(UZcut,'nLooseElectron==0 && nLoosePhoton==0 && nLooseMuon==2 && nTightMuon>0 && muonPdgId[0]+muonPdgId[1]==0 &&60<Mxx(muonPt[0],muonEta[0],muonPhi[0],0.106,muonPt[1],muonEta[1],muonPhi[1],0.106) && Mxx(muonPt[0],muonEta[0],muonPhi[0],0.106,muonPt[1],muonEta[1],muonPhi[1],0.106)<120'))),
     'dielectron'         : tAND(metFilter,tAND(presel,tAND(UZcut,'nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==2 && nTightElectron>0 && electronPdgId[0]+electronPdgId[1]==0 && 60<Mxx(electronPt[0],electronEta[0],electronPhi[0],0.0005,electronPt[1],electronEta[1],electronPhi[1],0.0005) && Mxx(electronPt[0],electronEta[0],electronPhi[0],0.0005,electronPt[1],electronEta[1],electronPhi[1],0.0005)<120'))),
    'pho'             : tAND(metFilter,tAND(presel,tAND(UAcut,'nLooseMuon==0 && nLoosePhoton==1 && nTightPhoton==1 && nLooseElectron==0'))), 
   'fkp'           : tAND(metFilter,tAND(presel,'nLooseMuon==0 && nLooseElectron<1 && nFakePhoton==1 && nLoosePhoton==0 && Alt$(DeltaPhi(combinedPhi(pfmet,pfmetphi,FakePho1Pt,FakePho1Phi),jetPhi[0]),99)>0.5 && Alt$(DeltaPhi(jetPhi[1],combinedPhi(pfmet,pfmetphi,FakePho1Pt,FakePho1Phi)),99)>0.5 && Alt$(DeltaPhi(jetPhi[2],combinedPhi(pfmet,pfmetphi,FakePho1Pt,FakePho1Phi)),99)>0.5 && Alt$(DeltaPhi(jetPhi[3],combinedPhi(pfmet,pfmetphi,FakePho1Pt,FakePho1Phi)),99)>0.5 && sqrt(pfmet*pfmet+FakePho1Pt*FakePho1Pt+2*cos(pfmetphi-FakePho1Phi)*pfmet*FakePho1Pt)>250 && fabs(pfmet-calomet)/sqrt(pfmet*pfmet+FakePho1Pt*FakePho1Pt+2*cos(pfmetphi-FakePho1Phi)*pfmet*FakePho1Pt)<0.5')),
    'hemsim'       : 'nTau==0 && jetNMBtags==0 && nLooseMuon==0 && nLooseElectron==0 && nLoosePhoton==0 && comFilter==1 && (trigger&32)==32 && hemlead>0 && combinedPt(0.925*(hempt+23),hemphi,pfmet,pfmetphi)>250 && pfmet<100 && hemveto==0' 
}

#metsf vs metsf2018
weights = {
  'signal'         : '%f*sf_pu*normalizedWeight*1/fabs(mcWeight)*metsf(pfmet)',     #0.2 for unblinding
  'singlemuon'              : '%f*sf_pu*normalizedWeight*1/fabs(mcWeight)*metsf(pfUWmag)*Alt$(muonSfTight[0],1)',
  'singleelectron'              :'%f*sf_pu*normalizedWeight*1/fabs(mcWeight)*sf_eleTrig*Alt$(electronSfTight[0],1)*Alt$(electronSfReco[0],1)',
  'dimuon'            : '%f*sf_pu*normalizedWeight*1/fabs(mcWeight)*metsf(pfUZmag)*tightIDsf(Alt$(muonSelBit[0],1),Alt$(muonSfTight[0],1),Alt$(muonSfLoose[0],1))*tightIDsf(Alt$(muonSelBit[1],1),Alt$(muonSfTight[1],1),Alt$(muonSfLoose[1],1))',
  'dielectron'      : '%f*sf_pu*normalizedWeight*1/fabs(mcWeight)*sf_eleTrig*tightIDsf(Alt$(electronSelBit[0],1),Alt$(electronSfTight[0],1),Alt$(electronSfLoose[0],1))*tightIDsf(Alt$(electronSelBit[1],1),Alt$(electronSfTight[1],1),Alt$(electronSfLoose[1],1))*Alt$(electronSfReco[0],1)*Alt$(electronSfReco[1],1)',
  'pho'          : '%f*sf_pu*normalizedWeight*1/fabs(mcWeight)*sfphonew(loosePho1Eta)*photrigsf(loosePho1Pt)',
  'fkp'       :  '%f/41500.*fr_pho(FakePho1Pt)*(nTau==0)',
  'gjk'      : 'sfgjets(genPhotonPt)',
  'wjk_2mu'      : 'sfwjets(genBosonPt)',
  'wjk_1mu'      : 'sfwjets(genBosonPt)',
  'wjk_2e'      : 'sfwjets(genBosonPt)',
  'wjk_1e'      : 'sfwjets(genBosonPt)',
  'zjk_2mu'      : 'sfzjets(genBosonPt)',
  'zjk_1mu'      : 'sfzjets(genBosonPt)',
  'zjk_2e'      : 'sfzjets(genBosonPt)',
  'zjk_1e'      : 'sfzjets(genBosonPt)',
  'wjk_p'     : 'sfwjets(genBosonPt)',
  'wjk_sig'   : 'sfwjets(genBosonPt)',
  'zjk_sig' : 'sfzjets(genBosonPt)',
  'wjk_sig_taucen': 'Alt$(1-tauSF(tauPt[0],0),1)*Alt$(1-tauSF(tauPt[1],0),1)',
  'wjk_sig_tauup': 'Alt$(1-tauSF(tauPt[0],1),1)*Alt$(1-tauSF(tauPt[1],1),1)',
  'wjk_sig_taudow': 'Alt$(1-tauSF(tauPt[0],2),1)*Alt$(1-tauSF(tauPt[1],2),1)',
  'wjk_sig_elecen': '(1-Alt$(electronSfLoose[0]*electronSfReco[0],0))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0))',
  'wjk_sig_eleup': '(1-Alt$(electronSfLoose[0]*electronSfReco[0],0)*(1+Alt$(eleunc(electronPt[0],electronEta[0]),0)))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0)*(1+Alt$(eleunc(electronPt[1],electronEta[1]),0)))',
  'wjk_sig_eledow': '(1-Alt$(electronSfLoose[0]*electronSfReco[0],0)*(1-Alt$(eleunc(electronPt[0],electronEta[0]),0)))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0)*(1-Alt$(eleunc(electronPt[1],electronEta[1]),0)))',
  'zjk_sig_taucen': 'Alt$(1-tauSF(tauPt[0],0),1)*Alt$(1-tauSF(tauPt[1],0),1)',
  'zjk_sig_tauup': 'Alt$(1-tauSF(tauPt[0],1),1)*Alt$(1-tauSF(tauPt[1],1),1)',
  'zjk_sig_taudow': 'Alt$(1-tauSF(tauPt[0],2),1)*Alt$(1-tauSF(tauPt[1],2),1)',
  'zjk_sig_elecen':'(1-Alt$(electronSfLoose[0]*electronSfReco[0],0))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0))',
  'zjk_sig_eleup': '(1-Alt$(electronSfLoose[0]*electronSfReco[0],0)*(1+Alt$(eleunc(electronPt[0],electronEta[0]),0)))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0)*(1+Alt$(eleunc(electronPt[1],electronEta[1]),0)))',
  'zjk_sig_eledow': '(1-Alt$(electronSfLoose[0]*electronSfReco[0],0)*(1-Alt$(eleunc(electronPt[0],electronEta[0]),0)))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0)*(1-Alt$(eleunc(electronPt[1],electronEta[1]),0)))',
  'wjk_1e_taucen':'Alt$(1-tauSF(tauPt[0],0),1)*Alt$(1-tauSF(tauPt[1],0),1)',
  'wjk_1e_tauup': 'Alt$(1-tauSF(tauPt[0],1),1)*Alt$(1-tauSF(tauPt[1],1),1)',
  'wjk_1e_taudow': 'Alt$(1-tauSF(tauPt[0],2),1)*Alt$(1-tauSF(tauPt[1],2),1)',
  'wjk_1e_elecen':'1',
  'wjk_1e_eleup': '1+Alt$(eleunc(electronPt[0],electronEta[0],1),0)',
  'wjk_1e_eledow': '1-Alt$(eleunc(electronPt[0],electronEta[0],1),0)',
  'wjk_1mu_taucen':'Alt$(1-tauSF(tauPt[0],0),1)*Alt$(1-tauSF(tauPt[1],0),1)',
  'wjk_1mu_tauup': 'Alt$(1-tauSF(tauPt[0],1),1)*Alt$(1-tauSF(tauPt[1],1),1)',
  'wjk_1mu_taudow': 'Alt$(1-tauSF(tauPt[0],2),1)*Alt$(1-tauSF(tauPt[1],2),1)',
  'wjk_1mu_elecen': '(1-Alt$(electronSfLoose[0]*electronSfReco[0],0))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0))',
  'wjk_1mu_eleup': '(1-Alt$(electronSfLoose[0]*electronSfReco[0],0)*(1+Alt$(eleunc(electronPt[0],electronEta[0]),0)))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0)*(1+Alt$(eleunc(electronPt[1],electronEta[1]),0)))',
  'wjk_1mu_eledow': '(1-Alt$(electronSfLoose[0]*electronSfReco[0],0)*(1-Alt$(eleunc(electronPt[0],electronEta[0]),0)))*(1-Alt$(electronSfLoose[1]*electronSfReco[1],0)*(1-Alt$(eleunc(electronPt[1],electronEta[1]),0)))',
  'wjk_sig_bcen': '(1-Alt$(sf_mb[0],0))*(1-Alt$(sf_mb[1],0))*(1-Alt$(sf_mb[2],0))*(1-Alt$(sf_mb[3],0))',
  'wjk_sig_bup': '(1-Alt$(sf_mb_up[0],0))*(1-Alt$(sf_mb_up[1],0))*(1-Alt$(sf_mb_up[2],0))*(1-Alt$(sf_mb_up[3],0))',
  'wjk_sig_bdow': '(1-Alt$(sf_mb_dow[0],0))*(1-Alt$(sf_mb_dow[1],0))*(1-Alt$(sf_mb_dow[2],0))*(1-Alt$(sf_mb_dow[3],0))',
  'zjk_sig_bcen': '(1-Alt$(sf_mb[0],0))*(1-Alt$(sf_mb[1],0))*(1-Alt$(sf_mb[2],0))*(1-Alt$(sf_mb[3],0))',
  'zjk_sig_bup': '(1-Alt$(sf_mb_up[0],0))*(1-Alt$(sf_mb_up[1],0))*(1-Alt$(sf_mb_up[2],0))*(1-Alt$(sf_mb_up[3],0))',
  'zjk_sig_bdow': '(1-Alt$(sf_mb_dow[0],0))*(1-Alt$(sf_mb_dow[1],0))*(1-Alt$(sf_mb_dow[2],0))*(1-Alt$(sf_mb_dow[3],0))',
  'wjk_1e_bcen': '(1-Alt$(sf_mb[0],0))*(1-Alt$(sf_mb[1],0))*(1-Alt$(sf_mb[2],0))*(1-Alt$(sf_mb[3],0))',
  'wjk_1e_bup': '(1-Alt$(sf_mb_up[0],0))*(1-Alt$(sf_mb_up[1],0))*(1-Alt$(sf_mb_up[2],0))*(1-Alt$(sf_mb_up[3],0))',
  'wjk_1e_bdow': '(1-Alt$(sf_mb_dow[0],0))*(1-Alt$(sf_mb_dow[1],0))*(1-Alt$(sf_mb_dow[2],0))*(1-Alt$(sf_mb_dow[3],0))',
  'wjk_1mu_bcen': '(1-Alt$(sf_mb[0],0))*(1-Alt$(sf_mb[1],0))*(1-Alt$(sf_mb[2],0))*(1-Alt$(sf_mb[3],0))',
  'wjk_1mu_bup': '(1-Alt$(sf_mb_up[0],0))*(1-Alt$(sf_mb_up[1],0))*(1-Alt$(sf_mb_up[2],0))*(1-Alt$(sf_mb_up[3],0))',
  'wjk_1mu_bdow': '(1-Alt$(sf_mb_dow[0],0))*(1-Alt$(sf_mb_dow[1],0))*(1-Alt$(sf_mb_dow[2],0))*(1-Alt$(sf_mb_dow[3],0))',
  'nom'       : '%f*sf_pu*normalizedWeight'
}


weights['gjk'] = tTIMES(weights['gjk'],weights['pho'])
weights['wjk_2mu'] = tTIMES(weights['wjk_2mu'],weights['dimuon'])
weights['wjk_1mu'] = tTIMES(weights['wjk_1mu'],weights['singlemuon'])
weights['wjk_2e'] = tTIMES(weights['wjk_2e'],weights['dielectron'])
weights['wjk_1e'] = tTIMES(weights['wjk_1e'],weights['singleelectron'])
weights['wjk_p'] = tTIMES(weights['wjk_p'],weights['pho'])
weights['wjk_sig'] = tTIMES(weights['wjk_sig'],weights['signal'])
weights['zjk_2mu'] = tTIMES(weights['zjk_2mu'],weights['dimuon'])
weights['zjk_1mu'] = tTIMES(weights['zjk_1mu'],weights['singlemuon'])
weights['zjk_2e'] = tTIMES(weights['zjk_2e'],weights['dielectron'])
weights['zjk_1e'] = tTIMES(weights['zjk_1e'],weights['singleelectron'])
weights['zjk_sig'] = tTIMES(weights['zjk_sig'],weights['signal'])

weights['zjk_sig_taucen'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_taucen'])
weights['zjk_sig_tauup'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_tauup'])
weights['zjk_sig_taudow'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_taudow'])
weights['zjk_sig_elecen'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_elecen'])
weights['zjk_sig_eleup'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_eleup'])
weights['zjk_sig_eledow'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_eledow'])

weights['wjk_sig_taucen'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_taucen'])
weights['wjk_sig_tauup'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_tauup'])
weights['wjk_sig_taudow'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_taudow'])
weights['wjk_sig_elecen'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_elecen'])
weights['wjk_sig_eleup'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_eleup'])
weights['wjk_sig_eledow'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_eledow'])

weights['wjk_1e_taucen'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_taucen'])
weights['wjk_1e_tauup'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_tauup'])
weights['wjk_1e_taudow'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_taudow'])
weights['wjk_1e_elecen'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_elecen'])
weights['wjk_1e_eleup'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_eleup'])
weights['wjk_1e_eledow'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_eledow'])

weights['wjk_1mu_taucen'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_taucen'])
weights['wjk_1mu_tauup'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_tauup'])
weights['wjk_1mu_taudow'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_taudow'])
weights['wjk_1mu_elecen'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_elecen'])
weights['wjk_1mu_eleup'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_eleup'])
weights['wjk_1mu_eledow'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_eledow'])


weights['wjk_1mu_bcen'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_bcen'])
weights['wjk_1mu_bup'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_bup'])
weights['wjk_1mu_bdow'] = tTIMES(weights['wjk_1mu'],weights['wjk_1mu_bdow'])
weights['wjk_1e_bcen'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_bcen'])
weights['wjk_1e_bup'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_bup'])
weights['wjk_1e_bdow'] = tTIMES(weights['wjk_1e'],weights['wjk_1e_bdow'])
weights['wjk_sig_bcen'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_bcen'])
weights['wjk_sig_bup'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_bup'])
weights['wjk_sig_bdow'] = tTIMES(weights['wjk_sig'],weights['wjk_sig_bdow'])
weights['zjk_sig_bcen'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_bcen'])
weights['zjk_sig_bup'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_bup'])
weights['zjk_sig_bdow'] = tTIMES(weights['zjk_sig'],weights['zjk_sig_bdow'])



#for x in ['dimuon','dielectron','singlemuon','singleelectron']:
#    if 'electron' in x:
##      if 'di' in x:
#        weights[x] = tTIMES(weights['z'], 'sf_eleTrig')
