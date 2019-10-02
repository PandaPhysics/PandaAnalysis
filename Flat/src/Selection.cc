#include "../interface/Selection.h"
#include "../interface/Common.h"

using namespace pa; 

static float best_recoil(const GeneralTree *gt, bool include_var) {
  include_var = true;
  int maxshift = include_var ? jes2i(shiftjes::kJESTotalDown)+1 : 1;
  float max_puppi = -1, max_pf = -1;
  for (int shift = 0; shift != maxshift; ++ shift) {
    max_puppi = std::max({max_puppi, gt->puppimet[shift], 
                          gt->puppiUZmag[shift], 
                          gt->puppiUWmag[shift], 
                          gt->puppiUAmag[shift]});
    max_pf = std::max({max_pf, gt->pfmet[shift], 
                          gt->pfUZmag[shift], 
                          gt->pfUWmag[shift], 
                          gt->pfUAmag[shift]});
  }
  return std::max(max_puppi, max_pf);
}


bool LeptonSel::do_accept() const 
{ 
  if (gt->nLooseLep < 1)
    return false;

  // lepton-photon pair selection
  if (gt->nLooseElectron == 1 && gt->electronPt[0] > 25 && gt->loosePho1Pt > 25){
    TLorentzVector vLepTemp;
    vLepTemp.SetPtEtaPhiM(gt->electronPt[0],gt->electronEta[0],gt->electronPhi[0],0.000510998928);
    TLorentzVector vPhoTemp;
    vPhoTemp.SetPtEtaPhiM(gt->loosePho1Pt,gt->loosePho1Eta,gt->loosePho1Phi,0);
    if (TMath::Abs((vLepTemp+vPhoTemp).M()-91.1876) < 15)
      return true;
  }
  if (gt->nLooseMuon == 1 && gt->muonPt[0] > 25 && gt->loosePho1Pt > 25){
    TLorentzVector vLepTemp;
    vLepTemp.SetPtEtaPhiM(gt->muonPt[0],gt->muonEta[0],gt->muonPhi[0],0.10566);
    TLorentzVector vPhoTemp;
    vPhoTemp.SetPtEtaPhiM(gt->loosePho1Pt,gt->loosePho1Eta,gt->loosePho1Phi,0);
    if (TMath::Abs((vLepTemp+vPhoTemp).M()-91.1876) < 15)
      return true;
  }

  // dilepton pair selection
  if (gt->nLooseLep < 2)
    return false;
  if (gt->nLooseElectron >= 2 && gt->electronPt[0] > 20 && gt->electronPt[1] > 20)
    return true;
  if (gt->nLooseMuon >= 2 && gt->muonPt[0] > 20 && gt->muonPt[1] > 20)
    return true;
  if (gt->nLooseMuon >= 1 && gt->nLooseElectron >= 1 && gt->muonPt[0] > 20 && gt->electronPt[0] > 20)
    return true;
  return false;
}


bool LeptonFakeSel::do_accept() const 
{ 
  bool passFakeTrigger = (gt->trigger & (1<<kMuFakeTrig)) != 0 || (gt->trigger & (1<<kEleFakeTrig)) != 0;
  if (!passFakeTrigger)
    return false;
  return (gt->nLooseLep==1 ||
          (gt->nLooseLep==2 && gt->diLepMass > 70));
}


bool RecoilSel::do_accept() const
{
  float bu = best_recoil(gt, vary_jes);
  return (bu > threshold); 
}

bool FJRecoilSel::do_accept() const
{
  bool base = RecoilSel::do_accept();
  return (base && gt->fjPt[0][0]>threshold);
}


bool MonotopSel::do_accept() const
{
  bool base = RecoilSel::do_accept();
  return (base && gt->nFatJet>0 && gt->fjPt[0][0]>200);
}


bool MonohiggsSel::do_accept() const
{
  bool base = RecoilSel::do_accept();

  return (base && ((gt->hbbpt[0]>50) || (gt->nFatJet>=1 && gt->fjPt[0][0]>350 && gt->fjMSD[0][0]>30)));
}


bool VHbbSel::do_accept() const
{
  float bestMet = -1, bestJet1 = -1, bestJet2 = -1;
  int maxshift = jes2i(shiftjes::kJESTotalDown)+1;
  for (int shift = 0; shift != maxshift; ++ shift) {
    bestMet = std::max(bestMet, gt->pfmet[shift]);
  }
  maxshift = jes2i(shiftjes::N);
  for (int shift = 0; shift != maxshift; ++ shift) {
    bestJet1 = std::max(bestJet1, gt->jotPt[shift][0]);
    bestJet2 = std::max(bestJet2, gt->jotPt[shift][1]);
  }

  // ZnnHbb
  if (bestMet>150 && (gt->hbbpt[0]>50 || gt->nFatJet>0)) 
  {
    return true;
  }

  // WlnHbb
  if (((gt->nTightElectron>0 && gt->electronPt[0]>30) ||
       (gt->nTightMuon>0 && gt->muonPt[0]>25)) &&
        gt->topWBosonPt>50 &&
        (gt->hbbpt[0]>50 || gt->nFatJet>0)) 
  {
    return true;
  }

  // zllhbb
  if (((gt->nLooseElectron>1 && gt->electronPt[0]>20 && gt->electronPt[1]>15) ||
       (gt->nLooseMuon>1 && gt->muonPt[0]>20 && gt->muonPt[1]>10) ||
       (gt->nLooseMuon>0 && gt->nLooseElectron>0 && (
         (gt->electronPt[0]>25 && gt->muonPt[0]>10) ||
         (gt->electronPt[0]>10 && gt->muonPt[0]>25) 
       ))
      ) &&
      (gt->hbbpt[0]>50 || gt->nFatJet>0))
  {
    return true;
  }

  return false;
}

bool ZllHbbSel::do_accept() const
{
  float bestMet = -1, bestJet1 = -1, bestJet2 = -1;
  int maxshift = jes2i(shiftjes::kJESTotalDown)+1;
  for (int shift = 0; shift != maxshift; ++ shift) {
    bestMet = std::max(bestMet, gt->pfmet[shift]);
  }
  maxshift = jes2i(shiftjes::N);
  for (int shift = 0; shift != maxshift; ++ shift) {
    bestJet1 = std::max(bestJet1, gt->jotPt[shift][0]);
    bestJet2 = std::max(bestJet2, gt->jotPt[shift][1]);
  }

  // zllhbb
  if (((gt->nLooseElectron>1 && gt->electronPt[0]>20 && gt->electronPt[1]>15) ||
       (gt->nLooseMuon>1 && gt->muonPt[0]>20 && gt->muonPt[1]>10) ||
       (gt->nLooseMuon>0 && gt->nLooseElectron>0 && (
         (gt->electronPt[0]>25 && gt->muonPt[0]>10) ||
         (gt->electronPt[0]>10 && gt->muonPt[0]>25) 
       ))
      ) &&
      (gt->nJot[0]==1))
  {
    return true;
  }

  return false;
}


bool VHbbSelTrigger::do_accept() const
{
  if (gt->nJet[0]<2)
    return false;
  if (gt->jotPt[0][0]>40 && gt->jotPt[0][1]>25 && gt->nLooseMuon>0 && gt->muonPt[0]>25 && gt->nTightMuon>=1)
    return true;
  else
    return false;
}


bool VBFGamma::do_accept() const
{

  int jUp = jes2i(shiftjes::kJESTotalUp);
  int jDown = jes2i(shiftjes::kJESTotalDown);
  int nominal = jes2i(shiftjes::kNominal);

  if (gt->nJot[jUp]<2 && gt->nJot[jDown]<2 && gt->nJot[nominal]<2)
    return false;

  if (gt->jot12Mass[jUp]<500 && gt->jot12Mass[jDown]<500 && gt->jot12Mass[nominal]<500)
    return false;

  if (gt->pfmet[jUp]<20 && gt->pfmet[jDown]<20 && gt->pfmet[nominal]<20 && gt->puppimet[jUp]<20 && gt->puppimet[jDown]<20 && gt->puppimet[nominal]<20){
    return false;
  }
  
  if (gt->loosePho1Pt>80 || gt->alterPho1Pt>80)
    return true;

  if (gt->nLooseMuon>=1)
    if (gt->muonPt[0]>80)
      return true;

  if (gt->nLooseElectron>=1)
    if (gt->electronPt[0]>80)
      return true;
    
  return false;

  /*if (gt->nLoosePhoton<1 && (gt->nLooseMuon<1 && gt->nLooseElectron<1))
    return false;

  
  if (gt->nLoosePhoton>=1) {
    if (gt->loosePho1Pt>=80)
      return true;
    //else
      return false;
  }
  else {
    if (gt->nLooseMuon>=1) {
      if (gt->muonPt[0]>30)
	return true;
      else {
	if (gt->nLooseElectron>=1)
	  if (gt->electronPt[0]>30)
	    return true;
	  else
	    return false;
	else
	  return false;
      }
    }
    if (gt->nLooseElectron>=1) {
      if (gt->electronPt[0]>30)
	return true;
      else {
	if (gt->nLooseMuon>=1)
	  if (gt->muonPt[0]>30)
	    return true;
	  else
	    return false;
	else
	  return false;
      }
    }
  }

  
  return true;
  */
}

bool Fakerates::do_accept() const
{
  int nominal = jes2i(shiftjes::kNominal);

  if (gt->nLooseElectron>0)
    return false;

  if (gt->nLooseMuon>0)
    return false;

  if (gt->pfmet[nominal]>100)
    return false;
      
  return true;
}
