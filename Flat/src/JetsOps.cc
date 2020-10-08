#include "../interface/JetsOps.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"
#include "PandaAnalysis/Utilities/interface/NeutrinoSolver.h"
#include "TSystem.h"

using namespace pa;
using namespace std;
using namespace panda;
namespace fj = fastjet;
using JECParams = JetCorrectorParameters;

inline float centralOnly(float x, float aeta, float def = -1)
{
  return  aeta < 2.4 ? x : -1;
}

void BaseJetOp::shiftMET(const RecoMet& met, TLorentzVector& v, shiftjes shift, bool metjer)
{
  float pt;
  float phi;
  TLorentzVector temp2;
  temp2.SetPtEtaPhiM(met.pt,0,met.phi,0);
  if(metjer) {temp2=temp2+met_jer_corr;}
  switch (shift) {
  case shiftjes::kNominal:
    pt = temp2.Pt();
    phi = temp2.Phi();
    break;
  case shiftjes::kJESTotalUp:
    pt=(temp2+met_corr).Pt();
    phi=(temp2+met_corr).Phi();
    break;
  case shiftjes::kJESTotalDown:
    pt=(temp2-met_corr).Pt();
    phi=(temp2-met_corr).Phi();
    break;
  default:
    logger.error("shiftMET", "Unknown JES type!");
    exit(1);
  }

  v.SetPtEtaPhiM(pt, 0, phi, 0);
}

void BaseJetOp::lowptshift(const CorrT1METJet& jet, bool smear)
{
  float pt = jet.pt;
  if (smear&&pt>0) {
      double smearFac=1, smearFacUp=1, smearFacDown=1;
      double tempt, temeta,temphi;
      tempt=-1;
      double minR=10;
      for (auto& gjet: event.GenJet){
        temeta=gjet.eta; temphi=gjet.phi;
        double dR=DeltaR2(temeta,temphi,jet.eta,jet.phi);
        if(dR<minR) {minR=dR; tempt=gjet.pt;}
      }
      if (minR<0.09){
        jer->getStochasticSmear(pt,jet.eta,event.fixedGridRhoFastjetAll,smearFac,smearFacUp,smearFacDown,1,tempt, analysis.year);
      }
      pt *= smearFac;

      if(jet.pt>15 && fabs(jet.eta)<5.191){
        TLorentzVector temp3,temp5;
        temp3.SetPtEtaPhiM(pt,0, jet.phi,0);
        temp5.SetPtEtaPhiM(pt/smearFac,0,jet.phi,0);
        met_jer_corr=met_jer_corr+temp5-temp3;
       }

    }
}

JetWrapper BaseJetOp::shiftJet(const Jet& jet, shiftjes shift, bool smear)
{
  float pt = jet.pt;
  //smear = false;
  if (smear) {
      double smearFac=1, smearFacUp=1, smearFacDown=1;
      double tempt, temeta,temphi;
      tempt=-1;
      double minR=10;
      for (auto& gjet: event.GenJet){
        temeta=gjet.eta; temphi=gjet.phi;
        double dR=DeltaR2(temeta,temphi,jet.eta,jet.phi);
        if(dR<minR) {minR=dR; tempt=gjet.pt;}
      } 
      if (minR<0.09){
	jer->getStochasticSmear(pt,jet.eta,event.fixedGridRhoFastjetAll,smearFac,smearFacUp,smearFacDown,1,tempt, analysis.year);
      }
      else {
	jer->getStochasticSmear(pt,jet.eta,event.fixedGridRhoFastjetAll,smearFac,smearFacUp,smearFacDown,0,-99, analysis.year);
      }
      pt *= smearFac;
      int njer=gt.nJer;
      if(shift == shiftjes::kNominal && jet.pt>15 && njer<NJET){
         gt.jesPt[njer]=jet.pt;
         gt.jerPt[njer]=pt;
         gt.jerEta[njer]=jet.eta;
         gt.jerPhi[njer]=jet.phi;
         gt.jercef[njer]=jet.chEmEF;
         gt.jernef[njer]=jet.neEmEF;
         gt.nJer++;
        
      } 
      if(shift == shiftjes::kNominal && jet.pt>15 && fabs(jet.eta)<5.191 && jet.chEmEF + jet.neEmEF < 0.9){
        TLorentzVector temp3,temp5,temp6,temp8;
        temp3.SetPtEtaPhiM(pt,0, jet.phi,0);
        temp5.SetPtEtaPhiM(pt/smearFac,0,jet.phi,0);
        temp6.SetPtEtaPhiM(pt*jet.l1fac,0, jet.phi,0);
        temp8.SetPtEtaPhiM(pt*jet.l1fac/smearFac,0,jet.phi,0);
        met_jer_corr=met_jer_corr+temp5-temp3+temp6-temp8;     
      }
  }
  if (shift != shiftjes::kNominal && analysis.rerunJES && !analysis.isData){
    int ishift = jes2i(shift);
    bool isUp = !(ishift % 2 == 0);
      (*scaleUnc)[ishift]->setJetPt(pt);
      (*scaleUnc)[ishift]->setJetEta(jet.eta);
      double relShift = (*scaleUnc)[ishift]->getUncertainty(isUp);
      if (!isUp)
        relShift = -relShift; 
      pt *= (1 + relShift);
      if(shift == shiftjes::kJESTotalUp){
        TLorentzVector temp2,temp4;
        temp2.SetPtEtaPhiM(pt,0, jet.phi,0);
        temp4.SetPtEtaPhiM(pt/(1+relShift),0,jet.phi,0);
        met_corr+=temp4-temp2;
      }
  }
  return JetWrapper(pt, jet);
}


void BaseJetOp::do_readData(TString dirPath)
{
  if (recalcJER) {
    jer.reset(new JERReader(dirPath+"/jec/"+jerV+"/"+jerV+"_MC_SF_"+jetType+".txt",
                            dirPath+"/jec/"+jerV+"/"+jerV+"_MC_PtResolution_"+jetType+".txt", analysis.year));
  }

  if (!analysis.rerunJES)
    return;

  TString jecVFull = jecReco+spacer+jecV;

  TString basePath = dirPath+"/jec/"+jecVFull+"/"+campaign+"_"+jecVFull;

  vector<JECParams> params = {
    JECParams((basePath+"_MC_L1FastJet_"+jetType+".txt").Data()),
    JECParams((basePath+"_MC_L2Relative_"+jetType+".txt").Data()),
    JECParams((basePath+"_MC_L3Absolute_"+jetType+".txt").Data()),
    JECParams((basePath+"_MC_L2L3Residual_"+jetType+".txt").Data())
  };
  scales["MC"].reset(new FactorizedJetCorrector(params));
  for (auto e : eraGroups) {
    basePath = dirPath+"/jec/"+jecVFull+"/"+campaign+"_"+jecReco+e+spacer+jecV;
    params = {
      JECParams((basePath+"_DATA_L1FastJet_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L2Relative_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L3Absolute_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L2L3Residual_"+jetType+".txt").Data())
    };
    scales["data"+e].reset(new FactorizedJetCorrector(params));
  }

  basePath = dirPath+"/jec/"+jecVFull+"/"+campaign+"_";
  setScaleUnc("MC", (basePath+jecVFull+"_MC_UncertaintySources_"+jetType+".txt").Data());
  for (auto e : eraGroups) {
    setScaleUnc("data"+e, (basePath+jecReco+e+spacer+jecV+"_DATA_UncertaintySources_"+jetType+".txt").Data());
  }
}


void BaseJetOp::setScaleUnc(TString tag, TString path)
{
  scaleUncs[tag] = std::vector<std::shared_ptr<JetCorrectionUncertainty>>(jes2i(shiftjes::N),nullptr);
  JESLOOP {
    if (shift % 2 == 0)
      continue;
    TString shiftName = jesName(i2jes(shift));
    shiftName.ReplaceAll("JES","");
    shiftName = shiftName(0,shiftName.Length()-2);
    scaleUncs[tag][shift] = make_shared<JetCorrectionUncertainty>(JECParams(path.Data(), shiftName.Data()));
    scaleUncs[tag][shift+1] = scaleUncs[tag][shift]; // Down = Up
  }
}

void JetOp::setupJES()
{
  if (!analysis.rerunJES || (scaleUnc != nullptr))
    return;

  if (analysis.isData) {
    TString thisEra = utils.eras->getEra(gt.runNumber);
    for (auto& iter : scaleUncs) {
      if (!iter.first.Contains("data"))
        continue;
      if (iter.first.Contains(thisEra)) {
        scaleUnc = &(scaleUncs[iter.first]);
        scale = scales[iter.first].get();
        return;
      }
    }
  } else {
    scaleUnc = &(scaleUncs["MC"]);
    scale = scales["MC"].get();
  }
}

void JetOp::varyJES()
{
    met_corr.SetPtEtaPhiM(0,0,0,0);
    met_jer_corr.SetPtEtaPhiM(0,0,0,0);

  if(analysis.LowPtJet){
    for(auto &ljet: event.CorrT1METJet){
      lowptshift(ljet,(recalcJER) && !analysis.isData);
    }
  }
 
  JESLOOP {
    auto& jets = (*jesShifts)[shift];
    jets.reserve(ak4Jets->size());
   
    std::vector<JetWrapper> all_presorted;
    all_presorted.reserve(ak4Jets->size());
    for (auto &j : *ak4Jets) {
      all_presorted.push_back(shiftJet(j, i2jes(shift), (recalcJER && !analysis.isData))); 
    }

    std::sort(all_presorted.begin(), all_presorted.end(),
	      [](const JetWrapper x, const JetWrapper y) { return x.pt > y.pt; });
    jets.all = all_presorted;    

  }

  for (size_t iJ = 0; iJ != (*jesShifts)[0].all.size(); ++iJ) {
    auto* nominal = &((*jesShifts)[0].all[iJ]);
    nominal->maxpt = 0;
    JESLOOP {
      auto* jet = &((*jesShifts)[shift].all[iJ]);
      jet->nominal = nominal;
      if (jet->pt > nominal->maxpt)
        nominal->maxpt = jet->pt;
    }
  }
}


void JetOp::do_execute()
{

  setupJES();

  varyJES();

  float maxJetEta = 4.7;
  int nJetDPhi = 50;
  float minJetPt = 20;

  JESLOOP {
    bool isNominal = (shift == jes2i(shiftjes::kNominal));
    bool metShift = (i2jes(shift) <= shiftjes::kJESTotalDown);
    JESHandler& jets = (*jesShifts)[shift];

    float max_pt=0;
    (*currentJES) = &jets;
    for (auto& jw : jets.all) {
      (*currentJet) = &jw;
      auto& jet = jw.get_base();
      float aeta = abs(jet.eta);
      float pt = jw.pt;
      if (analysis.year == 2016 || analysis.year == 2017) {
        if (isNominal && !isMatched(matchVeryLoosePhos.get(),0.16,jet.eta,jet.phi)) {
          // prefiring weights
          gt.sf_l1Prefire *= (1.0 - utils.getCorr(cL1PreFiring,jet.eta,pt));
        }
      }

      if (aeta > maxJetEta || jw.nominal->maxpt < minJetPt)
        continue;    
       
      if (isNominal) { // perform cleaning only on the nominal jet and save flags
        if (isMatched(matchLeps.get(),0.16,jet.eta,jet.phi))
          {jw.isLep = true;}
        if (!analysis.hbb && isMatched(matchPhos.get(),0.16,jet.eta,jet.phi))
          jw.isPho = true;
        if (analysis.hbb && jet.puId < utils.getCorr(cJetLoosePUID, aeta,min(39.99f,pt)))
          jw.isPileupJet = true;
      }

      if (jw.nominal->isLep || jw.nominal->isPho || jw.nominal->isPileupJet)
        continue;
     
      bool is_tight = false;
      is_tight=((jet.jetId>>1)%2==1);
      if(!is_tight && (analysis.year == 2017 || analysis.year == 2018))
        continue;

      if (jw.nominal->maxpt > minJetPt) {
        jets.cleaned.push_back(&jw);

        if (jets.cleaned.size() < 3) {
          if (utils.getCorr(cBadECALJets, jet.eta, jet.phi) > 0)
            gt.badECALFilter = 0;
        }

	int njet = jets.cleaned.size() - 1;
	
	if (pt > minJetPt && gt.nJet[shift]<NJET){
	  int nj=gt.nJet[shift];
	  gt.jetPt[shift][nj] = pt;
          if (isNominal) {
	    jw.cleaned_idx = njet; 
	    gt.jetEta[nj] = jet.eta;
	    gt.jetPhi[nj] = jet.phi;
	    gt.jetM[shift][nj] = jet.mass;
          }
	  gt.nJet[shift]++;         
       }
        if(gt.nJetMax<gt.nJet[shift]) gt.nJetMax=gt.nJet[shift];

        TLorentzVector vJet; vJet.SetPtEtaPhiM(pt, jet.eta, jet.phi, jet.mass);
        if (metShift && njet < nJetDPhi && pt > cfg.minJetPt) { // only do this for fully-correlated shifts
          gt.dphipfmet[shift]    = min(fabs(vJet.DeltaPhi(jets.vpfMET)),    (double)gt.dphipfmet[shift]);
          if (analysis.recoil) {
            gt.dphipfUA[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUA)), (double)gt.dphipfUA[shift]);
            gt.dphipfUW[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUW)), (double)gt.dphipfUW[shift]);
            gt.dphipfUZ[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUZ)), (double)gt.dphipfUZ[shift]);
          }
        }
      }
    }

    if (metShift) {
      switch (gt.whichRecoil) {
        case 0: // MET
          gt.dphipfU[shift] = gt.dphipfmet[shift];
          break;
        case -1: // photon
          gt.dphipfU[shift] = gt.dphipfUA[shift];
          break;
        case 1:
          gt.dphipfU[shift] = gt.dphipfUW[shift];
          break;
        case 2:
          gt.dphipfU[shift] = gt.dphipfUZ[shift];
          break;
        default: // c'est impossible !
          break;
      }
    }

    if (metShift) {
      jets.sort();
    }


    /*
    gt.avgDPhi[shift] = 0;
    std::vector<TLorentzVector> allVecs;
    TLorentzVector vJet(0,0,0,0); 
    for (int m = 0; m<gt.nJet[shift]; m++){
      TLorentzVector tmpJet(0,0,0,0); 
      tmpJet.SetPtEtaPhiM(gt.jetPt[shift][m], gt.jetEta[m], gt.jetPhi[m], gt.jetM[m]);
      vJet +=  tmpJet;
      allVecs.push_back(tmpJet);
      for (int n = m+1; n<gt.nJet[shift]; n++){
	gt.avgDPhi[shift] += 1./gt.nJet[shift]*SignedDeltaPhi(gt.jetPhi[m],gt.jetPhi[n]);
      }
    }
    gt.jetTotM[shift] = vJet.M();
    gt.sphericity[shift] = sphericity(2.,allVecs);
    */

  } // shift loop

  METLOOP {
    auto& jets = (*jesShifts)[shift];
    if (!analysis.puppiMet && analysis.year != 2017){
      shiftMET(event.MET, jets.vpfMET, i2jes(shift),analysis.METJER);
      gt.pfmet[shift] = jets.vpfMET.Pt();
      gt.pfmetphi[shift] = jets.vpfMET.Phi();
    }
    else if( !analysis.puppiMet && analysis.year == 2017){
      shiftMET(event.METFixEE2017, jets.vpfMET, i2jes(shift),analysis.METJER);
      gt.pfmet[shift] = jets.vpfMET.Pt();
      gt.pfmetphi[shift] = jets.vpfMET.Phi();
    }
    else{
      shiftMET(event.PuppiMET, jets.vpuppiMET, i2jes(shift),analysis.METJER);
      gt.puppimet[shift] = jets.vpuppiMET.Pt();
      gt.puppimetphi[shift] = jets.vpuppiMET.Phi();
    }

    jets.vpfMETNoMu.SetMagPhi(gt.pfmet[shift], gt.pfmetphi[shift]);
  }

     for (auto &gen : event.GenJet) {
            bool jetpass=true;
            if(gen.pt<30 || fabs(gen.eta)>2.4) continue;
            if (DeltaR2(gen.eta, gen.phi, gt.genPhotonEta, gt.genPhotonPhi) < 0.16 && gt.genPhotonPt>0) jetpass=false;
            for(auto &genP : event.GenPart){
               if((fabs(genP.pdgId)==11 || fabs(genP.pdgId)==13) && (genP.status)==1){
                   if (DeltaR2(gen.eta, gen.phi, genP.eta, genP.phi) < 0.16) jetpass=false;
                 }
            }          
         if(jetpass){ 
           gt.GenJetPt[gt.nGenJet]=gen.pt;
           gt.GenJetEta[gt.nGenJet]=gen.eta;
           gt.GenJetPhi[gt.nGenJet]=gen.phi;
           gt.nGenJet++;
         }
         if(gt.nGenJet>=10) break;
     }
     
     gt.genmet=event.GenMET.pt;
     gt.genmetphi = event.GenMET.phi;
}


