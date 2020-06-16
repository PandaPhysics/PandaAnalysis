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
//  cout <<"final jetmet=" << met_jer_corr.Pt() << " phi=" << met_jer_corr.Phi() << endl;
  switch (shift) {
  case shiftjes::kNominal:
    pt = temp2.Pt();
    phi = temp2.Phi();
//    cout << "met3=" << pt << endl;
//    cout << "jer_met=" << met_jer_corr.Pt() << endl;
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
        TLorentzVector temp3,temp5;
        temp3.SetPtEtaPhiM(pt,0, jet.phi,0);
        temp5.SetPtEtaPhiM(pt/smearFac,0,jet.phi,0);
        met_jer_corr=met_jer_corr+temp5-temp3;     
        
//        cout << "jetmet=" << (temp5-temp3).Pt() << " phi=" << (temp5-temp3).Phi() << " jetpt=" << pt << endl;
//        cout << "met_jer_corr=" << met_jer_corr.Pt() << " phi=" << met_jer_corr.Phi() << endl;
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
//  cout <<"in function" << endl;
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
    
    //for (auto &j : *ak4Jets) {
    //  jets.all.push_back(shiftJet(j, i2jes(shift), (analysis.hbb || analysis.darkg) && !analysis.isData));
    //}
    
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

  float maxJetEta = analysis.vbf ? 4.7 : 4.7;
  int nJetDPhi = analysis.vbf ? 4 : 5;
  float minMinJetPt = min(cfg.minJetPt, cfg.minBJetPt);
  TLorentzVector antihemJets, antihemchf, antihemnhf,antihemcef, antihemnef;

  
  if(analysis.LowPtJet && analysis.year == 2018 ){
    for(auto &ljet: event.CorrT1METJet){
      if(ljet.eta>-3.0 && ljet.eta<-1.3 && ljet.phi<1.57&&ljet.phi>0.87){
          TLorentzVector temP4;
          temP4.SetPtEtaPhiM(ljet.pt,ljet.eta,ljet.phi,ljet.mass);
          antihemJets += temP4;
      }
    }
  }

  int temmmm=0;
  JESLOOP {
    bool isNominal = (shift == jes2i(shiftjes::kNominal));
    bool metShift = (i2jes(shift) <= shiftjes::kJESTotalDown);
    JESHandler& jets = (*jesShifts)[shift];

    bool leadjet=false;
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

     if (isNominal) {
        if(pt>30 && jet.eta>-3.0 && jet.eta<-1.3 && jet.phi>-1.57 && jet.phi<-0.87) gt.hemveto++;
     }

// only for 2018 anti-hem region

      float CEMF = jet.chEmEF;
      float CHF = jet.chHEF;
      float NEMF = jet.neEmEF;
      float NHF = jet.neHEF;

    if(isNominal && analysis.year == 2018){
      if(pt>100&&aeta<2.4&&NHF<0.8&&CHF>0.1&& ((jet.jetId>>1)%2==1) && !(jet.eta>-1.3&&jet.eta<-3.0&&jet.phi<1.57&&jet.phi>0.87)) gt.hemlead=1;
      if(jet.eta>-3.0 && jet.eta<-1.3 && jet.phi<1.57&&jet.phi>0.87){
          TLorentzVector temP4;
          temP4.SetPtEtaPhiM(jet.pt,jet.eta,jet.phi,jet.mass);
          antihemJets += temP4;
          TLorentzVector temP4chf;
          temP4chf.SetPtEtaPhiM(jet.pt*CHF,jet.eta,jet.phi,jet.mass);
          antihemchf += temP4chf;
          TLorentzVector temP4nhf;
          temP4nhf.SetPtEtaPhiM(jet.pt*NHF,jet.eta,jet.phi,jet.mass);
          antihemnhf += temP4nhf;
          TLorentzVector temP4cef;
          temP4cef.SetPtEtaPhiM(jet.pt*CEMF,jet.eta,jet.phi,jet.mass);
          antihemcef += temP4cef;
          TLorentzVector temP4nef;
          temP4nef.SetPtEtaPhiM(jet.pt*NEMF,jet.eta,jet.phi,jet.mass);
          antihemnef += temP4nef;
      }
    }

      if (aeta > maxJetEta || jw.nominal->maxpt < minMinJetPt)
        continue;
    
       
      if (isNominal) { // perform cleaning only on the nominal jet and save flags
        if (isMatched(matchLeps.get(),0.16,jet.eta,jet.phi))
          {jw.isLep = true;}
        if (!analysis.hbb && isMatched(matchPhos.get(),0.16,jet.eta,jet.phi))
          jw.isPho = true;
        if (analysis.hbb && jet.puId < utils.getCorr(cJetLoosePUID, aeta,min(39.99f,pt)))
          jw.isPileupJet = true;
      }

//      float csv = centralOnly( (analysis.useDeepCSV? jet.btagDeepB : jet.btagCSVV2), aeta);
//      float cmva = centralOnly(jet.btagCMVA, aeta);

      if (jw.nominal->isLep || jw.nominal->isPho || jw.nominal->isPileupJet)
        continue;
     
      bool is_tight = false;

     
      if(isNominal && pt>100 && aeta<2.4) gt.GoodLeadJet=1;
/*
*/
      
      is_tight=((jet.jetId>>1)%2==1); //bit2
//      if ((analysis.vbf || analysis.hbb) && ((!((jet.jetId>>1)==1) && analysis.year == 2016) || (!is_tight && (analysis.year == 2017 || analysis.year == 2018))))
//        continue;


      float csv = centralOnly( (analysis.useDeepCSV? jet.btagDeepB :( analysis.useDeepJet? jet.btagDeepFlavB :jet.btagCSVV2)), aeta);
      float cmva = centralOnly(jet.btagCMVA, aeta);

      if(!is_tight) continue;

      if (pt > cfg.minBJetPt && aeta < 2.4) { // b jets
        if(csv>gt.leadbtag) gt.leadbtag=csv;
        if (csvLoose(csv)) {
          ++(gt.jetNBtags[shift]);
          if (csvMed(csv) && gt.jetNMBtags[shift] <NJET) {
             if (isNominal){
                gt.jetMBpt[gt.jetNMBtags[shift]]=pt;
                gt.jetMBeta[gt.jetNMBtags[shift]]=jet.eta;
                gt.jetMBphi[gt.jetNMBtags[shift]]=jet.phi;
                gt.jetMBflv[gt.jetNMBtags[shift]]=jet.partonFlavour;
             }
            ++(gt.jetNMBtags[shift]);
          }
        }
      
        if (isNominal && gt.jetNMBcand<NJET){
          jw.flavor=jet.partonFlavour;
          jets.bcand.push_back(&jw);
          gt.jetMBcandpt[gt.jetNMBcand]=pt;
          gt.jetMBcandeta[gt.jetNMBcand]=jet.eta;
          gt.jetMBcandphi[gt.jetNMBcand]=jet.phi;
          gt.jetMBcandflv[gt.jetNMBcand]=jet.partonFlavour;
          gt.jetNMBcand++;
       }
     }

      if(isNominal&&pt>max_pt){    // leadjet
         leadjet=false;
         max_pt=pt;
         if(pt>30&&aeta<2.4&&NHF<0.8&&CHF>0.1) {leadjet=true;}
      }      
      

      if (jw.nominal->maxpt > cfg.minJetPt) {
        // for H->bb, don't consider any jet past NJETSAVED, 
        // for other analyses, consider them, just don't save them

        jets.cleaned.push_back(&jw);

        if (jets.cleaned.size() < 3) {
          if (utils.getCorr(cBadECALJets, jet.eta, jet.phi) > 0)
            gt.badECALFilter = 0;
        }

        if (isNominal && analysis.fatjet) {
          isojet->execute();
          if (jw.iso && csvLoose(csv)) {
            ++(gt.isojetNBtags[shift]);
            if (csvMed(csv))
              ++(gt.isojetNMBtags[shift]);
          }
        }
        if (aeta < 2.4) {
          jets.central.push_back(&jw);

          int njet = jets.central.size() - 1;

          if (pt > cfg.minJetPt && gt.nJet[shift]<NJET){
            int nj=gt.nJet[shift];
            gt.jetPt[shift][nj] = pt;
          if (isNominal) {
              jw.central_idx = njet; 
              gt.jetEta[nj] = jet.eta;
              gt.jetPhi[nj] = jet.phi;
              gt.jetCSV[nj] = csv;
            //  gt.jetIsTight[njet] = jet.monojet ? 1 : 0;
              gt.jetIsIso[njet] = jw.iso ? 1 : 0;
          }
         gt.nJet[shift]++;
         
        }

       }
        if(gt.nJetMax<gt.nJet[shift]) gt.nJetMax=gt.nJet[shift];

        int njet = jets.cleaned.size() - 1;
/*
        if (njet < 2 || ((analysis.hbb || analysis.monoh) && njet < cfg.NJETSAVED)) {
          jw.cleaned_idx = njet; 
          gt.jotPt[shift][njet] = pt;
          if (isNominal) {
             double smearFac=1, smearFacUp=1, smearFacDown=1;
             double tempt, temeta,temphi;
             tempt=-1;
             double minR=10;
             for (auto& gjet: event.GenJet){
               temeta=gjet.eta; temphi=gjet.phi;
               double dR= DeltaR2(temeta,temphi,jet.eta,jet.phi);
               if(dR<minR) {minR=dR; tempt=gjet.pt;}
             }
            if (minR<0.15) gt.jotGenPt[njet] = tempt; 
            gt.jotSmear[njet] = jw.pt / jet.pt; // smeared / nominal
            gt.jotEta[njet] = jet.eta;
            gt.jotPhi[njet] = jet.phi;
            gt.jotM[njet] = jet.mass;
            gt.jotCSV[njet] = csv;
            gt.jotFlav[njet] = jw.flavor;
            gt.jotCMVA[njet] = cmva;
    //        gt.jotVBFID[njet] = (aeta < 2.4) ? (jet.monojet ? 1 : 0) : 1;
            gt.jotIso[njet] = jw.iso ? 1 : 0; 

            bjetreg->execute();
          }
        }
*/
     //     gt.nJot[shift]++;

        TLorentzVector vJet; vJet.SetPtEtaPhiM(pt, jet.eta, jet.phi, jet.mass);
        if (metShift && njet < nJetDPhi && pt > cfg.minJetPt) { // only do this for fully-correlated shifts
//          gt.dphipuppimet[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiMET)), (double)gt.dphipuppimet[shift]);
          gt.dphipfmet[shift]    = min(fabs(vJet.DeltaPhi(jets.vpfMET)),    (double)gt.dphipfmet[shift]);
          if (analysis.recoil) {
//            gt.dphipuppiUA[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUA)), (double)gt.dphipuppiUA[shift]);
//            gt.dphipuppiUW[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUW)), (double)gt.dphipuppiUW[shift]);
//            gt.dphipuppiUZ[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUZ)), (double)gt.dphipuppiUZ[shift]);
            gt.dphipfUA[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUA)), (double)gt.dphipfUA[shift]);
            gt.dphipfUW[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUW)), (double)gt.dphipfUW[shift]);
            gt.dphipfUZ[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUZ)), (double)gt.dphipfUZ[shift]);
          }
        }
      }
    }

    if(isNominal){
    if(leadjet) {gt.GoodLeadJet+=2;}
    gt.LeadJetPt=max_pt;}

//    gt.nJotMax = max(gt.nJot[shift], gt.nJotMax);
    if (metShift) {
      switch (gt.whichRecoil) {
        case 0: // MET
//          gt.dphipuppiU[shift] = gt.dphipuppimet[shift];
          gt.dphipfU[shift] = gt.dphipfmet[shift];
          break;
        case -1: // photon
//          gt.dphipuppiU[shift] = gt.dphipuppiUA[shift];
          gt.dphipfU[shift] = gt.dphipfUA[shift];
          break;
        case 1:
//          gt.dphipuppiU[shift] = gt.dphipuppiUW[shift];
          gt.dphipfU[shift] = gt.dphipfUW[shift];
          break;
        case 2:
//          gt.dphipuppiU[shift] = gt.dphipuppiUZ[shift];
          gt.dphipfU[shift] = gt.dphipfUZ[shift];
          break;
        default: // c'est impossible !
          break;
      }
    }

    // dijet system
    if (metShift) {
      jets.sort();
//      vbf->execute();
    }
//    hbb->execute();

    //if (isNominal)
    //  adjet->execute();

    /*
    gt.avgDPhi[shift] = 0;
    std::vector<TLorentzVector> allVecs;
    TLorentzVector vJet(0,0,0,0); 
    for (int m = 0; m<gt.nJot[shift]; m++){
      TLorentzVector tmpJet(0,0,0,0); 
      tmpJet.SetPtEtaPhiM(gt.jotPt[shift][m], gt.jotEta[m], gt.jotPhi[m], gt.jotM[m]);
      vJet +=  tmpJet;
      allVecs.push_back(tmpJet);
      for (int n = m+1; n<gt.nJot[shift]; n++){
	gt.avgDPhi[shift] += 1./gt.nJot[shift]*SignedDeltaPhi(gt.jotPhi[m],gt.jotPhi[n]);
      }
    }
    gt.jotTotM[shift] = vJet.M();
    gt.sphericity[shift] = sphericity(2.,allVecs);
    */
  } // shift loop
  gt.hempt = antihemJets.Pt();
  gt.hemphi= antihemJets.Phi();
  gt.hemchf = antihemchf.Pt()*cos(antihemchf.Phi()-antihemJets.Phi());
  gt.hemnhf = antihemnhf.Pt()*cos(antihemnhf.Phi()-antihemJets.Phi());
  gt.hemcef = antihemcef.Pt()*cos(antihemcef.Phi()-antihemJets.Phi());
  gt.hemnef = antihemnef.Pt()*cos(antihemnef.Phi()-antihemJets.Phi());

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


void JetFlavorOp::partonFlavor(JetWrapper& jw)
{
  auto& jet = jw.get_base();
  jw.flavor=0; jw.genpt=0;
  for (auto* genptr : *genP) {
    auto& gen = pToGRef(genptr);
    int apdgid = abs(gen.pdgId);
    if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
      continue;
    double dr2 = DeltaR2(jet.eta,jet.phi,gen.eta,gen.phi);
    if (dr2<0.09) {
      jw.genpt = gen.pt;
      if (apdgid==4 || apdgid==5) {
        jw.flavor=apdgid;
        return;
      } else {
        jw.flavor=0;
      }
    }
  }
}


void JetFlavorOp::clusteredFlavor(JetWrapper& jw)
{
  auto& jet = jw.get_base();
  for (auto &gen : event.GenJet) {
    if (DeltaR2(gen.eta, gen.phi, jet.eta, jet.phi) < 0.09) {
      int apdgid = abs(gen.partonFlavour);
      jw.genpt = gen.pt;
      if (apdgid == 4 || apdgid == 5) {
        jw.flavor = apdgid;
        return;
      } else {
        jw.flavor = 0;
      }
    }
  }
}


void JetFlavorOp::do_execute()
{
  for (auto& jw : (*jesShifts)[0].all) {
    if (analysis.jetFlavorPartons)
      partonFlavor(jw);
    else
      clusteredFlavor(jw);
    if (jw.cleaned_idx >= 0) {
//      gt.jotFlav[jw.cleaned_idx] = jw.flavor; 
      if (jw.central_idx >= 0) 
        gt.jetFlav[jw.central_idx] = jw.flavor;
    }
  }
}


void IsoJetOp::do_execute()
{
  auto& jw = **currentJet;
  auto& jets = **currentJES;
  auto& jet = jw.get_base();
  float maxIsoEta = analysis.monoh ? 4.5 : 2.5;
  bool isIsoJet = (
        gt.nFatJet == 0 ||
        (fabs(jet.eta) < maxIsoEta &&
         DeltaR2(gt.fjEta[0],gt.fjPhi[0],jet.eta,jet.phi) > cfg.FATJETMATCHDR2)
      );

  jw.iso = isIsoJet;

  if (isIsoJet)
    jets.iso.push_back(&jw);
}


const vector<double> BJetRegOp::Energies::dr2_bins {
  pow(0.05, 2), pow(0.1, 2), pow(0.2, 2), pow(0.3, 2), pow(0.4, 2)
};

/*
void BJetRegOp::do_execute()
{
  auto& jw = **currentJet;
  auto& jet = jw.get_base();
  auto& jets = **currentJES;

  int N = jets.cleaned.size() - 1;

  energies.clear();

  
  TLorentzVector vjet;
  TLorentzVector vraw; 
  vjet.SetPtEtaPhiM(jet.pt,jet.eta,jet.phi,jet.mass);
  vraw.SetPtEtaPhiM(jet.pt*(1-jet.rawFactor),jet.eta,jet.phi,jet.mass); 

  gt.jotEMF[N] = jet.chEmEF + jet.neEmEF;
  gt.jotHF[N] = jet.chHEF + jet.neHEF;
  gt.jotCEF[N] = jet.chEmEF;
  gt.jotNEF[N] = jet.neEmEF;
  gt.jotCHF[N] = jet.chHEF;
  gt.jotNHF[N] = jet.neHEF;
  gt.jotRawPt[N] = vraw.Pt();
  gt.jotRawMt[N] = vraw.Mt();
  gt.jotRawEt[N] = vraw.Et();
  gt.jotRawM[N] = vraw.M();
  gt.jotRawE[N] = vraw.E();
  gt.jotRho[N] = event.fixedGridRhoFastjetAll;
  gt.jotArea[N] = jet.area; 
  energies.jet_e = gt.jotRawE[N];

  // information about 10 hardest and 10 softest PF candidates

}

void VBFSystemOp::do_execute()
{
  auto& jets = **currentJES;

  int shift = jets.shift_idx;

  unsigned idx0=0, idx1=1;
  if (analysis.hbb) {
    if (analysis.fatjet && gt.nFatJet > 0 && event.FatJet[0].pt > 400) {
      const FatJet& fj = event.FatJet[0]; 
      int inc1 = 0;
      if (DeltaR2(jets.cleaned_sorted[idx1]->base->eta, 
                  jets.cleaned_sorted[idx1]->base->phi,
                  fj.eta, fj.phi) < cfg.FATJETMATCHDR2) {
        inc1++;
      }
      if (DeltaR2(jets.cleaned_sorted[idx0]->base->eta, 
                  jets.cleaned_sorted[idx0]->base->phi,
                  fj.eta, fj.phi) < cfg.FATJETMATCHDR2) {
        if (inc1 == 0) {
          inc1 = 1;
          idx0 = idx1; 
        } else {
          idx0 = idx1 + inc1;
          inc1++;
        }
      }
      idx1 += inc1;
    } else if (gt.hbbm[shift] > 0) {
      if (gt.hbbjtidx[shift][0] == 0 || gt.hbbjtidx[shift][1] == 0) {
        idx0++; idx1++;
      }
      if (gt.hbbjtidx[shift][0] == 1 || gt.hbbjtidx[shift][1] == 1) {
        if (idx0 == 0) 
          idx1++;
        else {
          idx0++; idx1++;
        }
      }
    }
  }
  if (jets.cleaned.size() > idx1) {
    TLorentzVector v0 = jets.cleaned_sorted[idx0]->p4();
    TLorentzVector v1 = jets.cleaned_sorted[idx1]->p4();
    gt.jot12Mass[shift] = (v0 + v1).M();
    gt.jot12DPhi[shift] = v0.DeltaPhi(v1);
    gt.jot12DEta[shift] = fabs(v0.Eta() - v1.Eta());
  }
  if (analysis.darkg){
//    auto& theVertex = event.vertices[0];
//    gt.vtxNTrk = theVertex.ntrk;
    gt.vtxScore =event.PV.score;
    gt.vtxChi2 = event.PV.chi2;
  }
}

void HbbSystemOp::do_execute()
{
  auto& jets = **currentJES;
  int shift = jets.shift_idx;

  
  if (shift == jes2i(shiftjes::kNominal)) {
    if (gt.nJot[shift] == 1 && abs(gt.jotEta[0])<2.5){
      (*hbbdJet) = jets.cleaned[0];
      auto& hbbdJetRef = **hbbdJet;
      deepreg->execute();
      gt.jotDeepBReg[0] = hbbdJetRef.breg;
   }
  }	
  

  if (jets.central.size() < 2)
    return;

  btagsorted = jets.central; // copy
  sort(btagsorted.begin(), btagsorted.end(),
       analysis.useCMVA ?
        [](const JetWrapper *x, const JetWrapper *y) { return x->base->btagCMVA > y->base->btagCMVA; } :
        (analysis.useDeepCSV ?
         [](const JetWrapper *x, const JetWrapper *y) { return x->base->btagDeepB > y->base->btagDeepB; } :
         [](const JetWrapper *x, const JetWrapper *y) { return x->base->btagCSVV2 > y->base->btagCSVV2; }
        )
      );
  map<const JetWrapper*, int> order; // needed for output indexing
  for (int i = 0; i != (int)jets.cleaned.size(); ++i)
    order[jets.cleaned[i]] = i;


  array<TLorentzVector,2> hbbd;
  for (int i = 0; i != 2; ++i)  {
    gt.hbbjtidx[shift][i] = order[btagsorted[i]];
    btagsorted[i]->p4(hbbd[i]);
  }

  TLorentzVector hbbsystem = hbbd[0] + hbbd[1];

  gt.hbbpt[shift] = hbbsystem.Pt();
  gt.hbbeta[shift] = hbbsystem.Eta();
  gt.hbbphi[shift] = hbbsystem.Phi();
  gt.hbbm[shift] = hbbsystem.M();

  array<TLorentzVector,2> hbbd_corr, hbbd_dcorr, hbbd_qcorr;
  if (gt.hbbm[shift] > 0) {
    for (int i = 0; i<2; i++) {
      int idx = gt.hbbjtidx[shift][i];
      (*hbbdJet) = jets.cleaned[idx];
      auto& hbbdJetRef = **hbbdJet;

      if (shift == jes2i(shiftjes::kNominal)) {
        deepreg->execute();
        gt.jotDeepBReg[i] = hbbdJetRef.breg;
        gt.jotDeepBRegWidth[i] = hbbdJetRef.bregwidth;
        gt.jotDeepBRegSampled[i] = (event.rng.normal() * hbbdJetRef.bregwidth)  + hbbdJetRef.breg;
      }
      auto scale_fn = [&](float x) { return 1 + (gt.jotPt[shift][idx] / hbbdJetRef.base->pt) * (x - 1); };
      hbbd_dcorr[i].SetPtEtaPhiM(
            gt.jotPt[shift][idx] * scale_fn(gt.jotDeepBReg[i]),
            gt.jotEta[idx],
            gt.jotPhi[idx],
            gt.jotM[idx]
          );
      hbbd_qcorr[i].SetPtEtaPhiM(
            gt.jotPt[shift][idx] * scale_fn(gt.jotDeepBRegSampled[i]),
            gt.jotEta[idx],
            gt.jotPhi[idx],
            gt.jotM[idx]
          );

      // Shifted values for the jet energies to perform the b-jet regression
      if (shift == jes2i(shiftjes::kNominal)) {
        bdtreg->execute();
        gt.jotBReg[i] = hbbdJetRef.breg;
      }
      hbbd_corr[i].SetPtEtaPhiM(
            gt.jotBReg[i] * gt.jotPt[shift][idx],
            gt.jotEta[idx],
            gt.jotPhi[idx],
            gt.jotM[idx]
          );
    }

    TLorentzVector hbbsystem_corr = hbbd_corr[0] + hbbd_corr[1];
    gt.hbbm_reg[shift] = hbbsystem_corr.M();
    gt.hbbpt_reg[shift] = hbbsystem_corr.Pt();

    TLorentzVector hbbsystem_dcorr = hbbd_dcorr[0] + hbbd_dcorr[1];
    gt.hbbm_dreg[shift] = hbbsystem_dcorr.M();
    gt.hbbpt_dreg[shift] = hbbsystem_dcorr.Pt();

    TLorentzVector hbbsystem_qcorr = hbbd_qcorr[0] + hbbd_qcorr[1];
    gt.hbbm_qreg[shift] = hbbsystem_qcorr.M();
    gt.hbbpt_qreg[shift] = hbbsystem_qcorr.Pt();

  } // regression

  if (gt.hbbm > 0) {
    gt.hbbCosThetaJJ[shift] = hbbsystem.CosTheta();
    float csj1;
    if (analysis.bjetBDTReg) {
      if (hbbd_corr[0].Pt() > hbbd_corr[1].Pt())
        csj1 = CosThetaCollinsSoper(hbbd_corr[0], hbbd_corr[1]);
      else
        csj1 = CosThetaCollinsSoper(hbbd_corr[1], hbbd_corr[0]);
    } else {
      if (hbbd[0].Pt() > hbbd[1].Pt())
        csj1 = CosThetaCollinsSoper(hbbd[0], hbbd[1]);
      else
        csj1 = CosThetaCollinsSoper(hbbd[1], hbbd[0]);
    }
    gt.hbbCosThetaCSJ1[shift] = csj1;
  }

  if (gt.hbbm > 0 && gt.nLooseLep > 0 && shift <= jes2i(shiftjes::kJESTotalDown)) {
    TLorentzVector leptonP4, metP4, nuP4, *jet0P4{nullptr}, *jet1P4{nullptr}, WP4, topP4;
    float dRJet0W, dRJet1W;
    bool jet0IsCloser;
    leptonP4.SetPtEtaPhiM((*looseLeps)[0]->pt,(*looseLeps)[0]->eta,(*looseLeps)[0]->phi,(*looseLeps)[0]->mass);

    metP4.SetPtEtaPhiM(gt.pfmet[shift], 0, gt.pfmetphi[shift], 0);
    nuP4 = getNu4Momentum( leptonP4, metP4 );
    WP4 = leptonP4 + nuP4;

    // If using b-jet regression, use the regressed jets for the top mass reconstruction
    // Otherwise, use the un regressed jets
    if (analysis.bjetBDTReg) {
      jet0P4 = &hbbd_corr[0]; jet1P4 = &hbbd_corr[1];
    } else {
      jet0P4 = &hbbd[0]; jet1P4 = &hbbd[1];
    }

    dRJet0W=jet0P4->DeltaR(leptonP4);
    dRJet1W=jet1P4->DeltaR(leptonP4);
    jet0IsCloser = (dRJet0W < dRJet1W);

    topP4 = jet0IsCloser ? (*jet0P4)+WP4 : (*jet1P4)+WP4;
    gt.topMassLep1Met[shift] = topP4.M();
    gt.topWBosonCosThetaCS[shift] = CosThetaCollinsSoper(WP4, jet0IsCloser ? *jet0P4 : *jet1P4);
    gt.topWBosonPt  = WP4.Pt();
    gt.topWBosonEta = WP4.Eta();
    gt.topWBosonPhi = WP4.Phi();
  }
  if (gt.nLooseLep > 1 && gt.hbbm > 0) {
    TLorentzVector HP4;
    if (analysis.bjetBDTReg)
      HP4.SetPtEtaPhiM(gt.hbbpt_reg[0], gt.hbbeta[0], gt.hbbphi[0], gt.hbbm_reg[0]);
    else
      HP4.SetPtEtaPhiM(gt.hbbpt[0], gt.hbbeta[0], gt.hbbphi[0], gt.hbbm[0]);
    TLorentzVector ZHP4 = (*dilep) + HP4;
    TLorentzVector lep1p4, lep2p4;
    lep1p4.SetPtEtaPhiM(looseLeps->at(0)->pt,looseLeps->at(0)->eta,looseLeps->at(0)->phi,looseLeps->at(0)->mass);
    lep2p4.SetPtEtaPhiM(looseLeps->at(1)->pt,looseLeps->at(1)->eta,looseLeps->at(1)->phi,looseLeps->at(1)->mass);
    gt.ZBosonLep1CosThetaStar = CosThetaStar(lep1p4, lep2p4, ZHP4);
  }
}

void AdJetOp::do_execute()
{
  if (gt.hbbjtidx[0][0]<0 || gt.hbbjtidx[0][1]<0)
    return; 
  auto& jets = **currentJES;
  auto* h0 = jets.cleaned[gt.hbbjtidx[0][0]]->base; 
  auto* h1 = jets.cleaned[gt.hbbjtidx[0][1]]->base; 
  for (auto* jw : jets.central) {
    if (jw->base == h0 || jw->base == h1)
     continue; 
    gt.adjetPt = max(gt.adjetPt, jw->pt);
    gt.adjetCMVA = max(gt.adjetCMVA, jw->base->btagCMVA);
  }
}
*/
