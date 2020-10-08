#include "PandaAnalysis/Utilities/interface/JetCorrector.h"
#include "../interface/JCorrOps.h"
#include "TVector2.h"

using namespace pa;
using namespace std;
using namespace panda;
namespace fj = fastjet;


std::pair<double,double> METXYCorr_Met_MetPhi(double pt, double phi, double METxcorr, double METycorr){
  
  std::pair<double,double>  TheXYCorr_Met_MetPhi(pt,phi);
  double CorrectedMET_x = pt*cos(phi)+METxcorr;
  double CorrectedMET_y = pt*sin(phi)+METycorr;
  
  double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  double CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;
  
  TheXYCorr_Met_MetPhi.first= CorrectedMET;
  TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;
  return TheXYCorr_Met_MetPhi;

}


void correctMETXY(int npv, int runnb, RecoMet& met, int year, bool isData)
{
  // XY correction
  double METxcorr(0.),METycorr(0.);
  if (year == 2018){
    if (isData){
      
	//Type 1 PFMET
      if(runnb >=315252 &&runnb<=316995) METxcorr = -(0.362865*npv -1.94505);
      if(runnb >=315252 &&runnb<=316995) METycorr = -(0.0709085*npv -0.307365);
      if(runnb >=316998 &&runnb<=319312) METxcorr = -(0.492083*npv -2.93552);
      if(runnb >=316998 &&runnb<=319312) METycorr = -(0.17874*npv -0.786844);
      if(runnb >=319313 &&runnb<=320393) METxcorr = -(0.521349*npv -1.44544);
      if(runnb >=319313 &&runnb<=320393) METycorr = -(0.118956*npv -1.96434);
      if(runnb >=320394 &&runnb<=325273) METxcorr = -(0.531151*npv -1.37568);
      if(runnb >=320394 &&runnb<=325273) METycorr = -(0.0884639*npv -1.57089);
      
      // Type 2
      /*
      if(runnb >=315252 &&runnb<=316995) METxcorr = -(0.362642*npv +-1.55094);
      if(runnb >=315252 &&runnb<=316995) METycorr = -(0.0737842*npv +-0.677209);
      if(runnb >=316998 &&runnb<=319312) METxcorr = -(0.485614*npv +-2.45706);
      if(runnb >=316998 &&runnb<=319312) METycorr = -(0.181619*npv +-1.00636);
      if(runnb >=319313 &&runnb<=320393) METxcorr = -(0.503638*npv +-1.01281);
      if(runnb >=319313 &&runnb<=320393) METycorr = -(0.147811*npv +-1.48941);
      if(runnb >=320394 &&runnb<=325273) METxcorr = -(0.520265*npv +-1.20322);
      if(runnb >=320394 &&runnb<=325273) METycorr = -(0.143919*npv +-0.979328);
      */
    }
    else{
      
      //Type 1 PFMET
      METxcorr = -(0.296713*npv -0.141506);
      METycorr = -(0.115685*npv +0.0128193);
      
      // Type 2
      //METxcorr = -(0.299448*npv +-0.13866);
      //METycorr = -(0.118785*npv +0.0889588);
    }

    std::pair<double,double> nominal = METXYCorr_Met_MetPhi(met.pt,met.phi,METxcorr,METycorr);
    met.pt = nominal.first;
    met.phi = nominal.second;
//    std::pair<double,double> up = METXYCorr_Met_MetPhi(met.ptCorrUp,met.phiCorrUp,METxcorr,METycorr);
//    met.ptCorrUp = up.first;
//    met.phiCorrUp = up.second;
//    std::pair<double,double> down = METXYCorr_Met_MetPhi(met.ptCorrDown,met.phiCorrDown,METxcorr,METycorr);
//    met.ptCorrDown = down.first;
//    met.phiCorrDown = down.second;
  }
}
/*
void JetCorrOp::lowPtClean(){
  for(auto& ljet : event.CorrT1METJet){
    bool clean=false;
    for (auto& jet: event.Jet){
     if(DeltaR(ljet.eta,ljet.phi,jet.eta,jet.phi)<0.15) {clean=true; break;}
    }
    if(!clean) continue;
    ljet.pt=ljet.rawPt;
    ljet.rawFactor=0;
    ljet.mass=0;
    ljet.eta=100;
  } 
}
*/
void JetCorrOp::do_execute()
{
/*
  if(analysis.LowPtJet){
     lowPtClean();
  }
*/
  gt.pfmetRaw = event.RawMET.pt;
  gt.calomet = event.CaloMET.pt;
  gt.sumETRaw = event.RawMET.sumEt;
  gt.trkmet = event.TkMET.pt;
  gt.trkmetphi = event.TkMET.phi;
  //gt.pfmetsig = event.MET.significance;
//  gt.puppimetsig = event.PuppiMET.significance;
//  cout << "met1=" << event.MET.pt << endl;
  if (analysis.isData) {
    TString thisEra = utils.eras->getEra(event.run);
    for (auto& iter : scaleUncs) {
      if (!iter.first.Contains("data")){
        continue;
      }
      if (iter.first.Contains(thisEra)) {
        scaleUnc = &(scaleUncs[iter.first]);
        scale = scales[iter.first].get();
        break;
      }
    }
  } else {
    scaleUnc = &(scaleUncs["MC"]);
    scale = scales["MC"].get();
  }

  panda::JetCollection *out_jets = 0;
  panda::Met *out_met = 0;

//  TLorentzVector met_corr;
  if (analysis.rerunJES){
    JetCorrector *jc = new JetCorrector();
    jc->SetYear(analysis.year);
    if ( analysis.puppiMet){
      jc->RunCorrection(analysis.isData,event.fixedGridRhoFastjetAll,&event.Muon,&event.Jet,&event.PuppiMET,event.run, scale);
      out_jets = jc->GetCorrectedJets();
      out_met = jc->GetCorrectedMet();
      event.PuppiMET.pt = out_met->pt;
      event.PuppiMET.phi = out_met->phi;
    }
    else if( analysis.year == 2017 ){
      jc->RunCorrection(analysis.isData,event.fixedGridRhoFastjetAll,&event.Muon,&event.Jet,&event.METFixEE2017,event.run, scale);
      out_jets = jc->GetCorrectedJets();
      out_met = jc->GetCorrectedMet();
        event.METFixEE2017.pt = out_met->pt;
        event.METFixEE2017.phi = out_met->phi;
        gt.pfmet_jes = out_met->pt;
        gt.pfmetjesphi = out_met->phi;

       
      if(analysis.LowPtJet){
        JetCorrector *jc2 = new JetCorrector();
        jc2->RunCorrectionLow(analysis.isData,event.fixedGridRhoFastjetAll,&event.Muon,&event.CorrT1METJet,&event.METFixEE2017,event.run, scale);
      }
    }
    else{
      jc->RunCorrection(analysis.isData,event.fixedGridRhoFastjetAll,&event.Muon,&event.Jet,&event.MET,event.run, scale);
      out_jets = jc->GetCorrectedJets();
      out_met = jc->GetCorrectedMet();

        event.MET.pt = out_met->pt;
        event.MET.phi = out_met->phi;
        gt.pfmet_jes = out_met->pt;
        gt.pfmetjesphi = out_met->phi;
      if(analysis.LowPtJet){
        JetCorrector *jc2 = new JetCorrector();
        jc2->RunCorrectionLow(analysis.isData,event.fixedGridRhoFastjetAll,&event.Muon,&event.CorrT1METJet,&event.MET,event.run, scale);
      }

    }
  }

  if (analysis.puppiMet){
    // Puppi Hack: taking photons out of puppimet
    TLorentzVector puppimet, puppimetUp, puppimetDown;
    puppimet.SetPtEtaPhiM(event.PuppiMET.pt,0.0,event.PuppiMET.phi,0.0);
//    puppimetUp.SetPtEtaPhiM(event.PuppiMET.pt,0.0,event.PuppiMET.phi,0.0);
//    puppimetDown.SetPtEtaPhiM(event.PuppiMET.pt,0.0,event.PuppiMET.phi,0.0);
//    puppimetUp.SetPtEtaPhiM(event.PuppiMET.ptCorrUp,0.0,event.puppiMet.phiCorrUp,0.0);
//    puppimetDown.SetPtEtaPhiM(event.puppiMet.ptCorrDown,0.0,event.puppiMet.phiCorrDown,0.0);
    /*
    for (auto& cand : event.pfCandidates) {
      if (abs(cand.pdgId()) == 22 && cand.pt()>80){
	TLorentzVector pho;
	pho.SetPtEtaPhiM((1-cand.puppiW())*cand.pt(),cand.eta(),cand.phi(),0);
	puppimet -= pho;
	puppimetUp -= pho;
	puppimetDown -= pho;
      }
    }
    */
    //std::cout << "Before photon corr" << std::endl;
    //std::cout << event.puppiMet.pt << std::endl;

    event.PuppiMET.pt = puppimet.Pt();
    event.PuppiMET.phi = puppimet.Phi();
/*
    event.PuppiMET.ptCorrUp = puppimetUp.Pt();
    event.PuppiMET.phiCorrUp = puppimetUp.Phi();
    event.PuppiMET.ptCorrDown = puppimetDown.Pt();
    event.PuppiMET.phiCorrDown = puppimetDown.Phi();
*/
    //std::cout << "After photon corr" << std::endl;
    //std::cout << event.puppiMet.pt << std::endl;

  }

//  if (!analysis.puppiMet) // xy due to dead tracker regions. PFMET includes CH from all vertices. Not needed for PuppiMET.
//    correctMETXY(event.PV.npvs,event.run,event.MET,analysis.year,analysis.isData);
  delete out_jets;
  delete out_met;
}
