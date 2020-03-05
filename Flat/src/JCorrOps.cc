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
    std::pair<double,double> up = METXYCorr_Met_MetPhi(met.ptCorrUp,met.phiCorrUp,METxcorr,METycorr);
    met.ptCorrUp = up.first;
    met.phiCorrUp = up.second;
    std::pair<double,double> down = METXYCorr_Met_MetPhi(met.ptCorrDown,met.phiCorrDown,METxcorr,METycorr);
    met.ptCorrDown = down.first;
    met.phiCorrDown = down.second;
  }
}




void shiftMET(const RecoMet& met, TLorentzVector& v, shiftjes shift)
{
  float pt;
  float phi;
  switch (shift) {
  case shiftjes::kNominal:
    pt = met.pt;
    phi = met.phi;
    break;
  case shiftjes::kJESTotalUp:
    pt = met.ptCorrUp;
    phi = met.phiCorrUp;
    break;
  case shiftjes::kJESTotalDown:
    pt = met.ptCorrDown;
    phi = met.phiCorrDown;
    break;
  default:
    logger.error("shiftMET", "Unknown JES type!");
    exit(1);
  }

  v.SetPtEtaPhiM(pt, 0, phi, 0);
}

void JetCorrOp::do_execute()
{

  gt.pfmetRaw = event.rawMet.pt;
  gt.calomet = event.caloMet.pt;
  gt.sumETRaw = event.pfMet.sumETRaw;
  gt.trkmet = event.trkMet.pt;
  gt.trkmetphi = event.trkMet.phi;
  gt.pfmetsig = event.pfMet.significance;
  gt.puppimetsig = event.puppiMet.significance;

  if (analysis.isData) {
    TString thisEra = utils.eras->getEra(event.runNumber);
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
  if (analysis.rerunJES){
    JetCorrector *jc = new JetCorrector();
    jc->SetYear(analysis.year);
    if (analysis.puppiJets && analysis.puppiMet){
      jc->RunCorrection(analysis.isData,event.rho,&event.muons,&event.puppiAK4Jets,&event.rawMet,&event.puppiMet,event.runNumber,scale,analysis.year,&event.pfCandidates);
      out_jets = jc->GetCorrectedJets();
      out_met = jc->GetCorrectedMet();
      event.puppiMet.pt = out_met->pt;
      event.puppiMet.phi = out_met->phi;
    }
    else if (!analysis.puppiJets && analysis.puppiMet){
      jc->RunCorrection(analysis.isData,event.rho,&event.muons,&event.chsAK4Jets,&event.rawMet,&event.puppiMet,event.runNumber,scale,analysis.year,&event.pfCandidates);
      out_jets = jc->GetCorrectedJets();
      out_met = jc->GetCorrectedMet();
      event.puppiMet.pt = out_met->pt;
      event.puppiMet.phi = out_met->phi;
    }
    else{
      jc->RunCorrection(analysis.isData,event.rho,&event.muons,&event.chsAK4Jets,&event.rawMet,&event.pfMet,event.runNumber,scale,analysis.year,&event.pfCandidates);
      out_jets = jc->GetCorrectedJets();
      out_met = jc->GetCorrectedMet();
      event.pfMet.pt = out_met->pt;
      event.pfMet.phi = out_met->phi;
    }
  }

  if (analysis.puppiMet){
    TLorentzVector puppimet, puppimetUp, puppimetDown;
    puppimet.SetPtEtaPhiM(event.puppiMet.pt,0.0,event.puppiMet.phi,0.0);
    puppimetUp.SetPtEtaPhiM(event.puppiMet.ptCorrUp,0.0,event.puppiMet.phiCorrUp,0.0);
    puppimetDown.SetPtEtaPhiM(event.puppiMet.ptCorrDown,0.0,event.puppiMet.phiCorrDown,0.0);

    event.puppiMet.pt = puppimet.Pt();
    event.puppiMet.phi = puppimet.Phi();
    event.puppiMet.ptCorrUp = puppimetUp.Pt();
    event.puppiMet.phiCorrUp = puppimetUp.Phi();
    event.puppiMet.ptCorrDown = puppimetDown.Pt();
    event.puppiMet.phiCorrDown = puppimetDown.Phi();

  }

  if (!analysis.puppiMet) // xy due to dead tracker regions. PFMET includes CH from all vertices. Not needed for PuppiMET.
    correctMETXY(event.npv,event.runNumber,event.pfMet,analysis.year,analysis.isData);

  METLOOP {
    auto& jets = (*jesShifts)[shift];
    if (!analysis.puppiMet){
    // PF
      shiftMET(event.pfMet, jets.vpfMET, i2jes(shift));
      gt.pfmet[shift] = jets.vpfMET.Pt();
      gt.pfmetphi[shift] = jets.vpfMET.Phi();
    }
    else{
    // Puppi
      shiftMET(event.puppiMet, jets.vpuppiMET, i2jes(shift));
      gt.puppimet[shift] = jets.vpuppiMET.Pt();
      gt.puppimetphi[shift] = jets.vpuppiMET.Phi();
      //if (shift == jes2i(shiftjes::kNominal)){
      //std::cout << "After filling" << std::endl;
      //std::cout << gt.puppimet[shift] << std::endl;
      //}
    }

    jets.vpfMETNoMu.SetMagPhi(gt.pfmet[shift], gt.pfmetphi[shift]);
  }
  delete out_jets;
  delete out_met;
}
