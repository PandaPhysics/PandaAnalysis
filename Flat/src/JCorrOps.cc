#include "PandaAnalysis/Utilities/interface/JetCorrector.h"
#include "../interface/JCorrOps.h"
#include "TVector2.h"

using namespace pa;
using namespace std;
using namespace panda;
namespace fj = fastjet;

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

  JetCorrector *jc = new JetCorrector();
  jc->SetYear(analysis.year);
  jc->RunCorrection(analysis.isData,event.rho,&event.muons,&event.chsAK4Jets,&event.rawMet,&event.pfMet,event.runNumber, scale);

  panda::JetCollection *out_jets = 0;
  out_jets = jc->GetCorrectedJets();
  panda::Met *out_met = 0;
  out_met = jc->GetCorrectedMet();

  //event.chsAK4Jets = *out_jets;
  event.pfMet.pt = out_met->pt;
  event.pfMet.phi = out_met->phi;

  METLOOP {
    auto& jets = (*jesShifts)[shift];
    // PF                                                                                                                                                                               
    shiftMET(event.pfMet, jets.vpfMET, i2jes(shift));
    gt.pfmet[shift] = jets.vpfMET.Pt();
    gt.pfmetphi[shift] = jets.vpfMET.Phi();

    // Puppi                                                                                                                                                                            
    shiftMET(event.puppiMet, jets.vpuppiMET, i2jes(shift));
    gt.puppimet[shift] = jets.vpuppiMET.Pt();
    gt.puppimetphi[shift] = jets.vpuppiMET.Phi();

    jets.vpfMETNoMu.SetMagPhi(gt.pfmet[shift], gt.pfmetphi[shift]);
  }
  delete out_jets;
  delete out_met;
}
