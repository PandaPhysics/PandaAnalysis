#include "../interface/JetCorrector.h"
#include "TLorentzVector.h"

using namespace pa; 
using namespace std;

JetCorrector::JetCorrector() 
{ 
}

JetCorrector::~JetCorrector()
{
	delete era;
}

void JetCorrector::SetYear(int year){
  era = new EraHandler(year);
}

void JetCorrector::RunCorrection(bool isData, float rho, panda::MuonCollection *muons_, panda::JetCollection *injets_, panda::Met *rawmet_, panda::Met *pfmet_, int runNumber, FactorizedJetCorrector *corrector, int year, panda::PFCandCollection *cands_)
{
        TVector2 new_met {pfmet_->v()};
        TVector2 met_correction {};

	TLorentzVector vBadJets;
	vBadJets.SetPtEtaPhiM(0,0,0,0);

	//TLorentzVector v_outmet;
	//v_outmet.SetPtEtaPhiM(rawmet_->pt,0,rawmet_->phi,0);
	outmet = new panda::Met();
	outjets = new panda::JetCollection();

	TLorentzVector v_j_in, v_j_in_corr;
	for (auto &j_in : *injets_) {
	  double jecFactor = 1;
	  v_j_in.SetPtEtaPhiM(j_in.rawPt,j_in.eta(),j_in.phi(),j_in.m());
	  v_j_in_corr.SetPtEtaPhiM(j_in.pt(),j_in.eta(),j_in.phi(),j_in.m());
	  if (fabs(j_in.eta())<5.191) {
 	    corrector->setJetEta(v_j_in.Eta());
 	    corrector->setJetPt (v_j_in.Pt());
 	    corrector->setJetPhi(v_j_in.Phi());
 	    corrector->setJetE  (v_j_in.E());
	    corrector->setRho(rho);
	    corrector->setJetA(j_in.area);
	    corrector->setJetEMF(j_in.cef + j_in.nef);
	    jecFactor = corrector->getCorrection();
	  }
	  auto new_pt_ini = jecFactor * j_in.rawPt;

 	  float old_pt_ini = j_in.pt();
 
 	  j_in.setPtEtaPhiM(new_pt_ini,j_in.eta(),j_in.phi(),j_in.m());

	  // 2017 EE fix  
	  if (year == 2017){
	    if(abs(j_in.eta()) > 2.650 && abs(j_in.eta()) < 3.139 && j_in.rawPt < 50 && j_in.pt() > 15 && (j_in.cef + j_in.nef) <= 0.9){
	    //if(abs(j_in.eta()) > 2.650 && abs(j_in.eta()) < 3.139 && j_in.rawPt < 50 && j_in.pt() > 15){
	      TLorentzVector vj;
	      //std::cout << j_in.pt() << std::endl; 
	      vj.SetPtEtaPhiM(j_in.pt(),j_in.eta(),j_in.phi(),j_in.m());
	      //vj.SetPtEtaPhiM(v_j_in.Pt(),j_in.eta(),j_in.phi(),v_j_in.M());
	      vBadJets = vBadJets + vj;
	    }
	  }

          if(TMath::Abs(old_pt_ini - new_pt_ini) < 0.01) continue;
 
          if(j_in.cef + j_in.nef > 0.9) continue;
 
          bool isMuonOverlap = false;
          for (auto& mu : *muons_) {
            if (mu.pt()>0 && (mu.global || mu.standalone) &&
               DeltaR2(mu.eta(),mu.phi(),v_j_in.Eta(),v_j_in.Phi()) < 0.16) {
               TLorentzVector muV; muV.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(),mu.m());
 	      v_j_in -= muV;
 	      v_j_in_corr -= muV;
 	      isMuonOverlap = true;
 	    }
          }
 	  if (fabs(v_j_in.Eta())<5.191 && isMuonOverlap == true) {
 	    corrector->setJetPt (v_j_in.Pt());
 	    corrector->setJetEta(v_j_in.Eta());
 	    corrector->setJetPhi(v_j_in.Phi());
 	    corrector->setJetE  (v_j_in.E());
 	    corrector->setRho(rho);
 	    corrector->setJetA(j_in.area);
 	    corrector->setJetEMF(j_in.cef + j_in.nef);
 	    jecFactor = corrector->getCorrection();
 	  }
 	  else if (isMuonOverlap == true) jecFactor = 1;
 
 	  auto new_pt = jecFactor * v_j_in.Pt();
 	  float old_pt = v_j_in_corr.Pt();

          if(new_pt <= 15) continue;

	  if (old_pt - new_pt >= 0)
	    met_correction.SetMagPhi(old_pt - new_pt, j_in.phi());
	  else
	    met_correction.SetMagPhi(new_pt - old_pt, 3.1415+j_in.phi());

	  new_met += met_correction;
	}

	// EEfix for 2017 
	if (year == 2017){
	  std::vector<int> clustered_idxs;
	  for (auto &j_in : *injets_) {
	    for (const auto& pf : j_in.constituents) {
	      clustered_idxs.push_back(pf.idx());
	    }
	  }

	  int counter = 0;
	  TLorentzVector vBadPFs;
	  vBadPFs.SetPtEtaPhiM(0,0,0,0);
	  for (auto& cand : *cands_) {
	    if (abs(cand.eta()) > 2.650 && abs(cand.eta()) < 3.139){
	      if (std::find(clustered_idxs.begin(),clustered_idxs.end(),counter) != clustered_idxs.end())
		continue;
	      else{
		TLorentzVector tmppf;
		tmppf.SetPtEtaPhiM(cand.pt(),cand.eta(),cand.phi(),cand.m());
		vBadPFs += tmppf;
	      }
	    }
	    counter++;
	  }
	  
	  TVector2 BADJETS,BADPFS;
	  BADJETS.SetMagPhi(vBadJets.Pt(), vBadJets.Phi());
	  BADPFS.SetMagPhi(vBadPFs.Pt(), vBadPFs.Phi());
	  TVector2 NEWMET;
	  NEWMET = BADPFS+BADJETS+new_met;	
	  outmet->setXY(NEWMET.X(), NEWMET.Y());
	  //outmet->setXY(new_met.X(), new_met.Y());
	}
	else
	  outmet->setXY(new_met.X(), new_met.Y());


	// reorder jets by pT
	//outjets->sort(panda::Particle::PtGreater);
	injets_->sort(panda::Particle::PtGreater);
}

panda::JetCollection *JetCorrector::GetCorrectedJets() { return outjets; outjets = 0; }

panda::Met *JetCorrector::GetCorrectedMet() { return outmet; outmet = 0; }

