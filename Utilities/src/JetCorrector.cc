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

void JetCorrector::RunCorrection(bool isData, float rho, panda::JetCollection *injets_, panda::Met *rawmet_, panda::Met *pfmet_, int runNumber, FactorizedJetCorrector *corrector)
{
        TVector2 new_met {pfmet_->v()};
        TVector2 met_correction {};

	//TLorentzVector v_outmet;
	//v_outmet.SetPtEtaPhiM(rawmet_->pt,0,rawmet_->phi,0);
	outmet = new panda::Met();
	outjets = new panda::JetCollection();

	TLorentzVector v_j_in, v_j_out;
	for (auto &j_in : *injets_) {
		double jecFactor = 1;
		v_j_in.SetPtEtaPhiM(j_in.rawPt,j_in.eta(),j_in.phi(),j_in.m());
		if (fabs(j_in.eta())<5.191) {
			corrector->setJetPt(j_in.rawPt);
			corrector->setJetEta(j_in.eta());
			corrector->setJetPhi(j_in.phi());
			corrector->setJetE(v_j_in.E());
			corrector->setRho(rho);
			corrector->setJetA(j_in.area);
			corrector->setJetEMF(j_in.cef + j_in.nef);
			jecFactor = corrector->getCorrection();
		}
		auto new_pt = jecFactor * j_in.rawPt;

		v_j_out.SetPtEtaPhiM(new_pt,j_in.eta(),j_in.phi(),j_in.m());
		
		panda::Jet &j_out = outjets->create_back();
		
		j_out.setPtEtaPhiM(new_pt,j_in.eta(),j_in.phi(),j_in.m());
		j_out.rawPt = j_in.rawPt;
		j_out.loose = j_in.loose;
		j_out.area = j_in.area;
		j_out.nhf = j_in.nhf;
		j_out.chf = j_in.chf;
		j_out.cef = j_in.cef;
		j_out.nef = j_in.nef;
		j_out.puid = j_in.puid;
		j_out.loose = j_in.loose;
		j_out.tight = j_in.tight;
		j_out.tightLepVeto = j_in.tightLepVeto;
		j_out.monojet = j_in.monojet;
		j_out.matchedGenJet = j_in.matchedGenJet;
		j_out.constituents = j_in.constituents;
		j_out.secondaryVertex = j_in.secondaryVertex;
		j_out.csv = j_in.csv;
		j_out.qgl = j_in.qgl;
		j_out.cmva = j_in.cmva;
		j_out.deepCSVudsg = j_in.deepCSVudsg;
		j_out.deepCSVb = j_in.deepCSVb;
		j_out.deepCSVbb = j_in.deepCSVbb;
		j_out.deepCSVc = j_in.deepCSVc;
		j_out.deepCSVcc = j_in.deepCSVcc;
		j_out.deepCMVAudsg = j_in.deepCMVAudsg;
		j_out.deepCMVAb = j_in.deepCMVAb;
		j_out.deepCMVAbb = j_in.deepCMVAbb;
		j_out.deepCMVAc = j_in.deepCMVAc;
		j_out.deepCMVAcc = j_in.deepCMVAcc;
		j_out.ptCorrUp = j_in.ptCorrUp;
		j_out.ptCorrDown = j_in.ptCorrDown;
		j_out.ptSmear = j_in.ptSmear;
		j_out.ptSmearUp = j_in.ptSmearUp;
		j_out.ptSmearDown = j_in.ptSmearDown;

		met_correction.SetMagPhi(j_in.pt() - new_pt, j_out.phi());		  

		new_met += met_correction;

		  //TLorentzVector v_diff = v_j_in - v_j_out;
		  //v_outmet += v_diff;
	}
	
	//outmet->pt = v_outmet.Pt();
	//outmet->phi = v_outmet.Phi();
	outmet->setXY(new_met.X(), new_met.Y());


	// reorder jets by pT
	outjets->sort(panda::Particle::PtGreater);
}

panda::JetCollection *JetCorrector::GetCorrectedJets() { return outjets; outjets = 0; }

panda::Met *JetCorrector::GetCorrectedMet() { return outmet; outmet = 0; }

