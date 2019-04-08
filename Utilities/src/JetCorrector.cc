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
		
		//panda::Jet &j_out = outjets->create_back();

		float old_pt = j_in.pt();

		j_in.setPtEtaPhiM(new_pt,j_in.eta(),j_in.phi(),j_in.m());

		met_correction.SetMagPhi(old_pt - new_pt, j_in.phi());

		new_met += met_correction;

	}
	
	outmet->setXY(new_met.X(), new_met.Y());

	// reorder jets by pT
	injets_->sort(panda::Particle::PtGreater);

}

panda::JetCollection *JetCorrector::GetCorrectedJets() { return outjets; outjets = 0; }

panda::Met *JetCorrector::GetCorrectedMet() { return outmet; outmet = 0; }

