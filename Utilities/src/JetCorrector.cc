#include "../interface/JetCorrector.h"
#include "TLorentzVector.h"
#define mPI 3.14159265

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

void JetCorrector::RunCorrection(bool isData, float rho, panda::MuonCollection *muons_, panda::JetCollection *injets_,  panda::Met *pfmet_, int runNumber, FactorizedJetCorrector *corrector)
{
        TVector2 new_met;
        new_met.SetMagPhi(pfmet_->pt,pfmet_->phi);// {pfmet_->v()};
        TVector2 met_correction {};

	outmet = new panda::Met();
        
	outjets = new panda::JetCollection();

	TLorentzVector v_j_in, v_j_in_corr;
	for (auto &j_in : *injets_) {
//          if(j_in.pt<0) continue;
	  double jecFactor = 1;
	  v_j_in.SetPtEtaPhiM(j_in.pt*(1-j_in.rawFactor),j_in.eta,j_in.phi,j_in.mass);
	  v_j_in_corr.SetPtEtaPhiM(j_in.pt,j_in.eta,j_in.phi,j_in.mass);
//          v_j_in.SetPtEtaPhiM(j_in.pt*(1-j_in.rawFactor)*(1-j_in.muonSubtrFactor),j_in.eta,j_in.phi,0);
//          double mu_pt = j_in.pt*(1-j_in.rawFactor)*j_in.muonSubtrFactor;

	  if (fabs(j_in.eta)<5.191) {
 	    corrector->setJetEta(v_j_in.Eta());
 	    corrector->setJetPt (v_j_in.Pt());
 	    corrector->setJetPhi(v_j_in.Phi());
 	    corrector->setJetE  (v_j_in.E());
	    corrector->setRho(rho);
	    corrector->setJetA(j_in.area);
	    corrector->setJetEMF(j_in.chEmEF + j_in.neEmEF);
	    jecFactor = corrector->getCorrection();
	  }
	  auto new_pt_ini = jecFactor * j_in.pt*(1-j_in.rawFactor);
//          auto new_pt_ini = jecFactor * j_in.pt*(1-j_in.rawFactor)*(1-j_in.muonSubtrFactor) + mu_pt;
 	  float old_pt_ini = j_in.pt;
 
// 	  j_in.setPtEtaPhiM(new_pt_ini,j_in.eta,j_in.phi,j_in.mass);
//          cout <<"debug" << endl;
//          cout << "og jet pt=" << j_in.pt << " JES jet pt=" << new_pt_ini << " raw jet pt=" << j_in.pt*(1-j_in.rawFactor) << " eta=" << j_in.eta << endl;
//          cout << "EMF=" << j_in.chEmEF + j_in.neEmEF << " jec=" << jecFactor << endl;
          j_in.pt=new_pt_ini;         // where to assign the JES 

          if(TMath::Abs(old_pt_ini - new_pt_ini) < 0.01) continue;
 
          if(j_in.chEmEF + j_in.neEmEF > 0.9) continue;

/*
          bool isMuonOverlap = false;
          for (auto& mu : *muons_) {
            if (mu.pt>0 && (mu.isGlobal) &&
               DeltaR2(mu.eta,mu.phi,v_j_in.Eta(),v_j_in.Phi()) < 0.16) {
               TLorentzVector muV; muV.SetPtEtaPhiM(mu.pt,mu.eta,mu.phi,mu.mass);
 	      v_j_in -= muV;
 	      v_j_in_corr -= muV;
 	      isMuonOverlap = true;
//              cout << "overlap mu pt=" << mu.pt << endl;
 	    }
          }

 	  if (fabs(v_j_in.Eta())<5.191 && isMuonOverlap == true) {
 	    corrector->setJetPt (v_j_in.Pt());
 	    corrector->setJetEta(v_j_in.Eta());
 	    corrector->setJetPhi(v_j_in.Phi());
 	    corrector->setJetE  (v_j_in.E());
 	    corrector->setRho(rho);
 	    corrector->setJetA(j_in.area);
 	    corrector->setJetEMF(j_in.chEmEF + j_in.neEmEF);
 	    jecFactor = corrector->getCorrection();
 	  }
 	  else if (isMuonOverlap == true) jecFactor = 1;

          
 	  auto new_pt = jecFactor * v_j_in.Pt();
 	  float old_pt = v_j_in_corr.Pt();
*/
  //        cout << "factor=" << jecFactor << " jet in pt=" << v_j_in.Pt() << " muon frac=" << j_in.muonSubtrFactor << endl;
 //         cout << "new pt=" << new_pt << " old_pt=" << old_pt << endl;

          if(new_pt_ini<=15) continue; 
//          if(new_pt <= 15) continue;

          float old_pt = old_pt_ini*(1-j_in.muonSubtrFactor);
          float new_pt = new_pt_ini*(1-j_in.muonSubtrFactor);
          if(old_pt - new_pt>0) met_correction.SetMagPhi(old_pt - new_pt, j_in.phi);
          else met_correction.SetMagPhi(new_pt-old_pt, j_in.phi+mPI);  //unfortunate it can't convert negative pt value
           
	  new_met += met_correction;

/*
         int ishift = jes2i(shiftjes::kJESTotalUp);
         bool isUp = true;
  
         (*scaleUnc)[ishift]->setJetPt(j_in.pt);
         (*scaleUnc)[ishift]->setJetEta(j_in.eta);
         double relShift = (*scaleUnc)[ishift]->getUncertainty(isUp);
    
         TLorentzVector temp2;
         temp2.SetPtEtaPhiM(-relShift*j_in.pt,0, j_in.phi,0);
         met_corr+=temp2;
*/
	}
	
        outmet->pt=sqrt(new_met.X()*new_met.X()+new_met.Y()*new_met.Y());
        outmet->phi=new_met.Phi();
	// reorder jets by pT
}

panda::JetCollection *JetCorrector::GetCorrectedJets() { return outjets; outjets = 0; }

panda::Met *JetCorrector::GetCorrectedMet() { return outmet; outmet = 0; }

void JetCorrector::RunCorrectionLow(bool isData, float rho, panda::MuonCollection *muons_, panda::CorrT1METJetCollection *injets_, panda::Met *pfmet_, int runNumber, FactorizedJetCorrector *corrector)
{
        TVector2 met_correction {};

        TLorentzVector v_j_in, v_j_in_corr;

        for (auto &j_in : *injets_) {
          if(j_in.pt<0) continue;
          double jecFactor = 1;
          v_j_in.SetPtEtaPhiM(j_in.rawPt,j_in.eta,j_in.phi,0);
          double mu_pt = j_in.rawPt*j_in.muonSubtrFactor;

          if (fabs(j_in.eta)<5.191) {
            corrector->setJetEta(v_j_in.Eta());
            corrector->setJetPt (v_j_in.Pt());
            corrector->setJetPhi(v_j_in.Phi());
            corrector->setJetE  (v_j_in.E());
            corrector->setRho(rho);
            corrector->setJetA(j_in.area);
            corrector->setJetEMF(0);
            jecFactor = corrector->getCorrection();
          }
          auto new_pt_ini = jecFactor * v_j_in.Pt()+mu_pt;

          j_in.pt=new_pt_ini;         // where to assign the JES 

        }
}


//TLorentzVector JetCorrector::GetCorrectedMetUnc() { return met_corr;}
