#ifndef PANDACORE_TOOLS_JetCorrector
#define PANDACORE_TOOLS_JetCorrector

#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include <map>
#include <string>
#include "TString.h"
#include "PandaTree/Objects/interface/Jet.h"
#include "PandaTree/Objects/interface/Met.h"
#include "PandaTree/Objects/interface/Muon.h"
#include "PandaTree/Objects/interface/CorrT1METJet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "PandaAnalysis/Flat/interface/Operator.h"
/**
 * \brief Corrects a jet collection and optionally propagates to Met 
 *
 */
namespace pa {
	class JetCorrector
	{
	public:
		JetCorrector();
		~JetCorrector();

		void RunCorrection(bool isData, float rho, panda::MuonCollection *muons_, panda::JetCollection *injets_, panda::Met *pfmet_=0,int runNumber = 0, FactorizedJetCorrector *corrector = 0);
                void RunCorrectionLow(bool isData, float rho, panda::MuonCollection *muons_, panda::CorrT1METJetCollection *injets_, panda::Met *pfmet_=0,int runNumber = 0, FactorizedJetCorrector *corrector = 0);
		panda::JetCollection *GetCorrectedJets();
		panda::Met *GetCorrectedMet();
 //               TLorentzVector GetCorrectedMetUnc();
		//		void SetMCCorrector(TString fpath);
		//void SetDataCorrector(TString fpath, TString iov = "all");
		void SetYear(int year);

	private:
		//FactorizedJetCorrector *mMCJetCorrector;	 
		//std::map<TString,FactorizedJetCorrector *> mDataJetCorrectors;	// map from era to corrector

			panda::JetCollection *outjets = 0;
			panda::Met *outmet = 0;
                        panda::Met *lowoutmet = 0;
			EraHandler *era = 0;	
	};
}
#endif
