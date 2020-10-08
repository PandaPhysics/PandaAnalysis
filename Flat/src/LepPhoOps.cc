#include "../interface/LepPhoOps.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"

using namespace pa;
using namespace std;
using namespace panda;

void SimpleLeptonOp::scaleFactors()
{
  if (analysis.isData)
    return; 
  // For the hadronic analyses, store a single branch for lepton ID,
  // isolation, and tracking, computed based on the 2 leading loose
  // leptons in the event. Cool guyz get per-leg scalefactors
  // computed in ComplicatedLeptonOp
  for (int iL=0; iL!=TMath::Min(gt.nLooseLep,2); ++iL) {
    auto* lep = (*looseLeps)[iL];
    float pt = lep->pt, eta = lep->eta, aeta = TMath::Abs(eta);
    Muon* mu = dynamic_cast<Muon*>(lep);
  }
}

void SimpleLeptonOp::do_execute() 
{
  // muons
  for (auto& mu : event.Muon) {
    float pt = mu.pt; float eta = mu.eta; float aeta = fabs(eta);
    if (!mu.looseId) 
       continue; // loose ID 
    if (pt<10 || aeta>2.4) 
      continue;
    if (mu.pfRelIso04_all > 0.25)
      continue; // loose iso
    //only keep loose muon
    bool isFake   = mu.tightId  && (mu.pfRelIso03_all < 0.4) && (mu.pfRelIso03_chg < 0.4); // what the hell is this ID?
    bool isMedium = mu.mediumId && (mu.pfRelIso03_all < 0.15);
    bool isTight  = mu.tightId  && (mu.pfRelIso03_all < 0.15)&& pt>20 && aeta<2.4;
    bool isDxyz   = MuonIP(mu.dxy,mu.dz);
    bool isMvaMedium     = (mu.mvaId==2);//mu.mvaMedium;
    bool isMvaTight      = mu.mvaId;//mu.mvaTight;
    bool isMiniIsoMedium = mu.miniIsoId;//mu.miniIsoMedium;
    bool isMiniIsoTight  = mu.miniIsoId;//mu.miniIsoTight;
    if (isTight) gt.nTightMuon++;
    int muSelBit                   = kLoose;
    if (isFake         ) muSelBit |= kFake;
    if (isMedium       ) muSelBit |= kMedium;
    if (isTight        ) muSelBit |= kTight;
    if (isDxyz         ) muSelBit |= kDxyz;
    if (isMvaMedium    ) muSelBit |= kMvaMedium;
    if (isMvaTight     ) muSelBit |= kMvaTight;
    if (isMiniIsoMedium) muSelBit |= kMiniIsoMedium;
    if (isMiniIsoTight ) muSelBit |= kMiniIsoTight;
    int iL=gt.nLooseMuon;
    gt.muonPt[iL]                   = pt;
    gt.muonEta[iL]                  = eta;
    gt.muonPhi[iL]                  = mu.phi;
    gt.muonSelBit[iL]               = muSelBit;
    gt.muonPdgId[iL]                = mu.charge*-13;
    looseLeps->push_back(&mu);
    matchLeps->push_back(&mu); 
    gt.nLooseMuon++;
    if (gt.nLooseMuon>=3)
      break;
  }
 //electrons
  for (auto& ele : event.Electron) {
    float pt = ele.pt; float eta = ele.eta; float aeta = fabs(eta);
//    if (!ele.veto) 
    if(!ele.mvaFall17V2Iso_WPL)  //loose selection
      continue; 
    if (pt<10 || aeta>2.5) 
      continue;
    if (!ElectronIP(ele.eta,ele.dxy,ele.dz)) 
      continue;
    int iL = gt.nLooseElectron;
    bool isFake   = false ;//ele.hltsafe; mark
    bool isMedium = ele.mvaFall17V2Iso_WP90;
    bool isTight  = ele.mvaFall17V2Iso_WP80 && pt>40 && aeta<2.5;
    bool isDxyz   = true; // already selected on this 
    if (isTight) gt.nTightElectron++;
    int eleSelBit            = kLoose;
    if (isFake  ) eleSelBit |= kFake;
    if (isMedium) eleSelBit |= kMedium;
    if (isTight ) eleSelBit |= kTight;
    if (isDxyz  ) eleSelBit |= kDxyz;
    if (ele.mvaFall17V2Iso_WP90) eleSelBit |= kEleMvaWP90;
    if (ele.mvaFall17V2Iso_WP80) eleSelBit |= kEleMvaWP80;
    gt.electronPt[iL]           = pt;
    gt.electronEta[iL]          = eta;
    gt.electronPhi[iL]          = ele.phi;
    gt.electronSelBit[iL]       = eleSelBit;
    gt.electronPdgId[iL]        = ele.charge*-11;
    looseLeps->push_back(&ele);
    matchLeps->push_back(&ele); 
    gt.nLooseElectron++;
    if (gt.nLooseElectron>=3) 
      break;
  }

  gt.nLooseLep = looseLeps->size();
  gt.nTightLep = gt.nTightElectron + gt.nTightMuon;

  if (gt.nLooseLep>0) {
    auto ptsort([](Lepton const* l1, Lepton const* l2)->bool {
      return l1->pt > l2->pt;
});

    int nToSort = TMath::Min(4,gt.nLooseLep);
    std::partial_sort(looseLeps->begin(),looseLeps->begin()+nToSort,looseLeps->end(),ptsort);

    Lepton* lep1 = (*looseLeps)[0];
    METLOOP {
      gt.mT[shift] = MT(lep1->pt,lep1->phi,gt.pfmet[shift],gt.pfmetphi[shift]);
    }
} 

  for (int i = 0; i != min(4, gt.nLooseLep); ++i) {
    Muon *mu = dynamic_cast<Muon*>((*looseLeps)[i]);
    if (mu != nullptr) {
      (*lepPdgId)[i] = mu->charge * -13;
    } else {
      Electron *ele = dynamic_cast<Electron*>((*looseLeps)[i]);
      (*lepPdgId)[i] = ele->charge * -11;
    }
  }

  if (gt.nLooseLep>1 && (*lepPdgId)[0]+(*lepPdgId)[1]==0) {
    TLorentzVector v1,v2;
    Lepton *lep1=(*looseLeps)[0], *lep2=(*looseLeps)[1];
    v1.SetPtEtaPhiM(lep1->pt,lep1->eta,lep1->phi,lep1->mass);
    v2.SetPtEtaPhiM(lep2->pt,lep2->eta,lep2->phi,lep2->mass);
    *dilep = v1+v2;
    gt.diLepMass = dilep->M();
  } else {
    gt.diLepMass = -1;
}

  scaleFactors();
}

void ComplicatedLeptonOp::do_readData(TString dirPath)
{
  rochesterCorrection.reset(new RoccoR(Form("%s/rcdata.2016.v3", dirPath.Data())));
}

void ComplicatedLeptonOp::do_execute()
{
  // muons
  int rocRNGIdx = 0;
  for (auto& mu : event.Muon) {
    float pt = mu.pt; float eta = mu.eta; float aeta = fabs(eta);
    if (!mu.looseId) continue; // loose ID   
    if (pt<2 || aeta>2.4) continue;
    double ptCorrection=1;
    if (analysis.isData) { // perform the rochester correction on the actual particle
      ptCorrection=rochesterCorrection->kScaleDT((int)mu.charge, pt, eta, mu.phi, 0, 0);
    } else if (pt > 0) { // perform the rochester correction to the simulated particle
      // attempt gen-matching to a final state muon
      bool muonIsTruthMatched=false; TLorentzVector genP4; GenPart genParticle;
      for (int iG = 0; iG != (int)event.GenPart.size() && !muonIsTruthMatched; ++iG) {
        genParticle = event.GenPart[iG];
        if (genParticle.status != 1) continue;
        if (genParticle.pdgId != ((int)mu.charge) * -13) continue;
//        genP4.SetPtEtaPhiM(genParticle.pt, genParticle.eta, genParticle.phi, 0.106);
        genP4.SetPtEtaPhiM(genParticle.pt, genParticle.eta, genParticle.phi, 0.106);
        TLorentzVector muP4;
        muP4.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, 0.106);
        double dR = genP4.DeltaR(muP4);
        if (dR < 0.3) muonIsTruthMatched=true;
      } 
      if (muonIsTruthMatched) { // correct using the gen-particle pt
        double random1 = gRandom->Rndm();//event.rng.uniform(rocRNGIdx);
        ptCorrection=rochesterCorrection->kScaleFromGenMC((int)mu.charge, 
                                                          pt, eta, mu.phi, 
                                                          mu.nTrackerLayers, 
                                                          genParticle.pt, 
                                                          random1, 0, 0);
      } else { // if gen match not found, correct the other way
        double random1 = gRandom->Rndm();// event.rng.uniform(rocRNGIdx); 
        double random2 = gRandom->Rndm();//event.rng.uniform(rocRNGIdx);
        ptCorrection=rochesterCorrection->kScaleAndSmearMC((int)mu.charge, 
                                                           pt, eta, mu.phi, 
                                                           mu.nTrackerLayers, 
                                                           random1, random2, 0, 0);
      }
    }
    pt *= ptCorrection;

    mu.pt = pt;
//    float miniRelIso;
    if (analysis.hbb) {
      if (pt<5 || aeta>2.4  || fabs(mu.dxy)>0.5 || fabs(mu.dz)>1.0) 
        continue;
//      miniRelIso = MiniRelIso(mu, pfCandsMap.get(), event.fixedGridRhoFastjetAll);
      if(mu.pfRelIso04_all>0.4 )//&& miniRelIso>0.4) 
        continue;
    } else {
      if (pt<10 || aeta>2.4 ) continue;
      if(mu.pfRelIso04_all>0.25 ) continue;
    }
//    mu.setPtEtaPhiM(pt,eta,mu.phi,0.106);
    matchLeps->push_back(&mu); // Muon jet cleaning and loose cuts are equivalent in HBB land

    bool isFake   = mu.tightId  && mu.pfRelIso03_all < 0.4 && mu.pfRelIso03_chg/mu.pt < 0.4;
    bool isMedium = mu.mediumId && mu.pfRelIso03_all < 0.15;
    bool isTight  = mu.tightId  && mu.pfRelIso04_all < 0.15 && pt>20 && aeta<2.4 ;
    bool isDxyz   = MuonIP(mu.dxy,mu.dz);
    bool isMvaMedium     = (mu.mvaId==2);//mu.mvaMedium;
    bool isMvaTight      = (mu.mvaId==3);//mu.mvaTight;
    bool isMiniIsoMedium = (mu.miniIsoId==2);//mu.miniIsoMedium;
    bool isMiniIsoTight  = (mu.miniIsoId==3);//mu.miniIsoTight;
    if (isTight) gt.nTightMuon++;
    int muSelBit                   = kLoose;
    if (isFake         ) muSelBit |= kFake;
    if (isMedium       ) muSelBit |= kMedium;
    if (isTight        ) muSelBit |= kTight;
    if (isDxyz         ) muSelBit |= kDxyz;
    if (isMvaMedium    ) muSelBit |= kMvaMedium;
    if (isMvaTight     ) muSelBit |= kMvaTight;
    if (isMiniIsoMedium) muSelBit |= kMiniIsoMedium;
    if (isMiniIsoTight ) muSelBit |= kMiniIsoTight;
    int iL=gt.nLooseMuon;
    gt.muonPt[iL]                   = pt;
    gt.muonEta[iL]                  = eta;
    gt.muonPhi[iL]                  = mu.phi;
    gt.muonD0[iL]                   = mu.dxy;
    gt.muonDZ[iL]                   = mu.dz;
    double varXY[2] = {TMath::Abs(mu.eta), mu.pt};
    if(analysis.year == 2017 || analysis.year == 2018) {
     varXY[0] = mu.pt;
     varXY[1] = TMath::Abs(mu.eta);
    }
    gt.muonSfLoose[iL]              = utils.getCorr(cMuLooseID , varXY[0], varXY[1]) *
                                      utils.getCorr(cMuLooseIso, varXY[0], varXY[1]);
    gt.muonSfMedium[iL]             = utils.getCorr(cMuMediumID,  varXY[0], varXY[1]) *
                                      utils.getCorr(cMuMediumIso, varXY[0], varXY[1]);
    gt.muonSfTight[iL]              = utils.getCorr(cMuTightID,  varXY[0], varXY[1]) *
                                      utils.getCorr(cMuTightIso, varXY[0], varXY[1]);
    gt.muonSfUnc[iL]                = utils.getError(cMuMediumID,  varXY[0], varXY[1]) *
                                      utils.getError(cMuMediumIso, varXY[0], varXY[1]);
    gt.muonSfReco[iL] = 1.0;
    gt.muonSelBit[iL]               = muSelBit;
    gt.muonPdgId[iL]                = mu.charge*-13;
    
    looseLeps->push_back(&mu);
    gt.nLooseMuon++;

    gt.muTrigMatch[iL]=0;
//trigger matching
  for(auto &trgobj : event.TrigObj){
      if(trgobj.id==13 && (trgobj.filterBits>>1)%2==1) {
         if(DeltaR2(eta,mu.phi,trgobj.eta,trgobj.phi)<0.09){
           gt.muTrigMatch[iL]=1;
         }
      }
  }

   if (isTight) gt.sf_lepID*=gt.muonSfTight[iL];
   else gt.sf_lepID*=gt.muonSfLoose[iL];

    if (gt.nLooseMuon>=NLEP) 
      break;
  }

//electron
  for (auto& ele : event.Electron) {
    float pt = ele.pt; float eta = ele.eta; float aeta = fabs(eta);
  //  float sc_eta=eta+ele.deltaEtaSC;
    if (analysis.hbb) { 
      if (pt<7 || aeta>2.4 || fabs(ele.dxy)>0.05 || fabs(ele.dz)>0.2)
        continue;
      if(ele.dr03TkSumPt>0.4*pt) 
        continue;
    } else {
      if (pt<10 || aeta>2.5 || ele.cutBased<=0) 
        continue;
    }
    
    if(fabs(eta)<1.479 && !(fabs(ele.dxy)<0.05 && fabs(ele.dz)<0.1)) continue;
    if(fabs(eta)>=1.479 && !(fabs(ele.dxy)<0.1 && fabs(ele.dz)<0.2)) continue; 
     
  //  pt = ele.smearedPt;
 //   ele.setPtEtaPhiM(pt,eta,ele.phi,511e-6);
    matchLeps->push_back(&ele);
    bool isFake   = false;//ele.hltsafe;
    bool isLoose = (ele.cutBased>=2);
    bool isMedium = (ele.cutBased>=3);
    bool isTight  = (ele.cutBased==4) &&pt>40 &&aeta<2.5;
    bool isDxyz = ElectronIP(ele.eta,ele.dxy,ele.dz);
    bool eleMVAPresel = pt > 15 && ele.mvaFall17V2Iso_WPL;
    int eleSelBit            = 0;
    if (isLoose) eleSelBit |=kLoose;
    if (isFake  ) eleSelBit |= kFake;
    if (isMedium) eleSelBit |= kMedium;
    if (isTight ) eleSelBit |= kTight;
    if (isDxyz  ) eleSelBit |= kDxyz;
    if (ele.mvaFall17V2Iso_WP90 && eleMVAPresel) eleSelBit |= kEleMvaWP90;
    if (ele.mvaFall17V2Iso_WP80 && eleMVAPresel) eleSelBit |= kEleMvaWP80;

    if (isTight) gt.nTightElectron++;
    int iL=gt.nLooseElectron;
    gt.electronPt[iL]           = pt;
    gt.electronEta[iL]          = eta;
    gt.electronPhi[iL]          = ele.phi;
    gt.electronD0[iL]           = ele.dxy;
    gt.electronDZ[iL]           = ele.dz;

    
    gt.electronSfLoose[iL]      = utils.getCorr(cEleVeto, eta, pt);   //only for monojet, loose ele use veto id
    gt.electronSfMedium[iL]     = utils.getCorr(cEleMedium, eta, pt);
    gt.electronSfTight[iL]      = utils.getCorr(cEleTight, eta, pt);
    gt.electronSfMvaWP90[iL]    = utils.getCorr(cEleMvaWP90, eta, pt);
    gt.electronSfMvaWP80[iL]    = utils.getCorr(cEleMvaWP80, eta, pt);
    gt.electronSfUnc[iL]        = utils.getError(cEleMedium, eta, pt);
    if(analysis.year == 2017){
       if(pt>20) gt.electronSfReco[iL]       = utils.getCorr(cEleReco, eta, pt);
       else if(pt>10) gt.electronSfReco[iL] = utils.getCorr(cEleLowReco, eta, pt); 
    }
    else gt.electronSfReco[iL]   = utils.getCorr(cEleReco, eta, pt);
    gt.electronSelBit[iL]       = eleSelBit;
    gt.electronPdgId[iL]        = ele.charge*-11;
    gt.electronCombIso[iL] = ele.dr03TkSumPt;
    looseLeps->push_back(&ele);

   gt.eleTrigMatch[iL]=0;
//trigger matching
    for(auto &trgobj : event.TrigObj){
       if(trgobj.id==11 && (trgobj.filterBits>>1)%2==1) {
          if(DeltaR2(eta,ele.phi,trgobj.eta,trgobj.phi)<0.09){
            gt.eleTrigMatch[iL]=1;
          }
       }
    }

   if (isTight) gt.sf_lepID*=gt.electronSfTight[iL];
   else gt.sf_lepID*=gt.electronSfLoose[iL];

    gt.nLooseElectron++; 
    if (gt.nLooseElectron>=NLEP) 
      break;
    }
   gt.nLooseLep = looseLeps->size();

}

void SimplePhotonOp::do_execute()
{
  for (auto& pho : event.Photon) {
//    if (!pho.loose || !pho.csafeVeto)
//      continue;
    float pt = pho.pt;
    if (pt<1) 
      continue;
    float eta = pho.eta, phi = pho.phi;
    if (pt<15 || fabs(eta)>2.5)
      continue;
    
    if (isMatched(matchLeps.get(),0.25,pho.eta,pho.phi))
      continue;

    if(analysis.year==2016) pho.cutBasedBitmap=pho.cutBased;
    if(pho.cutBasedBitmap==0) continue;

    if(!pho.electronVeto) continue;
  
    loosePhos->push_back(&pho);

    gt.nLoosePhoton++;
    if (gt.nLoosePhoton==1) {
      gt.loosePho1Pt = pt;
      gt.loosePho1Eta = eta;
      gt.loosePho1Phi = phi;
      if(pho.electronVeto) gt.loosePhoeveto=1;
      else gt.loosePhoeveto=0;
    }
    if ((pho.cutBasedBitmap&2)==2 && pt>230 && fabs(eta)<1.4442) { // apply eta cut offline
      if (gt.nLoosePhoton==1)
        gt.loosePho1IsTight=1;
      gt.nTightPhoton++;
      tightPhos->push_back(&pho);
    }
  }
  
//  cout << "nPho=" << gt.nLoosePhoton << endl;
  for (auto& pho : event.Photon) {
    float pt = pho.pt;
    if (pt<1) 
      continue; 
    float eta = pho.eta, phi = pho.phi;
    if (pt<230 || fabs(eta)>1.4442)
      continue;

    if (isMatched(matchLeps.get(),0.25,pho.eta,pho.phi))
      continue;

    if(!pho.electronVeto) continue;

    int mask1 = 42;         // 00000000101010
    int mask2 = 10752;      // 10101000000000

    bool id_noiso=((pho.vidNestedWPBitmap&mask1)==mask1);
    bool inv_iso =((pho.vidNestedWPBitmap&mask2)!=mask2); 

    if(!(id_noiso&&inv_iso)) continue;
    gt.nFakePhoton++;
    if (gt.nFakePhoton==1) {
      gt.FakePho1Pt = pt;
      gt.FakePho1Eta = eta;
      gt.FakePho1Phi = phi;
    }
  }

  gt.genPhotonPt = -1;
  gt.genPhotonEta = -99;
  gt.genPhotonPhi = -99;

  TLorentzVector v1,v2;
  int indx=0;
  bool bosonFind=false;
  for (auto& genP : event.GenPart) {
//    if(event.event==18599484) cout << "sampe genID=" << genP.pdgId << " status=" << genP.status << "  pt=" << genP.pt << endl;
    if(genP.pdgId==22 && gt.genPhotonPt<genP.pt && (genP.statusFlags&1)==1 && genP.status==1){
      gt.genPhotonPt=genP.pt;
      gt.genPhotonEta = genP.eta;
      gt.genPhotonPhi = genP.phi;    
    }

    if((fabs(genP.pdgId)==23||(fabs(genP.pdgId)==24)) && genP.status==62 && !bosonFind){
       gt.genBosonPt=genP.pt;
       gt.genBosonEta=genP.eta;
       gt.genBosonPhi=genP.phi;
       gt.genBosonMass=genP.mass;
       bosonFind=true;
    }
    if(fabs(genP.pdgId)>=11 && fabs(genP.pdgId)<=16 && (genP.statusFlags&128)==128){
      if(indx==0) v1.SetPtEtaPhiM(genP.pt,genP.eta,genP.phi,genP.mass);
      if(indx==1) v2.SetPtEtaPhiM(genP.pt,genP.eta,genP.phi,genP.mass);
      indx++;
    }
  }
 
  if(!bosonFind) {
    gt.genBosonPt=(v1+v2).Pt(); 
    gt.genBosonEta=(v1+v2).Eta();
    gt.genBosonPhi=(v1+v2).Phi();
    gt.genBosonMass=(v1+v2).M();
  }
   scaleFactors();  
//  if(!mask) cout << "broken even=" << event.event << " lumi=" << event.luminosityBlock << endl; 
//
  gt.lheminr=10;
  gt.lhedr=0;
  float maxlhepho=-1;
  for (auto& lheP : event.LHEPart) {
    if(lheP.pdgId!=22) continue;
     if(lheP.pt>maxlhepho){
        gt.lhephopt=lheP.pt;
        gt.lhephoeta=lheP.eta;
        gt.lhephophi=lheP.phi;
        maxlhepho=lheP.pt;
     }
     for (auto& genP : event.LHEPart) {
      if(abs(genP.pdgId)<7 || abs(genP.pdgId)==9||abs(genP.pdgId)==21){
        if(DeltaR2(genP.eta, genP.phi, lheP.eta, lheP.phi) < 0.16){ gt.lhedr=1;
      }
        if(sqrt(DeltaR2(genP.eta, genP.phi, lheP.eta, lheP.phi)) < gt.lheminr) gt.lheminr = sqrt(DeltaR2(genP.eta, genP.phi, lheP.eta, lheP.phi));
    }
  }
 }
  
 scaleFactors();
}

void SimplePhotonOp::scaleFactors()
{
  if (analysis.isData)
     return; 
 // prefiring weights
  gt.sf_l1Prefire = 1.0;
  if (analysis.year == 2016 || analysis.year == 2017) {
    for (auto& pho : event.Photon) {
      float pt = pho.pt, eta = pho.eta;
      if (pt<20) 
        continue;
      veryLoosePhos->push_back(&pho);
      gt.sf_l1Prefire *= (1.0 - utils.getCorr(cL1PhotonPreFiring,eta,pt));
    }
  }

  if (gt.nLoosePhoton < 1)
    return;
  float pt = gt.loosePho1Pt, eta = gt.loosePho1Eta;
  if(!analysis.isData)
    gt.sf_pho = utils.getCorr(cPho,eta,pt);
  if (analysis.isData && pt>175) {
    gt.sf_phoPurity = utils.getCorr(cPhoFake, pt);
  }
}


void ComplicatedPhotonOp::do_execute()
{
  for (auto& pho : event.Photon) {
    float pt = pho.pt;
    float eta = pho.eta, phi = pho.phi;
    if (pt<25 || fabs(eta)>2.5)
      continue;
    
    if (isMatched(matchLeps.get(),0.25,pho.eta,pho.phi))
      continue;

    float nhiso_barrel = 1.189 + 0.01512*pho.pt + 2.259e-5*pho.pt*pho.pt;
    float phiso_barrel = 2.08 + 0.004017*pho.pt;
    float chiso_barrel = 1.141;
    float sieie_barrel = 0.01015;
    float hovere_barrel = 0.02197;

    float nhiso_barrel_loose = 24.032 + 0.01512*pho.pt + 2.259e-05*pho.pt*pho.pt;
    float phiso_barrel_loose = 2.876 + 0.004017*pho.pt;

    if (analysis.year == 2016){

      nhiso_barrel = 2.725+0.0148*pho.pt+0.000017*pho.pt*pho.pt;
      phiso_barrel = 2.571+0.0047*pho.pt;
      chiso_barrel = 0.441;
      hovere_barrel = 0.0396;
      sieie_barrel = 0.01022;

      nhiso_barrel_loose = 10.910+0.0148*pho.pt+0.000017*pho.pt*pho.pt;
      phiso_barrel_loose = 3.630+0.0047*pho.pt;

    }
    

    bool pho_medium_barrel = (abs(pho.eta) < 1.479 && pho.sieie < sieie_barrel && pho.hoe < hovere_barrel && pho.pfRelIso03_chg < chiso_barrel && (pho.pfRelIso03_all-pho.pfRelIso03_chg) < nhiso_barrel && pho.pfRelIso03_all < phiso_barrel);
    bool pho_medium_barrel_alter = (abs(pho.eta) < 1.479 && pho.sieie < sieie_barrel && pho.hoe < hovere_barrel && pho.pfRelIso03_chg > chiso_barrel && pho.pfRelIso03_chg < 11. && (pho.pfRelIso03_all-pho.pfRelIso03_chg) < nhiso_barrel_loose && pho.pfRelIso03_all < phiso_barrel_loose);

    float nhiso_endcap = 2.718 + 0.0117*pho.pt + 2.3e-5*pho.pt*pho.pt;
    float phiso_endcap = 3.867 + 0.0037*pho.pt;
    float chiso_endcap = 1.051;
    float sieie_endcap = 0.0272;
    float hovere_endcap = 0.0326;
    
    float nhiso_endcap_loose = 19.722 + 0.0117*pho.pt + 2.3e-05*pho.pt*pho.pt;
    float phiso_endcap_loose = 4.162 + 0.0037*pho.pt;
    
    if (analysis.year == 2016){
      nhiso_endcap = 1.715+0.0163*pho.pt+0.000014*pho.pt*pho.pt;
      phiso_endcap = 3.863+0.0034*pho.pt;
      chiso_endcap = 0.442;
      hovere_endcap = 0.0219;
      sieie_endcap = 0.03001;
      
      nhiso_endcap_loose = 5.931+0.0163*pho.pt+0.000014*pho.pt*pho.pt;
      phiso_endcap_loose = 6.641+0.0034*pho.pt;
    }

    bool pho_medium_endcap = (abs(pho.eta) > 1.479 && pho.sieie < sieie_endcap && pho.hoe < hovere_endcap && pho.pfRelIso03_chg < chiso_endcap && (pho.pfRelIso03_all-pho.pfRelIso03_chg) < nhiso_endcap && pho.pfRelIso03_all < phiso_endcap);
    bool pho_medium_endcap_alter = (abs(pho.eta) > 1.479 && pho.sieie < sieie_endcap && pho.hoe < hovere_endcap && pho.pfRelIso03_chg > chiso_endcap && pho.pfRelIso03_chg < 11. && (pho.pfRelIso03_all-pho.pfRelIso03_chg) < nhiso_endcap_loose && pho.pfRelIso03_all < phiso_endcap_loose);


    if (!(pho_medium_barrel || pho_medium_endcap)){
      if (pho_medium_barrel_alter || pho_medium_endcap_alter){
     /*
	if (gt.alterPho1Pt < pt){
	  gt.alterPho1Pt = pt;
	  gt.alterPho1Eta = eta;
	  gt.alterPho1Phi = phi;
	  gt.alterPho1sieie = pho.sieie;
	  gt.alterPho1r9 = pho.r9;
	  gt.alterPho1hOverE = pho.hoe;
	  gt.alterPho1chIso = pho.pfRelIso03_chg;
	  gt.alterPho1phIso = pho.pfRelIso03_all;
	  gt.alterPho1nhIso = pho.pfRelIso03_all-pho.pfRelIso03_chg;
	  int phoSelBit = 0;
// this is always true as of now, but safer to have it like this

	  if (pho_medium_barrel||pho_medium_endcap) phoSelBit |= pMedium;
	  if (pho.mvaID_WP80)                  phoSelBit |= pTight;
//	  if (pho.highpt)                 phoSelBit |= pHighPt;
//	  if (pho.csafeVeto)              phoSelBit |= pCsafeVeto;
	  if (pho.pixelSeed)              phoSelBit |= pPixelVeto;
//	  if (!pfChargedPhotonMatch(pho)) phoSelBit |= pTrkVeto;
	  if (pho_medium_barrel_alter||pho_medium_endcap_alter) phoSelBit |= pMediumNM1;
	  if (pho.mvaID_WP90)                 phoSelBit |= pMediumPanda; 
	  gt.alterPho1SelBit = phoSelBit;
	}
        */
	continue;
      }
      continue;
    }
    loosePhos->push_back(&pho);
    gt.nLoosePhoton++;
    if (gt.loosePho1Pt < pt) {
      gt.loosePho1Pt = pt;
      gt.loosePho1Eta = eta;
      gt.loosePho1Phi = phi;
      gt.loosePho1sieie = pho.sieie;
      gt.loosePho1r9 = pho.r9;
      gt.loosePho1hOverE = pho.hoe;
      gt.loosePho1chIso = pho.pfRelIso03_chg;
      gt.loosePho1phIso = pho.pfRelIso03_all;
      gt.loosePho1nhIso = pho.pfRelIso03_all-pho.pfRelIso03_chg;
      int phoSelBit = 0;
// this is always true as of now, but safer to have it like this

      if (pho_medium_barrel||pho_medium_endcap) phoSelBit |= pMedium;
      if (pho.mvaID_WP80)                  phoSelBit |= pTight;
//      if (pho.highpt)                 phoSelBit |= pHighPt;
//      if (pho.csafeVeto)              phoSelBit |= pCsafeVeto;
      if (!pho.pixelSeed)              phoSelBit |= pPixelVeto;
//      if (!pfChargedPhotonMatch(pho)) phoSelBit |= pTrkVeto;
      if (pho_medium_barrel_alter||pho_medium_endcap_alter) phoSelBit |= pMediumNM1;
      if (pho.mvaID_WP90)                 phoSelBit |= pMediumPanda; 
      gt.loosePho1SelBit = phoSelBit;
      if (pho.mvaID_WP90 && /*pho.csafeVeto &&*/ !pho.pixelSeed) gt.loosePho1IsTight = 1;
      else                                              gt.loosePho1IsTight = 0;
    }
    if ( (pho_medium_barrel || pho_medium_endcap) && /*pho.csafeVeto &&*/ !pho.pixelSeed) { // apply eta cut offline
      gt.nTightPhoton++;
      tightPhos->push_back(&pho);
//matchPhos->push_back(&pho);
     }
   }
  scaleFactors();
  }
/*
bool ComplicatedPhotonOp::pfChargedPhotonMatch(const Photon& photon)
{
  double matchedRelPt = -1.;

  for (auto& cand : event.Jet) {
    if (cand.charge == 0) 
      continue;

    double dr(cand.dR(photon));
    double rawPt = photon.pt;
    double relPt(cand.pt / rawPt);
    if (dr < 0.1 && relPt > matchedRelPt) {
      matchedRelPt = relPt;
    }

  }

  return (matchedRelPt > 0.6);
}
*/
void TauOp::do_execute()
{
  for (auto& tau : event.Tau) {

    if (analysis.vbf) {
      if (!tau.idDecayMode || !tau.idDecayModeNewDMs)
        continue;
    } else {
      if ((!tau.idDecayMode) || (tau.idMVAoldDM2017v2&2)!=2)
        continue;
    }
    if (tau.pt<18 || fabs(tau.eta)>2.3)
      continue;

    if (isMatched(matchLeps.get(),0.16,tau.eta,tau.phi))
      continue;

    int iL=gt.nTau;
    gt.tauPt[iL]                   = tau.pt;
    gt.tauEta[iL]                  = tau.eta;
    gt.tauPhi[iL]                  = tau.phi;

    gt.nTau++;
    
  }
} 

/*void GenLepOp::do_execute()
{
  gt.genTauPt = -1;
  gt.genElectronPt = -1;
  gt.genMuonPt = -1;
  GenPart *tau = NULL;
  bool foundTauLeptonic = false; 
  for (auto* genptr : *genP) {
    auto& gen = pToGRef(genptr);
    int apdgid = abs(gen.pdgId);
    float pt = gen.pt;
    bool isEmu = false; 

    if (apdgid == 11 && pt > gt.genElectronPt) {
      gt.genElectronPt = pt; 
      gt.genElectronEta = gen.eta; 
      isEmu = true; 
    }
    
    if (apdgid == 13 && pt > gt.genMuonPt) {
      gt.genMuonPt = pt; 
      gt.genMuonEta = gen.eta; 
      isEmu = true; 
    }

    if (isEmu && !foundTauLeptonic && tau) {
      const GenPart *parent = &gen;
      while (parent->genPartIdxMother>=0) {
        if (parent == tau) {
          foundTauLeptonic = true; 
          gt.genTauPt = -1; 
          gt.genTauEta = -1;
          break;
        }
      }
    }

    if (!foundTauLeptonic && apdgid == 15 && pt > gt.genTauPt
        && ((gen.statusFlags>>7)%2==1 
            || (gen.statusFlags >> 11)%2==1
            || ((gen.statusFlags >> 1)%2==1
                && (gen.statusFlags >>8)%2==1
                )
            )
        ) 
    {
      gt.genTauPt = pt; 
      gt.genTauEta = gen.eta;
    }
  }
}
*/
