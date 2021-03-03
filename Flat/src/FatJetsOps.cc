#include "../interface/FatJetsOps.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"
#include "fastjet/contrib/Njettiness.hh"

using namespace pa;
using namespace std;
using namespace panda;
using namespace fastjet;

JetWrapper BaseJetOp::shiftJet(const FatJet& jet, shiftjes shift, bool smear)
{
  float pt = jet.pt;
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
      if (minR<0.15){
        jer->getStochasticSmear(pt,jet.eta,event.fixedGridRhoFastjetAll,smearFac,smearFacUp,smearFacDown,1,tempt, analysis.year);
      }
      else {
        jer->getStochasticSmear(pt,jet.eta,event.fixedGridRhoFastjetAll,smearFac,smearFacUp,smearFacDown,0,-99, analysis.year);
      }
      pt *= smearFac;
  }
  if (shift != shiftjes::kNominal) {
    int ishift = jes2i(shift);
    bool isUp = !(ishift % 2 == 0);
      (*scaleUnc)[ishift]->setJetPt(pt);
      (*scaleUnc)[ishift]->setJetEta(jet.eta);
      double relShift = (*scaleUnc)[ishift]->getUncertainty(isUp);
      if (!isUp)
        relShift = -relShift;
      pt *= (1 + relShift);
  }
  return JetWrapper(pt, jet);
}


void FatJetOp::setupJES()
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
        return;
      }
    }
  } else {
    scaleUnc = &(scaleUncs["MC"]);
  }
}

void FatJetOp::do_execute()
{

  vector<fastjet::PseudoJet> finalTracks;
  vector<TLorentzVector> tracksP4;

  gt.HTTot = 0;

  for (auto &pfcand : event.PFCands){
    //std::cout << pfcand.pt << std::endl;                                                                                                                                                
    if (pfcand.trkPt>1.){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(pfcand.trkPt,pfcand.trkEta,pfcand.trkPhi,pfcand.mass);
      tracksP4.push_back(tmp);
      fastjet::PseudoJet curtrk(tmp.Px(),tmp.Py(),tmp.Pz(),tmp.E());
      if (!std::isinf(curtrk.perp())){
	finalTracks.emplace_back(curtrk);
	gt.HTTot += tmp.Pt();
      }
    }
  }

  fastjet::JetDefinition *jetDef = new fastjet::JetDefinition(fastjet::kt_algorithm,1.5);
  fastjet::ClusterSequence seq(sorted_by_pt(finalTracks), *jetDef);
  vector<fastjet::PseudoJet> allJets(sorted_by_pt(seq.inclusive_jets(200)));

  gt.nFatJet=0;
  gt.nFatJetTrunc = 0;

  int max_mult = 0;
  int idx_mult = 0;

  for (auto& trackJet : allJets) {
    gt.nFatJet++;

    if ((int)trackJet.constituents().size()>max_mult){
      max_mult = trackJet.constituents().size();
      idx_mult = gt.nFatJet-1;
    }

    if (gt.nFatJetTrunc < nMaxFJ) {
      gt.nFatJetTrunc++;
    }
    else
      continue;

    gt.fjEta[gt.nFatJetTrunc-1] = trackJet.eta();
    gt.fjPhi[gt.nFatJetTrunc-1] = trackJet.phi();
    gt.fjNconst[gt.nFatJetTrunc-1] = trackJet.constituents().size();
    
    JESLOOP {
      gt.fjPt[shift][gt.nFatJetTrunc-1] = trackJet.perp();
      gt.fjM[shift][gt.nFatJetTrunc-1] = trackJet.m();
    }
  }  

  // Three algos for three different SUEP candidates
  // SUEP candidates

  if (gt.nFatJet>0){

  // highest multiplicity

    JESLOOP {
      gt.SUEP_mult_pt[shift] = allJets[idx_mult].perp();
      gt.SUEP_mult_m[shift] = allJets[idx_mult].m();
    }
    
    gt.SUEP_mult_eta = allJets[idx_mult].eta();
    gt.SUEP_mult_phi = allJets[idx_mult].phi();
    gt.SUEP_mult_nconst = allJets[idx_mult].constituents().size();
    
    TLorentzVector SUEP_mult; SUEP_mult.SetPtEtaPhiM(allJets[idx_mult].perp(),allJets[idx_mult].eta(),allJets[idx_mult].phi(),allJets[idx_mult].m());
    TVector3 SUEP_mult_boost = SUEP_mult.BoostVector();
    vector<fastjet::PseudoJet> constituents_mult = allJets[idx_mult].constituents();
    vector<TLorentzVector> tracksP4_mult;
    for (unsigned int j=0; j<constituents_mult.size(); j++){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(constituents_mult[j].perp(),constituents_mult[j].eta(),constituents_mult[j].phi(),constituents_mult[j].m());
      tmp.Boost(-SUEP_mult_boost);
      tracksP4_mult.push_back(tmp);
    }
    
    gt.SUEP_mult_spher = sphericity(2.,tracksP4_mult);

    // highest pT

    JESLOOP {
      gt.SUEP_pt_pt[shift] = gt.fjPt[shift][0];
      gt.SUEP_pt_m[shift] = gt.fjM[shift][0];
    }

    gt.SUEP_pt_eta = gt.fjEta[0];
    gt.SUEP_pt_phi = gt.fjPhi[0];
    gt.SUEP_pt_nconst = gt.fjNconst[0];  

    TLorentzVector SUEP_pt; SUEP_pt.SetPtEtaPhiM(allJets[0].perp(),allJets[0].eta(),allJets[0].phi(),allJets[0].m());
    TVector3 SUEP_pt_boost = SUEP_pt.BoostVector();
    vector<fastjet::PseudoJet> constituents_pt = allJets[0].constituents();
    vector<TLorentzVector> tracksP4_pt;
    for (unsigned int j=0; j<constituents_pt.size(); j++){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(constituents_pt[j].perp(),constituents_pt[j].eta(),constituents_pt[j].phi(),constituents_pt[j].m());
      tmp.Boost(-SUEP_pt_boost);
      tracksP4_pt.push_back(tmp);
    }
    gt.SUEP_pt_spher = sphericity(2.,tracksP4_pt);
    //std::cout << "Spher " <<  gt.SUEP_pt_spher << std::endl;

  }

  // ISR removal (remove N hardest tracks)

  TLorentzVector SUEP_isr(0.,0.,0.,0.);
  int N = 10;

  if ((int)finalTracks.size()>N){ // sanity cut

    vector<fastjet::PseudoJet> sorted_tracks = sorted_by_pt(finalTracks);
    sorted_tracks.erase(sorted_tracks.begin(),sorted_tracks.begin()+N);    

    for (auto &trk : sorted_tracks){
      TLorentzVector tmp; tmp.SetPtEtaPhiM(trk.perp(),trk.eta(),trk.phi(),trk.m());
      SUEP_isr += tmp;
    } 
    
    gt.SUEP_isr_pt = SUEP_isr.Pt();
    gt.SUEP_isr_eta = SUEP_isr.Eta();
    gt.SUEP_isr_phi = SUEP_isr.Phi();
    gt.SUEP_isr_m = SUEP_isr.M();
    gt.SUEP_isr_nconst = (int)finalTracks.size()-N;
    
    TVector3 SUEP_isr_boost = SUEP_isr.BoostVector();
    vector<TLorentzVector> tracksP4_isr;
    for (auto &trk : sorted_tracks){
      TLorentzVector tmp; tmp.SetPtEtaPhiM(trk.perp(),trk.eta(),trk.phi(),trk.m());
      tmp.Boost(-SUEP_isr_boost);
      tracksP4_isr.push_back(tmp);
    }
    gt.SUEP_isr_spher = sphericity(2.,tracksP4_isr);
  }
}

const GenPart* FatJetMatchingOp::matchGen(double eta, double phi, double radius, int pdgid) const
{
  const GenPart* found=NULL;
  double r2 = radius*radius;
  pdgid = abs(pdgid);

  for (auto iG=genObjects.begin();
      iG!=genObjects.end(); ++iG) {
    if (found!=NULL)
      break;
    if (pdgid!=0 && abs(iG->first->pdgId)!=pdgid)
      continue;
    if (DeltaR2(eta,phi,iG->first->eta,iG->first->phi)<r2)
      found = iG->first;
  }

  return found;
}


void FatJetMatchingOp::do_execute()
{
//  if (fjPtrs->size() == 0)
//    return; 

  int pdgidTarget=0;
  if (!analysis.isData && analysis.processType>=kTT && analysis.processType<=kSignal) {
    switch(analysis.processType) {
      case kTop:
      case kTT:
      case kSignal:
        pdgidTarget=6;
        break;
      case kV:
        pdgidTarget=24;
        break;
      case kH:
        pdgidTarget=25;
        break;
      default:
        // analysis.processType>=kTT means we should never get here
        logger.error("FatJetMatchingOp::do_execute","Reached an unknown process type");
    }

    std::vector<int> targets;

    int nGen = genP->size();
    for (int iG=0; iG!=nGen; ++iG) {
      auto& part = pToGRef((*genP)[iG]);
      int pdgid = part.pdgId;
      int abspdgid = abs(pdgid);
      if (abspdgid == pdgidTarget)
        targets.push_back(iG);
    } //looking for targets

    for (int iG : targets) {
      auto& part = pToGRef((*genP)[iG]);

      // check there is no further copy:
      bool isLastCopy=true;
      for (int jG : targets) {
        auto& temG=event.GenPart[pToGPtr((*genP)[jG])->genPartIdxMother];
        if (&temG == &part) {
          isLastCopy=false;
          break;
        }
      }
      if (!isLastCopy)
        continue;

      // (a) check it is a hadronic decay and if so, (b) calculate the size
      if (analysis.processType==kTop||analysis.processType==kTT) {

        // first look for a W whose parent is the top at iG, or a W further down the chain
        const GenPart* lastW(0);
        for (int jG=0; jG!=nGen; ++jG) {
          const GenPart& partW = pToGRef((*genP)[jG]);
          if (TMath::Abs(partW.pdgId)==24 && partW.pdgId*part.pdgId>0) {
            // it's a W and has the same sign as the top
            if (!lastW && &event.GenPart[partW.genPartIdxMother] == &part) {
              lastW = &partW;
            } else if (lastW && &event.GenPart[partW.genPartIdxMother] == lastW) {
              lastW = &partW;
            }
          }
        } // looking for W
        if (!lastW) {// ???
          continue;
        }
        auto& partW(*lastW);

        // now look for b or W->qq
        int iB=-1, iQ1=-1, iQ2=-1;
        double size=0, sizeW=0;
        for (int jG=0; jG!=nGen; ++jG) {
          auto& partQ = pToGRef((*genP)[jG]);
          int pdgidQ = partQ.pdgId;
          int abspdgidQ = TMath::Abs(pdgidQ);
          if (abspdgidQ>5)
            continue;
          if (abspdgidQ==5 && iB<0 && &event.GenPart[partQ.genPartIdxMother] == &part) {
            // only keep first copy
            iB = jG;
            size = TMath::Max(DeltaR2(part.eta,part.phi,partQ.eta,partQ.phi),size);
          } else if (abspdgidQ<5 && &event.GenPart[partQ.genPartIdxMother] == &partW) {
            if (iQ1<0) {
              iQ1 = jG;
              size = TMath::Max(DeltaR2(part.eta,part.phi,partQ.eta,partQ.phi),
                  size);
              sizeW = TMath::Max(DeltaR2(partW.eta,partW.phi,partQ.eta,partQ.phi),
                  sizeW);
            } else if (iQ2<0) {
              iQ2 = jG;
              size = TMath::Max(DeltaR2(part.eta,part.phi,partQ.eta,partQ.phi),
                  size);
              sizeW = TMath::Max(DeltaR2(partW.eta,partW.phi,partQ.eta,partQ.phi),
                  sizeW);
            }
          }
          if (iB>=0 && iQ1>=0 && iQ2>=0)
            break;
        } // looking for quarks


        bool isHadronic = (iB>=0 && iQ1>=0 && iQ2>=0); // all 3 quarks were found
        if (isHadronic)
          genObjects[&part] = size;

        bool isHadronicW = (iQ1>=0 && iQ2>=0);
        if (isHadronicW)
          genObjects[&partW] = sizeW;

      } else { // these are W,Z,H - 2 prong decays

        int iQ1=-1, iQ2=-1;
        double size=0;
        for (int jG=0; jG!=nGen; ++jG) {
          auto& partQ = pToGRef((*genP)[jG]);
          int pdgidQ = partQ.pdgId;
          int abspdgidQ = TMath::Abs(pdgidQ);
          if (abspdgidQ>5)
            continue;
          if (&event.GenPart[partQ.genPartIdxMother] == &part) {
            if (iQ1<0) {
              iQ1=jG;
              size = TMath::Max(DeltaR2(part.eta,part.phi,partQ.eta,partQ.phi),
                  size);
            } else if (iQ2<0) {
              iQ2=jG;
              size = TMath::Max(DeltaR2(part.eta,part.phi,partQ.eta,partQ.phi),
                  size);
            }
          }
          if (iQ1>=0 && iQ2>=0)
            break;
        } // looking for quarks

        bool isHadronic = (iQ1>=0 && iQ2>=0); // both quarks were found

        // add to collection
        if (isHadronic)
          genObjects[&part] = size;
      }

    } // loop over targets
  } // process is interesting

  int iFJ = -1; 
  for (auto& fj : event.FatJet) {
    ++iFJ; 
    // first see if jet is matched
    auto* matched = matchGen(fj.eta,fj.phi,1.5,pdgidTarget);
    if (matched!=nullptr) {
      gt.fjIsMatched[iFJ] = 1;
      gt.fjGenPt[iFJ] = matched->pt;
      gt.fjGenSize[iFJ] = genObjects[matched];
    } else {
      gt.fjIsMatched[iFJ] = 0;
    }
    if (pdgidTarget==6) { // matched to top; try for W
      auto* matchedW = matchGen(fj.eta,fj.phi,1.5,24);
      if (matchedW!=nullptr) {
        gt.fjIsWMatched[iFJ] = 1;
        gt.fjGenWPt[iFJ] = matchedW->pt;
        gt.fjGenWSize[iFJ] = genObjects[matchedW];
      } else {
        gt.fjIsWMatched[iFJ] = 0;
      }
    }

    bool found_b_from_g=false;
    int bs_inside_cone=0;
    int has_gluon_splitting=0;
    const GenPart* first_b_mo(0);
    // now get the highest pT gen particle inside the jet cone
    for (auto* genptr : *genP) {
      auto& gen = pToGRef(genptr);
      float pt = gen.pt;
      int pdgid = gen.pdgId;
      if (pt>(gt.fjHighestPtGenPt[iFJ])
          && DeltaR2(gen.eta,gen.phi,fj.eta,fj.phi)<cfg.FATJETMATCHDR2) {
        gt.fjHighestPtGenPt[iFJ] = pt;
        gt.fjHighestPtGen[iFJ] = pdgid;
      }

      if (gen.genPartIdxMother>=0 &&  event.GenPart[gen.genPartIdxMother].pdgId==gen.pdgId)
        continue;

      //count bs and cs
      int apdgid = abs(pdgid);
      if (apdgid!=5 && apdgid!=4)
        continue;

      if (DeltaR2(gen.eta,gen.phi,fj.eta,fj.phi)<cfg.FATJETMATCHDR2) {
        gt.fjNHF[iFJ]++;
        if (apdgid==5) {
          if (gen.genPartIdxMother>=0 && event.GenPart[gen.genPartIdxMother].pdgId==21 && event.GenPart[gen.genPartIdxMother].pt>20) {
            if (!found_b_from_g) {
              found_b_from_g=true;
              first_b_mo=&event.GenPart[gen.genPartIdxMother];
              bs_inside_cone+=1;
            } else if (&event.GenPart[gen.genPartIdxMother]==first_b_mo) {
              bs_inside_cone+=1;
              has_gluon_splitting=1;
            } else {
              bs_inside_cone+=1;
            }
          } else {
            bs_inside_cone+=1;
          }
        }
      }
    }

    gt.fjNbs[iFJ]=bs_inside_cone;
    gt.fjgbb[iFJ]=has_gluon_splitting;

    if (analysis.btagSFs && iFJ == 0) {
      // now get the subjet btag SFs
      vector<btagcand> sj_btagcands;
      vector<double> sj_sf_cent, sj_sf_bUp, sj_sf_bDown, sj_sf_mUp, sj_sf_mDown;
      int nSJ = (fj.subJetIdx1>=0)+(fj.subJetIdx2>=0);
      for (int iSJ=0; iSJ!=nSJ; ++iSJ) {
        auto& subjet = (iSJ==0) ? event.SubJet[fj.subJetIdx1] : event.SubJet[fj.subJetIdx2];
        int flavor=0;
        for (auto* genptr : *genP) {
          auto& gen = pToGRef(genptr);
          int apdgid = abs(gen.pdgId);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(subjet.eta,subjet.phi,gen.eta,gen.phi);
          if (dr2<0.09) {
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the subjet flavor

        float pt = subjet.pt;
        float btagUncFactor = 1;
        float eta = subjet.eta;
        double eff(1),sf(1),sfUp(1),sfDown(1);
        if (flavor==5) {
          eff = utils.getCorr(cCSVBL, pt, fabs(eta));
        } else if (flavor==4) {
          eff = utils.getCorr(cCSVCL, pt, fabs(eta));
        } else {
          eff = utils.getCorr(cCSVLL, pt, fabs(eta));
        }
        if (analysis.hbb)
          utils.btag->calcSF(bSubJetM,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
        else
          utils.btag->calcSF(bSubJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
        sj_btagcands.push_back(btagcand(iSJ,flavor,eff,sf,sfUp,sfDown));
        sj_sf_cent.push_back(sf);
        if (flavor>0) {
          sj_sf_bUp.push_back(sfUp); sj_sf_bDown.push_back(sfDown);
          sj_sf_mUp.push_back(sf); sj_sf_mDown.push_back(sf);
        } else {
          sj_sf_bUp.push_back(sf); sj_sf_bDown.push_back(sf);
          sj_sf_mUp.push_back(sfUp); sj_sf_mDown.push_back(sfDown);
        }

      } // loop over subjets
      utils.btag->evalSF(sj_btagcands,sj_sf_cent,GeneralTree::bCent,GeneralTree::bSubJet,true);
      utils.btag->evalSF(sj_btagcands,sj_sf_bUp,GeneralTree::bBUp,GeneralTree::bSubJet,true);
      utils.btag->evalSF(sj_btagcands,sj_sf_bDown,GeneralTree::bBDown,GeneralTree::bSubJet,true);
      utils.btag->evalSF(sj_btagcands,sj_sf_mUp,GeneralTree::bMUp,GeneralTree::bSubJet,true);
      utils.btag->evalSF(sj_btagcands,sj_sf_mDown,GeneralTree::bMDown,GeneralTree::bSubJet,true);
    }

  }
}

