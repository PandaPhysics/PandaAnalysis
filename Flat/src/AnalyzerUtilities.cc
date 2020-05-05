#include "../interface/AnalyzerUtilities.h"
#include <cassert>

namespace fj = fastjet;
using namespace std;
using namespace pa;
using namespace panda;


bool pa::isAncestor(const GenPart& child, const GenPart& ancestor, Event& event)
{
/*  auto* parent = &child;
  while (parent->parent.isValid()) {
    parent = parent->parent.get();
    if (parent == &ancestor)
      return true;
  }
  return false;
*/
  if(child.genPartIdxMother>=0){
    auto& parent = event.GenPart[child.genPartIdxMother];
    if( &parent  == &ancestor)
      return true;
  }
  return false;
}


bool pa::hasChild(const GenPart& parent, const vector<Particle*>& genP, Event& event, bool isHard)
{
  for (auto* pptr : genP) {
    auto& child = pToGRef(pptr);
    if (child.pdgId != parent.pdgId)
      continue;
    if (isHard && !hard(child))
      continue; 
    if (child.genPartIdxMother>=0){
      auto& ancester = event.GenPart[child.genPartIdxMother];
      if( &ancester  == &parent)
      return true;
  
    }
  }
  return false;
}

void pa::downloadData(TString url, TString outpath, bool force, TString opts) 
{
  if (!force && !gSystem->AccessPathName(outpath))
    return; // "Attention, bizarre convention of return value!!" -ROOT docs for this function
  TString cmd = "wget ";
  cmd += opts + " ";
  cmd += "-O " + outpath;
  cmd += " \"" + url + "\"";
  gSystem->Exec(cmd);
}

////////////////////////////////////////////////////////////////////////////////////

JetTree::Node::Node(fj::PseudoJet& pj_):
  _pj(pj_)
{
  fj::PseudoJet dau1, dau2;
  if (_pj.has_parents(dau1, dau2)) {
    l = new Node(dau1);
    r = new Node(dau2);
  }
}

void JetTree::Node::GetTerminals(vector<int>& terminals_)
{
  if (l && r) {
    l->GetTerminals(terminals_);
    r->GetTerminals(terminals_);
  } else {
    terminals_.push_back(_pj.user_index());
  }
}

////////////////////////////////////////////////////////////////////////////////////

ParticleGridder::ParticleGridder(unsigned etaN, unsigned phiN, float etaMax):
  _etaMax(etaMax),
  _phiMax(TMath::Pi() - 0.000001),
  _etaBin(-_etaMax, _etaMax, etaN),
  _phiBin(-_phiMax, _phiMax, phiN)
{
}

void ParticleGridder::clear() 
{
  for (auto &v : _collections) {
    v.second.clear();
  }
  _collections.clear();
  _particles.clear();
  _gridded.clear();
}

void ParticleGridder::add(const panda::Particle& p) 
{
    TLorentzVector temp4;
    temp4.SetPtEtaPhiM(p.pt,p.eta,p.phi,p.mass);

  _particles.push_back(temp4);
}

vector<TLorentzVector>& ParticleGridder::get() 
{
  for (auto &p : _particles) {
    float eta = p.Eta();
    float phi = p.Phi();
    if (fabs(eta) >= _etaMax)
      continue;
    phi = bound(phi, -_phiMax, _phiMax);
    int iEta = _etaBin.find(eta);
    int iPhi = _phiBin.find(phi);
    pair<int,int> i(iEta, iPhi);
    auto iter = _collections.find(i);
    if (iter == _collections.end()) {
      _collections[i] = {&p};
    } else {
      iter->second.push_back(&p);
    }      
  }

  // grid is filled
  _gridded.clear();
  for (auto &iter : _collections) {
    const pair<int,int> &bin = iter.first;
    const vector<TLorentzVector*> coll = iter.second;
    int iEta = bin.first, iPhi = bin.second;
    if (coll.size() > 0) {
      if (_on) {
        TLorentzVector vSum;
        for (auto *p : coll) {
          vSum += *p;
        }
        if (_etaphi) {
          float eta = _etaBin.center(iEta);
          float phi = _phiBin.center(iPhi);
          vSum.SetPtEtaPhiM(vSum.Pt(), eta, phi, vSum.M()); // => no spatial resolution within the cell
        }
//        logger.debug("ParticleGridder out",
//               Form("pt=%.3f,eta=%.3f,phi=%.3f,m=%.3f in %i,%i", vSum.Pt(), eta, phi, vSum.M(), iEta, iPhi));
//        logger.debug("ParticleGridder out",
//               Form("eta = [%.3f,%.3f], phi=[%.3f,%.3f]", 
//                    _etaBin.left(iEta),
//                    _etaBin.left(iEta+1),
//                    _phiBin.left(iPhi),
//                    _phiBin.left(iPhi+1))); 
        _gridded.push_back(vSum);
      } else {
        for (auto *p : coll) {
          _gridded.push_back(*p);
        }
      }
    } else {
      logger.error("ParticleGridder::get",Form("Bin %i,%i is supposed to be non-empty, but found otherwise!",
                                         iEta, iPhi));
    }
  }
  return _gridded;
}

////////////////////////////////////////////////////////////////////////////////////

JetRotation::JetRotation(float x1, float y1, float z1, float x2, float y2, float z2)
{
  TVector3 axis1(x1, y1, z1); // this axis gets rotated onto the z-axis
  TVector3 axis2(x2, y2, z2); // this axis will get rotated into the x-z plane 
  TVector3 axisz(0, 0, 1);
  TVector3 axisx(1, 0, 0);

  r_toz.Rotate(axis1.Angle(axisz), axis1.Cross(axisz)); 
  assert((r_toz*axis1).Angle(axisz) < 0.0001); // allow some rounding

  axis2 = r_toz * axis2;      // first rotate it as before 
  axis2.SetZ(0);              // zero-out the z-component 
  r_inxy.Rotate(axis2.Angle(axisx), axis2.Cross(axisx));
  assert((r_inxy*axis2).Angle(axisx) < 0.0001);
}

void JetRotation::Rotate(float& x, float& y, float& z) 
{
  TVector3 v(x, y, z);
  v = r_toz * v;
  v = r_inxy * v;
  x = v.x(); y = v.y(); z = v.z();
}

////////////////////////////////////////////////////////////////////////////////////

VPseudoJet pa::convertPFCands(const vector<const panda::PFParticle*> &incoll, bool puppi, double minPt) 
{
  VPseudoJet vpj;
  vpj.reserve(incoll.size());
  int idx = -1;
  for (auto *incand : incoll) {
    double factor = 1; //puppi ? incand->puppiW() : 1;
    idx++;
    if (factor*incand->pt<minPt)
      continue;
    TLorentzVector temp4;
    temp4.SetPtEtaPhiM(incand->pt,incand->eta,incand->phi,incand->mass);
    vpj.emplace_back(factor*temp4.Px(),factor*temp4.Py(),
                     factor*temp4.Pz(),factor*temp4.E());
    vpj.back().set_user_index(idx);
  }
  return vpj;
}

VPseudoJet pa::convertPFCands(const panda::RefVector<panda::PFParticle> &incoll, bool puppi, double minPt) 
{
  vector<const panda::PFParticle*> outcoll;
  outcoll.reserve(incoll.size());
  for (auto incand : incoll)
    outcoll.push_back(incand.get());

  return pa::convertPFCands(outcoll, puppi, minPt);
}


////////////////////////////////////////////////////////////////////////////////////

bool pa::ElectronIP(double eta, double dxy, double dz) 
{
  double aeta = fabs(eta);
  if (aeta<1.4442) {
    return (dxy < 0.05 && dz < 0.10) ;
  } else {
    return (dxy < 0.10 && dz < 0.20);
  }
}

bool pa::MuonIP(double dxy, double dz) 
{
  return (dxy < 0.02 && dz < 0.10);
}


