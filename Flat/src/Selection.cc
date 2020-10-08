#include "../interface/Selection.h"
#include "../interface/Common.h"

using namespace pa; 

bool VBFGamma::do_accept() const
{

//    if(gt->nTightMuon>=1 || gt->nTightElectron>=1 ||gt->nTightPhoton==1 || gt->nFakePhoton==1)
//      if((gt->nTightMuon>=1 && gt->nLooseMuon==2) || (gt->nTightElectron>=1 && gt->nLooseElectron))
//       if(gt->genBosonPt>=100|| gt->pfUZmag[0]>170)
//      if(gt->genPhotonPt>=130 || gt->lhephopt>=130)
     if((gt->whichRecoil==2 && gt->pfUZmag[0]>200) || (gt->whichRecoil==1 && gt->pfUWmag[0]>200) || (gt->whichRecoil==-1 && gt->pfUAmag[0]>200) || (gt->whichRecoil==0 && gt->pfmet[0]>200))
//      if( gt->hemlead>0  && gt->hempt>80 && gt->trigger>=16)
      return true;

  return false;

}


bool SUEPSel::do_accept() const
{
  return true;
}

