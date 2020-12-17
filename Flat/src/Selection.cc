#include "../interface/Selection.h"
#include "../interface/Common.h"

using namespace pa; 

bool VBFGamma::do_accept() const
{
  return false;
}


bool SUEPSel::do_accept() const
{
  if((gt->HLT_PFHT1050==1 || gt->HLT_PFJet500) && gt->nFatJet>0)
    return true;
  else
    return false;
}

