#include "../interface/Selection.h"
#include "../interface/Common.h"

using namespace pa; 

bool VBFGamma::do_accept() const
{

     if(gt->nJetMax>=2)
      return true;

  return false;

}

