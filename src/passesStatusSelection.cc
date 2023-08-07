#include "TauAnalysis/Entanglement/interface/passesStatusSelection.h"

#include <algorithm> // std::find()

bool
passesStatusSelection(int status, const vint& selection)
{
  if ( std::find(selection.begin(), selection.end(), status) != selection.end() ) return true;
  else                                                                            return false;
}
