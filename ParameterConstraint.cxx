//********************************
//
//********************************

#include "ParameterConstraint.h"

ParameterConstraint::ParameterConstraint(Int_t constrainedParam, vector<Int_t> constrainedSystems, const vector<PairSystem*> &pairSystems)
{
  fParameter = constrainedParam;
  fConstrainedSystems = constrainedSystems;
  SortConstraints(pairSystems);
}

ParameterConstraint::~ParameterConstraint()
{

}

void ParameterConstraint::SortConstraints(const vector<PairSystem*> &pairSystems)
{
  // Sort the constraints so that they follow the same order as the pair systems

  // Fill a vector with system types in the order that they appear
  // in the pair system vector
  vector<Int_t> systemOrder;
  for(UInt_t iSys = 0; iSys < pairSystems.size(); iSys++)
  {
    systemOrder.push_back(pairSystems[iSys]->GetSystemType());
  }

  vector<Int_t> tempOrder(fConstrainedSystems);
  Int_t entry = 0;
  // Put the systems into fConstrainedSystem in the same order
  // as they appear in systemorder
  for(UInt_t iSys = 0; iSys < systemOrder.size(); iSys++)
  {
    for(UInt_t iTemp = 0; iTemp < tempOrder.size(); iTemp++)
    {
      if(systemOrder[iSys] == tempOrder[iTemp])
      {
	fConstrainedSystems[entry] = tempOrder[iTemp];
	entry++;
      }
    }    
  }
}
