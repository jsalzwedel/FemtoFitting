//********************************
//
//********************************

#include "ParameterConstraint.h"

ParameterConstraint::ParameterConstraint(Int_t constrainedParam, vector<Int_t> constrainedSystems)
{
  fParameter = constrainedParam;
  fConstrainedSystems = constrainedSystems;
}

ParameterConstraint::~ParameterConstraint()
{

}
