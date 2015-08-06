//********************************
//
//********************************

#ifndef ParameterContraint_H
#define ParameterContraint_H

#include <vector>
#include "Rtypes.h"
#include "PairSystem.h"

using std::vector;

class ParameterConstraint{
 public:
  ParameterConstraint(Int_t constrainedParam, vector<Int_t> constrainedSystems, const vector<PairSystem*> &pairSystems);
  virtual ~ParameterConstraint();
  Int_t GetConstrainedParam() {return fParameter;};
  vector<Int_t> GetConstrainedSystems() {return fConstrainedSystems;};
  

 private:
  void SortConstraints(const vector<PairSystem*> &pairSystems);
  Int_t fParameter; // The index of the parameter that is constrained
  vector<Int_t> fConstrainedSystems; // Indices of the systems that share the constraint
};

#endif
