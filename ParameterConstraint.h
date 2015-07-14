//********************************
//
//********************************

#ifndef ParameterContraint_H
#define ParameterContraint_H

class ParameterConstraint{
 public:
  ParameterConstraint(Int_t constrainedParam, vector<Int_t> constrainedSystems);
  virtual ~ParameterConstraint();
  Int_t GetConstrainedParam() {return fParameter;};
  vector<Int_t> GetConstrainedSystems() {return fConstrainedSystems};
  

 private:
  Int_t fParameter; // The index of the parameter that is constrained
  vector<Int_t> fConstrainedSystems; // Indices of the systems that share the constraint
};
