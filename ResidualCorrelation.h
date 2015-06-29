//*********************************************
// Residual correlation class
//
//*********************************************

#ifndef ResidualCor_H
#define ResidualCor_H

#include LednickyEqn.h

class ResidualCorrelation{
 public:
  ResidualCorrelation();
  virtual ~ResidualCorrelation();
  TH2D *GetTransformMatrix();
  TGraph *GetLednickyGraph();
  void SetLednickyParameters();
  
 private:
  void SetTransformMatrix(/* */);

  TH2D *fTransform;
  LednickyEqn fLedEqn;

};
