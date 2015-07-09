//************************************************
// Class to contain all the fitting info for a 
// given pair type / centrality combination.
// One of these will be created by the Fitter
// for each pair type / centrality being fit.
//************************************************

#ifndef PairSystem_H
#define PairSystem_H

#include "LednickyEqn.h"

class PairSystem{
 public:
  PairSystem(TH1D *cfData/*, */);
  ~PairSystem();
  Double_t CalculateFitChisquare(vector<Double_t> pars);
  void SetCFData(TH2D* correlationFunction);
  void CreateNewLednickyEqn(Bool_t isIdentical, TH2D *transformMatrix);

 private:
  TString fPairTypeName; // e.g. "LambdaLambda" or "LambdaAntilambda"
  TString fCentralityName; // e.g. "0-10%"
  vector<*LednickyEqn> fLednickyEqns; // List of primary and secondary L&L eqn objects
  TH1D *fCF; // Correlation function data for the pair type
  //... fit range (Double? Bin integer?)
  
  void SetLednickyParams(vector<Double_t> pars);
  
  

};
