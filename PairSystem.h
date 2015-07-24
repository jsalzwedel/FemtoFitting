//************************************************
// Class to contain all the fitting info for a 
// given pair type / centrality combination.
// One of these will be created by the Fitter
// for each pair type / centrality being fit.
//************************************************

#ifndef PairSystem_H
#define PairSystem_H

#include "TH1D.h"
#include "TGraph.h"

#include <vector>

using std::vector;

class LednickyEqn;
/* class TH1D; */

class PairSystem{
 public:
  PairSystem(TH1D *cfData, const vector<LednickyInfo> &ledInfo, TString pairTypeName, Int_t sysIndex);
  ~PairSystem();
  Double_t CalculateFitChisquare();
  Int_t GetSystemIndex() const {return fSystemIndex;};
  void SetCFData(TH2D* correlationFunction);
  void SetLednickyParameters(vector<Double_t> pars);
  void SetLowFitBin(Int_t bin) {fLowFitBin = bin;};
  void SetHighFitBin(Int_t bin) {fHighFitBin = bin;};
  void SetSystemIndex(Int_t sysIndex) {fSystemIndex = sysIndex;};
  TGraph *GetCombinedTGraph();
  TH1D *GetCF() {return fCF;};

 private:
  /* void CreateDefaultLednickyEqns(); */
  /* void CreateNewLednickyEqn(TString name, Bool_t isIdentical, TH2D *transformMatrix); */
  /* void ReadInLambdaParams(); */
  TString fPairTypeName; // e.g. "LambdaLambda" or "LambdaAntilambda"
  TString fCentralityName; // e.g. "0-10"
  vector<LednickyEqn*> fLednickyEqns; // List of primary and secondary L&L eqn objects
  vector<Double_t> fLambdaParameters; // Lambda parameters associated with the L&L eqn objects
  TH1D *fCF; // Correlation function data for the pair type
  //... fit range (Double? Bin integer?)
  Double_t fNorm; // Normalization parameter for combined Lednicky graph
  Bool_t fUseEstimatedLambdaParams;
  Int_t fNBins; // Number of kstar bins in correlation function
  Double_t fBinWidth; // Width of each kstar bin
  Int_t fLowFitBin; // Lowest bin that will be fit
  Int_t fHighFitBin; // Highest bin that will be fit
  Int_t fSystemIndex; // Index of the system (based on enum in Main.cxx)

};


#endif
