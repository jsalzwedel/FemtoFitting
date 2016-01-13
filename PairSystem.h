//************************************************
// Class to contain all the fitting info for a 
// given pair type / centrality combination.
// One of these will be created by the Fitter
// for each pair type / centrality being fit.
//************************************************

#ifndef PairSystem_H
#define PairSystem_H

#include "LednickyInfo.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"

#include <vector>

using std::vector;

class LednickyEqn;

class PairSystem{
 public:
  PairSystem(TH1D *cfData, const vector<LednickyInfo> &ledInfo, TString pairTypeName, Int_t sysType);
  PairSystem(vector<TH1D*> numHists, vector<TH1D*> denHists, const vector<LednickyInfo> &ledInfo, TString pairTypeName, Int_t sysType);
  ~PairSystem();
  Double_t CalculateFitChisquare();
  TH1D*    GetCF() {return fCF;};
  TGraph*  GetCombinedTGraph();
  TH1D*    GetLednickyHist();
  TCanvas* GetResidualComponentCanvas();
  Int_t    GetSystemType() const {return fSystemType;};
  void SetCFData(TH2D* correlationFunction);
  void SetLednickyParameters(vector<Double_t> pars);
  void SetLowFitBin(Int_t bin) {fLowFitBin = bin;};
  void SetHighFitBin(Int_t bin) {fHighFitBin = bin;};
  void SetLowBkgFitBin(Int_t bin) {fLowFitBin = bin;};
  void SetHighBkgFitBin(Int_t bin) {fHighFitBin = bin;};
  void SetSystemType(Int_t sysType) {fSystemType = sysType;};
  void SetUseQuadraticBkg(Bool_t shouldUse) {fUseQuadraticBkg = shouldUse;};

 private:
  void AddToGraph(TGraph *graph, Double_t d);
  void ScaleGraph(TGraph *graph, Double_t d);
  /* void CreateDefaultLednickyEqns(); */
  /* void CreateNewLednickyEqn(TString name, Bool_t isIdentical, TH2D *transformMatrix); */
  /* void ReadInLambdaParams(); */
  TString fPairTypeName; // e.g. "LambdaLambda" or "LambdaAntilambda"
  /* TString fCentralityName; // e.g. "0-10" */
  Int_t fSystemType; // Index of the system (based on enum in Main.cxx)
  vector<LednickyEqn*> fLednickyEqns; // List of primary and secondary L&L eqn objects
  vector<Double_t> fLambdaParameters; // Lambda parameters associated with the L&L eqn objects
  TH1D *fCF; // Correlation function data for the pair type
  vector<TH1D*> fNumHists; // Collection of numerator histograms for likelihood fitting
  vector<TH1D*> fDenHists; // Collection of denominator histograms for likelihood fitting
  UInt_t fNHists; // Number of num-den hist pairs that are being fit
  vector<Double_t> fNorms; // Normalization parameters for combined Lednicky graph.  Each CF or Num/Den pair has a norm parameter.
  Double_t fQuadraticBkg; // Quadratic background parameter
  Bool_t fUseEstimatedLambdaParams;
  Bool_t fUseQuadraticBkg;
  Bool_t fUseLogLikelihood;
  Int_t fNBins; // Number of kstar bins in correlation function
  Double_t fBinWidth; // Width of each kstar bin
  Int_t fLowFitBin; // Lowest bin that will be fit
  Int_t fHighFitBin; // Highest bin that will be fit
  Int_t fLowBkgFitBin; // Lowest background bin that will be fit
  Int_t fHighBkgFitBin; // Highest background bin that will be fit
  

};


#endif
