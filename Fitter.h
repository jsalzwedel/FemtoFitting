//**********************************************************
// Fitting correlation functions
//**********************************************************

#ifndef Fitter_H
#define Fitter_H

#include "LednickyInfo.h"
#include "PairSystem.h"
#include "ParameterConstraint.h"


#include "TStopwatch.h"
#include "Rtypes.h"
#include "TMinuit.h"


#include <vector>

using std::vector;

class LednickyInfo;


class Fitter{
 public:
  Fitter();
  virtual ~Fitter();
  void AddPairAnalysisChisquareFit(TString simpleName, TString fileName, TString cfName, Int_t sysType, const vector<LednickyInfo> &ledInfo, vector<Double_t> initParams, vector<Double_t> minParams, vector<Double_t> maxParams, vector<Bool_t> fixParams);
  void AddPairAnalysisLogFit(TString simpleName, TString fileName, vector<TString> numNames, vector<TString> denNames, Int_t sysType, const vector<LednickyInfo> &ledInfo, vector<Double_t> initParams, vector<Double_t> minParams, vector<Double_t> maxParams, vector<Bool_t> fixParams);
  void DoFitting();
  Int_t GetNMinuitParams() const {return fMinuitParNames.size();};
  Int_t GetNSystems() const {return fNSystems;};
  /* Int_t GetNParams() const {return fNParams;}; */
  Int_t GetTotalParams() const {return fNParamsTotal;};
  void InitializeMinuitParameters(TMinuit *minuit);
  void SaveOutputPlots();
  void SetDisplayResidualComponents(Bool_t shouldDraw) {fDisplayResidualComponents = shouldDraw;};
  void SetHighFitBin(Int_t bin);
  void SetLowFitBin(Int_t bin);
  void SetMinuitVerbosity(Int_t verbosity);
  void SetOutputString(TString outStr) {fOutputString = outStr;};
  void SetParametersAndFit(Int_t& i, Double_t &totalChisquare, Double_t *par);
  void SetOutputLednickyHists(Bool_t shouldOutput) {fOutputLednickyHist = shouldOutput;};
  void SetStartingStepSize(Double_t step) {fStepSize = step;};
  void SetupConstraint(Int_t param, vector<Int_t> systems);
  void SetupInitialParameterVectors();
  void SetUseMINOS(Bool_t shouldUse) {fUseMINOS = shouldUse;};

 private:
  void CreatePairSystemChisquare(TString simpleName, TString fileName, TString cfName, Int_t sysType, const vector<LednickyInfo> &ledInfo);
  void CreatePairSystemLog(TString simpleName, TString fileName, vector<TString> numNames, vector<TString> denNames, Int_t sysType, const vector<LednickyInfo> &ledInfo);
  Int_t GetConstrainedParamIndex(const Int_t currentSys, const Int_t currentPar);
  Double_t GetChisquarePerNDF();
  Double_t GetPvalue();
  Bool_t IsParameterConstrained(const Int_t sys, const Int_t par);
  void PushBackParams(TString simpleName, vector<Double_t> initParams, vector<Double_t> minParams, vector<Double_t> maxParams, vector<Bool_t> fixParams, UInt_t nNormParams);
  void SaveResidualComponentPlot(Int_t sysIndex);
  void Timer();

  TMinuit *fMinuit;
  vector<TString>  fMinuitParNames;
  vector<Double_t> fMinuitParInitial;
  vector<Double_t> fMinuitParMinimum; 
  vector<Double_t> fMinuitParMaximum;
  vector<Bool_t>   fMinuitParIsFixed;

  vector<ParameterConstraint*> fParamConstraints;
  vector<Double_t> fNParamsPerSystem;
  Int_t fNParamsTotal;
  Int_t fFixedParams;
  Int_t fFitBins;
  Int_t fHighFitBin;
  Int_t fLowFitBin;
  vector<TString> fParamNames;
  Double_t fNSystems;
  Double_t fNLogLikelihoodDataSets;
  Int_t fMaxMinuitCalls;
  Double_t fStepSize;
  Bool_t fUseEstimatedLambdaParams;
  Bool_t fUseMINOS;
  Bool_t fDisplayResidualComponents;
  Bool_t fUseChisquareFitting;
  Bool_t fUseLogLikelihoodFitting;
  Bool_t fOutputLednickyHist;
  TString fOutputString; // Optional suffix for saved object names


  // PairSystem information
  vector<PairSystem*> fPairSystems;
  vector<TString> fSystemNames;
  vector<vector<Double_t> > fInitParams;
  vector<vector<Double_t> > fMinParams;
  vector<vector<Double_t> > fMaxParams;
  vector<vector<Bool_t> > fFixParams;

  Int_t fFitCalls;
  Double_t fChisquare;
  
  // For debug use
  TStopwatch *fTimer;
  Int_t fMinuitVerbosity;
};

#endif
