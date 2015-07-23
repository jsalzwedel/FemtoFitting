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
  
  enum ParamType  {kRad    = 0,
  		   kF0Real = 1,
  		   kF0Imag = 2,
  		   kD0     = 3,
  		   kNorm   = 4,
  		   kLambda = 5};
  Fitter();
  virtual ~Fitter();
  
  /* void CreateAllPairSystems(Int_t configuration); */
  void CreatePairSystem(TString simpleName, TString fileName, TString histName, const vector<LednickyInfo> &ledInfo, vector<Double_t> initParams, vector<Double_t> minParams, vector<Double_t> maxParams, vector<Bool_t> fixParams);
  /* void CreateMinuit(); */
  void DoFitting(TMinuit *minuit);
  Int_t GetNMinuitParams() const {return fMinuitParNames.size();};
  Int_t GetNSystems() const {return fNSystems;};
  Int_t GetNParams() const {return fNParams;};
  Int_t GetTotalParams() const {return fNSystems * fNParams;};
  void InitializeMinuitParameters(TMinuit *minuit);
  void PassLednickyParameters(Double_t pairSystemPars);
  void SetHighFitBin(Int_t bin);
  void SetLowFitBin(Int_t bin);
  void SaveOutputPlots();
  /* void SetFitOptions(); */
  void SetParametersAndFit(Int_t& i, Double_t &totalChisquare, Double_t *par);

  /* void SetUseEstimatedLambdaParams(Bool_t useParam); */

 private:
  void ConstrainF0D0();
  void ConstrainRadii();
  Int_t GetConstrainedParamIndex(const Int_t currentSys, const Int_t currentPar);
  Double_t GetChisquarePerNDF();
  /* void GetHistConfiguration(Int_t config, vector<TString> &fileNames, vector<TString> &histNames); */
  Bool_t IsParameterConstrained(const Int_t sys, const Int_t par);
  void SetupInitialParameters();
  void SetupParameterConstraints(const Int_t config);
  void Timer();
  /* void SetupParameterVectors(); */

  TMinuit *fMinuit;
  vector<TString>  fMinuitParNames;
  vector<Double_t> fMinuitParInitial;
  vector<Double_t> fMinuitParMinimum; 
  vector<Double_t> fMinuitParMaximum;
  /* vector<Double_t> fMinuitParCurrent; */
  vector<Bool_t>   fMinuitParIsFixed;

  vector<ParameterConstraint*> fParamConstraints;
  Double_t fNParams;
  Int_t fFixedParams;
  Int_t fFitBins;
  Int_t fHighFitBin;
  Int_t fLowFitBin;
  vector<TString> fParamNames;
  Double_t fNSystems;
  /* vector<Bool_t> fAllowImagF0; // Do we allow non-zero ImF0?  For each system */
  Int_t fMaxMinuitCalls;
  Bool_t fUseEstimatedLambdaParams;

  // PairSystem information
  vector<PairSystem*> fPairSystems;
  vector<TString> fSystemNames;
  vector<vector<Double_t> > fInitParams;
  vector<vector<Double_t> > fMinParams;
  vector<vector<Double_t> > fMaxParams;
  vector<vector<Bool_t> > fFixParams;

  Int_t fFitCalls;
  Double_t fChisquare;
  
  //Debug
  TStopwatch *fTimer;
};

#endif
