//**********************************************************
// Fitting correlation functions
//**********************************************************

#ifndef Fitter_H
#define Fitter_H

enum SystemType {kLL010  = 1, 
		 kLL1030 = 2, 
		 kLL3050 = 4,
		 kLA010  = 8,
		 kLA1030 = 16,
		 kLA3050 = 32};

enum ParamType {kF0Real, kF0Imag, kD0, kRad, kNorm, kLambda};


class Fitter{
 public:
  Fitter();
  virtual ~Fitter();
  
  void CreateAllPairSystems(Int_t configuration);
  void CreateMinuit();
  void DoFitting();
  void SaveOutputPlots();
  void SetFitOptions();
  
 private:
  void GetHistConfiguration(Int_t config, vector<TString> &fileNames, vector<TString> &histNames);
  void SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag);
  void InitializeParameters(TMinuit *minuit);
  void SetupParameters();
  void SetupParameterConstraints(const Int_t config);
  void ConstrainF0D0();
  void ConstrainRadii();
  TMinuit *fMinuit;
  vector<*PairSystem> fPairSystems;
  vector<TString>  fParNames;
  vector<Double_t> fParInitial;
  vector<Double_t> fParMinimum;
  vector<Double_t> fParMaximum;
  vector<Double_t> fParCurrent;
  vector<Bool_t>   fParIsFixed;
  vector<*ParameterConstraint> fParamConstraints;
  Double_t fNParams;
  Double_t fNSystems;
  Int_t fMaxMinuitCalls;
}

#endif
