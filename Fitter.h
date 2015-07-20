//**********************************************************
// Fitting correlation functions
//**********************************************************

#ifndef Fitter_H
#define Fitter_H




class Fitter{
 public:
  enum SystemType {kLL010  = 1 << 0, 
		   kLL1030 = 1 << 1, 
		   kLL3050 = 1 << 2,
		   kLA010  = 1 << 3,
		   kLA1030 = 1 << 4,
		   kLA3050 = 1 << 5};
  
  enum ParamType  {kRad    = 0, 
		   kF0Real = 1,
		   kF0Imag = 2,
		   kD0     = 3,
		   kNorm   = 4,
		   kLambda = 5};
  Fitter();
  virtual ~Fitter();
  
  void CreateAllPairSystems(Int_t configuration);
  void CreateMinuit();
  void DoFitting();
  void SaveOutputPlots();
  void SetFitOptions();
  void SetUseEstimatedLambdaParams(Bool_t useParam);

 private:
  void ConstrainF0D0();
  void ConstrainRadii();
  void GetHistConfiguration(Int_t config, vector<TString> &fileNames, vector<TString> &histNames);
  void InitializeParameters(TMinuit *minuit);
  Bool_t IsParameterConstrained(const SystemType sys, const ParamType par);
  void SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag);
  void SetupInitialParameters()
  void SetupParameterConstraints(const Int_t config);
  void SetupParameterVectors();
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
  Bool_t fUseEstimatedLambdaParams;
}

#endif
