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
  
  /* void CreateAllPairSystems(Int_t configuration); */
  void CreatePairSystem(TString simpleName, TString fileName, TString histName, Bool_t isPrimaryIdentical, vector<Double_t> initParams, vector<Bool_t> fixParams)
  void CreateMinuit();
  void DoFitting();
  Int_t GetNSystems() const {return fNSystems;};
  void SaveOutputPlots();
  void SetFitOptions();
  void SetUseEstimatedLambdaParams(Bool_t useParam);

 private:
  void ConstrainF0D0();
  void ConstrainRadii();
  Double_t GetConstrainedParamIndex(const SystemType sys, const ParamType par);
  void GetHistConfiguration(Int_t config, vector<TString> &fileNames, vector<TString> &histNames);
  void InitializeParameters(TMinuit *minuit);
  Bool_t IsParameterConstrained(const SystemType sys, const ParamType par);
  void SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag);
  void SetupInitialParameters()
  void SetupParameterConstraints(const Int_t config);
  void SetupParameterVectors();

  TMinuit *fMinuit;
  vector<TString>  fMinuitParNames;
  vector<Double_t> fMinuitParInitial;
  /* vector<Double_t> fMinuitParMinimum; */
  /* vector<Double_t> fMinuitParMaximum; */
  /* vector<Double_t> fMinuitParCurrent; */
  vector<Bool_t>   fMinuitParIsFixed;

  vector<*ParameterConstraint> fParamConstraints;
  Double_t fNParams;
  vector<TString> fParamNames;
  Double_t fNSystems;
  /* vector<Bool_t> fAllowImagF0; // Do we allow non-zero ImF0?  For each system */
  Int_t fMaxMinuitCalls;
  Bool_t fUseEstimatedLambdaParams;

  // PairSystem information
  vector<*PairSystem> fPairSystems;
  vector<TString> fSystemNames;
  vector<vector<Double_t> > fInitParams;
  vector<vector<Bool_t> > fFixParams;
  /* Bool_t fFixRadius; */
  /* Bool_t fFixF0Real; */
  /* Bool_t fFixF0Imag; */
  /* Bool_t fFixD0; */
  /* Bool_t fFixNorm; */
  /* Bool_t fFixLambda; */
}

#endif
