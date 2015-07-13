//**********************************************************
// Fitting correlation functions
//**********************************************************

#ifndef Fitter_H
#define Fitter_H



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
  TMinuit *fMinuit;
  vector<*PairSystem> fPairSystems;
  vector<TString>  fParNames;
  vector<Double_t> fParInitial;
  vector<Double_t> fParMinimum;
  vector<Double_t> fParMaximum;
  vector<Double_t> fParCurrent;
  vector<Bool_t>   fParIsFixed;
  Double_t fNParams;
  Double_t fNSystems;
  Int_t fMaxMinuitCalls;
}
