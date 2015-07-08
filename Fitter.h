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
  void SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag);
  void GetHistConfiguration(Int_t config, vector<TString> &fileNames, vector<TString> &histNames);
  TMinuit *fMinuit;
  vector<*PairSystem> fPairSystems;
  vector<TString>  parNames;
  vector<Double_t> parInitial;
  vector<Double_t> parMinimum;
  vector<Double_t> parMaximum;
  vector<Bool_t>   parIsFixed;
  Int_t fMaxMinuitCalls;
}
