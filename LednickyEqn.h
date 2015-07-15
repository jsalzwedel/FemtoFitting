//*********************************************
// Base class for Lednicky correlation function
// equation objects. Contains the parameters for
// a particular Lednicky Eqn, and can return a 
// TGraph of the eqn.
// Maybe can get and use a residual correlation?
//*******************************************

#ifndef LednickyEqn_H
#define LednickyEqn_H

// #include ...


class LednickyEqn{
 public:
  LednickyEqn(TString name, Bool_t isIdentical, TH2D *transformMatrix);
  virtual ~LednickyEqn();
  TGraph *GetLednickyGraph();
  void SetParameters(const vector<Double_t> &pars);
  Double_t GetLambdaParam() {return fLambda;};
  //Other setters/getters?

 private:
  TString  fName;   // Type of pair this interaction describes
  Double_t fF0Real; // Real part of the scattering length
  Double_t fF0Imag; // Imaginary part of the scattering length
  Double_t fD0;     // Effective range of interaction
  Double_t fRadius; // Homogeneity radius
  Bool_t   fIsIdentical; // Are these pairs identical particles?
  Int_t    fNBins;  // Number of kstar bins in correlation function
  Double_t fBinWidth; // Width of each kstar bin
  TH2D*    fTransformMatrix; // Transform matrix for residual correlations.  Will be null for primary correlation

  Double_t GetLednickyF1(Double_t z); // Calculate the F1 function
  Double_t GetLednickyF2(Double_t z); // Calculate the F2 function
  TGraph *GetBaseLednickyGraph(); //Calculate Lednicky in parent k* frame
  TGraph *TransformLednickyGraph(TGraph *base);
  

};


#endif
