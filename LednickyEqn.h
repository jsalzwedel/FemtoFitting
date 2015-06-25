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
  LednickyEqn( /*Input Parameters*/ );
  virtual ~LednickyEqn();
  TGraph *GetLednickyGraph();
  void SetParameters(/*Parameters*/);
  
  //Setters?

 private:
  Double_t fF0Real; // Real part of the scattering length
  Double_t fF0Imag; // Imaginary part of the scattering length
  Double_t fD0;     // Effective range of interaction
  Double_t fRadius; // Homogeneity radius
  Double_t fNorm;   // Normalization parameter
  Double_t fLambda; // Lambda (pair purity) parameter 
  Bool_t   fIsIdentical; // Are these pairs identical particles?
  Double_t fNBins;  // Number of bins in correlation function

  Double_t GetLednickyF1(Double_t z); // Calculate the F1 function


};
