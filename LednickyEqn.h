//*********************************************
// Base class for Lednicky correlation function
// equation objects. Contains the parameters for
// a particular Lednicky Eqn, and can return a 
// TGraph of the eqn.
// Maybe can get and use a residual correlation?
//*******************************************

#ifndef LednickyEqn_H
#define LednickyEqn_H

#include "TH2D.h"
#include "TGraph.h"
#include <vector>

class LednickyInfo;

using std::vector;

class LednickyEqn{
 public:
  /* LednickyEqn(TString name, Bool_t isIdentical, TH2D *transformMatrix, Int_t nBins, Double_t binWidth); */
  LednickyEqn(const LednickyInfo &info, Int_t nBins, Double_t binWidth);
  virtual ~LednickyEqn();
  TGraph *GetLednickyGraph();
  void SetParameters(const vector<Double_t> &pars);
  //Other setters/getters?

 private:
  TString  fName;   // Type of pair this interaction describes
  Bool_t   fFixScatterParams; // If on, fitter can't change scatter param values
  Double_t fF0Real; // Real part of the scattering length
  Double_t fF0Imag; // Imaginary part of the scattering length
  Double_t fD0;     // Effective range of interaction
  Bool_t   fFixRadius; // if true, fitter can't change radius
  Double_t fRadius; // Homogeneity radius
  Bool_t   fIsIdentical; // Are these pairs identical particles?
  Bool_t   fUseRootSScaling; // Should scatter amp scale by sqrt s
  Int_t    fNBins;  // Number of kstar bins in correlation function
  Double_t fBinWidth; // Width of each kstar bin
  Double_t fBaseMass1; // Mass of base particle for use with scatter-amp scaling
  Double_t fBaseMass2; // Mass of base particle for use with scatter-amp scaling
  Double_t fActualMass1; // Mass of particle for scaled interaction
  Double_t fActualMass2; // Mass of particle for scaled interaction
  TH2D*    fTransformMatrix; // Transform matrix for residual correlations.  Will be null for primary correlation

  Double_t GetLednickyF1(Double_t z); // Calculate the F1 function
  Double_t GetLednickyF2(Double_t z); // Calculate the F2 function
  TGraph *GetBaseLednickyGraph(); //Calculate Lednicky in parent k* frame
  TGraph *TransformLednickyGraph(TGraph *base);
  Double_t HbarC() const {return 0.197327;};
  

};


#endif
