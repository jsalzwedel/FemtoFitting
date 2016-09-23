#ifndef LednickyInfo_H
#define LednickyInfo_H

#include "TH2D.h"

class LednickyInfo {
 public: 
  LednickyInfo(TString systemName, Double_t lambdaParam, TH2D *transformMatrix, Bool_t isIdenticalPair, Bool_t useRootSScaling, Double_t baseMass1 = 0, Double_t baseMass2 = 0, Double_t actualMass1 = 0, Double_t actualMass2 = 0, Bool_t fixScatterParams = kFALSE, Double_t reF0 = 0., Double_t imF0 = 0., Double_t d0 = 0., Bool_t fixRadius = kFALSE, Double_t radius = 0.);
  ~LednickyInfo();
  TString GetSystemName() const {return fSystemName;};
  Double_t GetLambdaParam() const {return fLambdaParam;};
  TH2D *GetTransformMatrix() const {return fTransformMatrix;};
  Bool_t GetIsIdenticalPair() const {return fIsIdenticalPair;};
  Bool_t GetUseRootSScaling() const {return fUseRootSScaling;};
  Double_t GetBaseMass1() const {return fBaseMass1;};
  Double_t GetBaseMass2() const {return fBaseMass2;};
  Double_t GetActualMass1() const {return fActualMass1;};
  Double_t GetActualMass2() const {return fActualMass2;};
  Bool_t   GetFixScatterParams() const {return fFixScatterParams;};
  Double_t GetReF0() const {return fReF0;};
  Double_t GetImF0() const {return fImF0;};
  Double_t GetD0() const {return fD0;};
  Double_t GetFixRadius() const {return fFixRadius;};
  Double_t GetRadius() const {return fRadius;};
  

 private:
  TString fSystemName;
  Double_t fLambdaParam;
  TH2D* fTransformMatrix;
  Bool_t fIsIdenticalPair;
  // The variables below are used for sqrt(s) scaling of scatt length
  Bool_t fUseRootSScaling;
  Double_t fBaseMass1; 
  Double_t fBaseMass2;
  Double_t fActualMass1;
  Double_t fActualMass2;
  Bool_t   fFixScatterParams;
  Double_t fReF0; // These are only used if you are manually fixing
  Double_t fImF0; // the scattering parameters for a particular
  Double_t fD0;   // residual correlation
  Bool_t   fFixRadius;
  Double_t fRadius;
  
};



#endif
