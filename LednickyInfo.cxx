
#include "LednickyInfo.h"

LednickyInfo::LednickyInfo(TString systemName, Double_t lambdaParam, TH2D *transformMatrix, Bool_t isIdenticalPair, Bool_t useRootSScaling, Double_t baseMass1 /*= 0*/, Double_t baseMass2 /*= 0*/, Double_t actualMass1 /*= 0*/, Double_t actualMass2 /*= 0*/, Bool_t fixScatterParams /*= kFALSE*/,  Double_t reF0 /*= 0.*/, Double_t imF0 /*= 0.*/, Double_t d0 /*= 0.*/, Bool_t fixRadius /*=kFALSE*/, Double_t radius /*=0.*/):
  fSystemName(systemName),
  fLambdaParam(lambdaParam),
  fTransformMatrix(transformMatrix),
  fIsIdenticalPair(isIdenticalPair),
  fUseRootSScaling(useRootSScaling),
  fBaseMass1(baseMass1),
  fBaseMass2(baseMass2),
  fActualMass1(actualMass1),
  fActualMass2(actualMass2),
  fFixScatterParams(fixScatterParams),
  fReF0(reF0),
  fImF0(imF0),
  fD0(d0),
  fFixRadius(fixRadius),
  fRadius(radius)
{
}

LednickyInfo::~LednickyInfo()
{
  //Transform matrices will be deleted by their lednicky eqns
}
