
#include "LednickyInfo.h"

LednickyInfo::LednickyInfo(TString systemName, Double_t lambdaParam, TH2D *transformMatrix, Bool_t isIdenticalPair, Double_t baseMass1, Double_t baseMass2, Double_t actualMass1, Double_t actualMass2):
  fSystemName(systemName),
  fLambdaParam(lambdaParam),
  fTransformMatrix(transformMatrix),
  fIsIdenticalPair(isIdenticalPair),
  fBaseMass1(baseMass1),
  fBaseMass2(baseMass2),
  fActualMass1(actualMass1),
  fActualMass2(actualMass2)
{
}

LednickyInfo::~LednickyInfo()
{
  //Transform matrices will be deleted by their lednicky eqns
}
