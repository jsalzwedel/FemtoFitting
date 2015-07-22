
#include "SystemInfo.h"

SystemInfo::SystemInfo(TString systemName, Double_t lambdaParam, TH2D *transformMatrix, Bool_t isIdenticalPair):
  fSystemName(systemName),
  fLambdaParam(lambdaParam),
  fTransformMatrix(transformMatrix),
  fIsIdenticalPair(isIdenticalPair)
{
}

SystemInfo::~SystemInfo()
{
  //Transform matrices will be deleted by their lednicky eqns
}
