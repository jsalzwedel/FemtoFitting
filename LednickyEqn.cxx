//*********************************************
// Base class for Lednicky correlation function
// equation objects. Contains the parameters for
// a particular Lednicky Eqn, and can return a 
// TGraph of the eqn.
// Maybe can get and use a residual correlation?
//*******************************************

#include "LednickyEqn.h"
#include "Faddeeva.cc"

LednickyEqn::LednickyEqn(Bool_t isIdentical, TH2D *transformMatrix)
{
  fIsIdentical = isIdentical;
  //If transform matrix is null, this is a primary correlation
  fTransformMatrix = transformMatrix;
 
}

LednickyEqn::~LednickyEqn()
{
  if(fTransformMatrix){
    delete fTransformMatrix;
    fTransformMatrix = NULL;
  }
}

TGraph* LednickyEqn::GetBaseLednickyGraph()
{
  // Do all the math to calculate a base lednicky equation
  // in the k* frame of the parent particles. Return a
  // TGraph of the function.

  TGraph * baseGraph; // = new TGraph(bins, xArray, yArray);

  //..

  return baseGraph;
}

TGraph *LednickyEqn::GetLednickyGraph()
{
  // Used by Fitter.  Get a TGraph of the Lednicky Eqn
  // (in the case of residual correlations, transform it)

  TGraph *lednickyGraph = GetBaseLednickyGraph();
  lednickyGraph = TransformLednickyGraph(lednickyGraph);
  return lednickyGraph;
}

void LednickyEqn::SetParameters(vector<Double_t> pars)
{
  // Set all the fit parameters, like f0 and d0
  //...
}

TGraph* LednickyEqn::TransformLednickyGraph(TGraph *base)
{
  // If the LednickyEqn object is a residual correlation,
  // return a transformed version of the graph.
  if(!fTransformMatrix) return base;

  // Calculate new TGraph using the transformation math
  TGraph *transformedGraph;

  //...

  delete base;
  return transformedGraph;
}


Double_t LednickyEqn::GetLednickyF1(Double_t z)
{
  double result = (1./z)*Faddeeva::Dawson(z);
  return result;
}

