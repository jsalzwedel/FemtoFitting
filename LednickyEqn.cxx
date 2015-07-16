//*********************************************
// Base class for Lednicky correlation function
// equation objects. Contains the parameters for
// a particular Lednicky Eqn, and can return a 
// TGraph of the eqn.
// Maybe can get and use a residual correlation?
//*******************************************

#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "LednickyEqn.h"
#include "Faddeeva.hh"
#include <iostream>

using namespace std;

LednickyEqn::LednickyEqn(TString name, Bool_t isIdentical, TH2D *transformMatrix, Int_t nBins, Double_t binWidth)
{
  fName = name;
  fIsIdentical = isIdentical;
  //If transform matrix is null, this is a primary correlation
  fTransformMatrix = transformMatrix;
  fF0Real = 0.;
  fF0Imag = 0.;
  fD0 = 0.;
  fRadius = 0.;
  fNBins = nBins;
  fBinWidth = binWidth;
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
  
  // Arrays to store X and Y values of correlation function
  Double_t *kstar = new Double_t[fNBins];
  Double_t *cf    = new Double_t[fNBins];

  for(Int_t iBin = 0; iBin < fNBins; iBin++)
  {
    // Calculate the value of the correlation function at each bin

    // Use the center of the bin
    kstar[iBin] = fBinWidth * (1.*iBin + 0.5); 

    // Calculate the denominator of the scattering amplitude
    Double_t scatterAmpDen = 0;
    scatterAmpDen += pow((1. + fF0Imag*kstar[iBin]/HbarC()),2)
                   + pow((fF0Real * kstar[iBin] / HbarC()),2);
    scatterAmpDen += pow(kstar[iBin]/HbarC(),4) * pow(fD0,2) 
                   * (pow(fF0Real,2) + pow(fF0Imag,2)) / 4.;
                  
    scatterAmpDen += fF0Real * fD0 * pow(kstar[iBin] / HbarC(), 2);
    
    // Calculate the real and imaginary parts of the scattering
    // amplitude numerator
    Double_t scatterAmpRe = fF0Real
      + 0.5 * fD0 * pow(kstar[iBin]/HbarC(),2) * (pow(fF0Real,2) + pow(fF0Imag,2));
    scatterAmpRe /= scatterAmpDen;
    Double_t scatterAmpIm = fF0Imag 
      + (pow(fF0Real,2) + pow(fF0Imag,2)) * kstar[iBin] / HbarC();
    scatterAmpIm /= scatterAmpDen;
    Double_t scatterAmpMagSqr = (pow(scatterAmpRe,2) + pow(scatterAmpIm,2));
    
    // Finally, calculate the cf value
    cf[iBin] = 0.5 * (scatterAmpMagSqr/pow(fRadius,2)) 
      * (1. - fD0/(2. * sqrt(TMath::Pi()) * fRadius));
    cf[iBin] += 2. * scatterAmpRe 
      * GetLednickyF1(2. * kstar[iBin] * fRadius / HbarC())
      / (sqrt(TMath::Pi())*fRadius);
    cf[iBin] -= scatterAmpIm 
      * GetLednickyF2(2. * kstar[iBin] * fRadius / HbarC())
      / fRadius;
    if(fIsIdentical) {
      cf[iBin] *= 0.5;
      cf[iBin] -= 0.5 * exp(-4. * pow(kstar[iBin] * fRadius/HbarC(),2));
    }
    cf[iBin] += 1.;
  }
  
  // Now make a TGraph of the correlation function
  TGraph * baseGraph = new TGraph(fNBins, kstar, cf);
  delete[] kstar;
  delete[] cf;
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

void LednickyEqn::SetParameters(const vector<Double_t> pars)
{
  // Set all the fit parameters
  fRadius = pars[0]; 
  fF0Real = pars[1];
  fF0Imag = pars[2];
  fD0     = pars[3];     
}

TGraph* LednickyEqn::TransformLednickyGraph(TGraph *base)
{
  // If the LednickyEqn object is a residual correlation,
  // return a graph that has been transformed into the 
  // relative momentum space of the primary correlations.

  if(!fTransformMatrix) return base;

  TGraph *transformedGraph = (TGraph*) base->Clone("transformedGraph");

  const Int_t nBins = transformedGraph->GetN();
  
  // DaughterBins are in the relative momentum space of the 
  // primary correlations
  for(Int_t daughterBin = 0; daughterBin < nBins; daughterBin++)
  {
    Double_t weightSum = 0.;
    Double_t valueSum = 0.;
    // ParentBins are in the relative momentum space of the
    // residual correlation
    for(Int_t parentBin = 0; parentBin < nBins; parentBin++)
    {
      Double_t weight = fTransformMatrix->GetBinContent(daughterBin+1, parentBin+1);
      weightSum += weight;
      valueSum += weight * base->GetY()[parentBin];
    }
    if(weightSum < 0.99) transformedGraph->GetY()[daughterBin] = 0;
    else transformedGraph->GetY()[daughterBin] = valueSum/weightSum;
  }
  
  delete base;
  return transformedGraph;
}


Double_t LednickyEqn::GetLednickyF1(Double_t z)
{
  double result = (1./z)*Faddeeva::Dawson(z);
  return result;
}

Double_t LednickyEqn::GetLednickyF2(Double_t z)
{
  TF1 lednickyF2("lednickyf2","(1-exp(-x*x))/(x)");
  return lednickyF2.Eval(z);
}

