//*********************************************
// Base class for Lednicky correlation function
// equation objects. Contains the parameters for
// a particular Lednicky Eqn, and can return a 
// TGraph of the eqn.
// Maybe can get and use a residual correlation?
//*******************************************

#include "LednickyEqn.h"
#include "Faddeeva.cc"

LednickyEqn::LednickyEqn(TString name, Bool_t isIdentical, TH2D *transformMatrix)
{
  fName = name;
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
    scatterAmpDen += pow((1. + fF0Im*kstar[iBin]/hbarc),2)
                   + pow((fF0Re * kstar[iBin] / hbarc),2);
    scatterAmpDen += pow(kstar,4) * pow(fD0,2) *
                   * (pow(fF0Re,2) + pow(fF0Im,2)) 
                   / (4. * pow(hbarc,4));
    scatterAmpDen += fF0Re * fD0 * pow(kstar[iBin] / hbarc, 2);
    
    // Calculate the real and imaginary parts of the scattering
    // amplitude numerator
    Double_t scatterAmpRe = fF0Real
      + 0.5 * fD0 * pow(kstar[iBin]/hbarc,2) * (pow(fF0Re,2) + pow(fF0Im,2));
    Double_t scatterAmpIm = fF0Im 
      + (pow(fF0Re,2) + pow(fF0Im,2)) * kstar[iBin] / hbarc;
    Double_t scatterAmpMagSqr = (pow(scatterAmpRe,2) + pow(scatterAmpIm,2)) / pow(scatterAmpDen,2);
    
    // Finally, calculate the cf value
    cf[iBin] = 0.5 * (scatterAmpMagSqr/pow(fRadius,2)) 
      * (1. - fD0/(2. * sqrt(TMath::Pi()) * fRadius));
    cf[iBin] += 2. * scatterAmpRe 
      * GetLednickyF1(2. * kstar[iBin] * fRadius)
      / (sqrt(TMath::Pi())*fRadius);
    cf[iBin] -= scatterAmpIm 
      * GetLednickyF2(2. * kstar[iBin] * fRadius)
      / fRadius;
    if(fIsIdentical) {
      cf[iBin] *= 0.5;
      cf[iBin] -= 0.5 * exp(-4. * pow(kstar[iBin] * fRadius,2));
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

void LednickyEqn::SetParameters(const vector<Double_t> &pars)
{
  // Set all the fit parameters
  fLambda = pars[0];
  fRadius = pars[1]; 
  fF0Real = pars[2];
  fF0Imag = pars[3];
  fD0     = pars[4];     
  fNorm   = pars[5];
}

TGraph* LednickyEqn::TransformLednickyGraph(TGraph *base)
{
  // If the LednickyEqn object is a residual correlation,
  // return a graph that has been transformed into the 
  // relative momentum space of the primary correlations.

  if(!fTransformMatrix) return base;

  TGraph *transformedGraph = base->Clone("transformedGraph");

  const Int_t nBins = transformedGraph->GetN();
  
  // DaughterBins are in the relative momentum space of the 
  // primary correlations
  for(Int_t daughterBin = 0; daughterBin < nBin; daughterBin++)
  {
    Double_t weightSum = 0.;
    Double_t valueSum = 0.;
    // ParentBins are in the relative momentum space of the
    // residual correlation
    for(Int_t parentBin = 0; parentBin < nBin; parentBin++)
    {
      Double_t weight = fTransformMatrix->GetBinContent(daughterBin+1, parentBin+1);
      weightSum += weight;
      valueSum += weight * base->GetBinContent(parentBin);
    }
    if(weightSum < 0.99) transformedGraph->GetY[daughterBin] = 0;
    else transformedGraph->GetY[daughterBin] = valueSum/weightSum;
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
  TF1 lednickyF2("lednickyf2","(1-exp(-[0]*[0]*x*x))/([0]*x)");
  const Double_t hbarc = 0.197327;
  lednickyF2.SetParameter(0, 2.*fRadius/hbarc);
  return lednickyF2.Eval(z);
}
