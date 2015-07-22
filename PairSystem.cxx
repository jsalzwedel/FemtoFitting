#include "LednickyInfo.h"
#include "LednickyEqn.h"

Double_t PairSystem::CalculateFitChisquare()
{
  // Calculate the chisquare difference between the
  // correlation function data and the combined LednickyEqns
  
  // Update the fit parameters and get the full Lednicky Eqn graph
  TGraph *combinedGraph = GetCombinedTGraph();

  //Calculate the total chisquare difference over the fit range.  
  Double_t chi2 = 0.;

  // Find the difference between the fit and the data.  Note that
  // TGraph bins start counting at 0, while histograms bins start
  // counting at 1 (the 0 bin is underflow)
  for(Int_t iBin = fLowFitBin; iBin < fHighFitBin; iBin++)
  {
    Double_t diff = combinedGraph->GetY()[iBin] - fCF->GetBinContent(iBin+1);
    Double_t err = fCf->GetBinError(iBin+1);
    chi2 += pow(diff,2)/pow(err,2);
  }
  delete combinedGraph;
  return chi2;
}

// void PairSystem::CreateDefaultLednickyEqns()
// {

  
//   // Make all the residual correlation lednicky eqns
//   TH2D *transformLamSig; // = ... 
//   CreateNewLednickyEqn("LamSig", kFALSE, transformLamSig);

//   TH2D *transformLamXiC; // = ... 
//   CreateNewLednickyEqn("LamXiC", kFALSE, transformLamXiC);

//   TH2D *transformLamXi0; // = ... 
//   CreateNewLednickyEqn("LamXi0", kFALSE, transformLamXi0);

//   TH2D *transformSigSig; // = ... 
//   CreateNewLednickyEqn("SigSig", isIdenticalPrimary, transformSigSig);  // This is identical unless we have particle-antiparticle system

//   TH2D *transformSigXiC; // = ... 
//   CreateNewLednickyEqn("SigXiC", kFALSE, transformSigXiC);
 
//   TH2D *transformSigXi0; // = ... 
//   CreateNewLednickyEqn("SigXi0", kFALSE, transformSigXi0);

// }

// void PairSystem::CreateNewLednickyEqn(TString name, Bool_t isIdentical, TH2D *transformMatrix)
// {
//   // Called by Fitter
//   LednickyEqn *lednicky = new LednickyEqn(name, isIdentical, transformMatrix, fNBins, fBinWidth);
//   fLednickyEqns.push_back(lednicky);
// }

TGraph* PairSystem::GetCombinedTGraph()
{
  // Make a TGraph object that sums all the LednickyEqns with
  // lambda parameters.  

  // Set up the graph's dimensions
  TGraph *combinedLednicky = fLednickyEqns[0]->GetLednickyGraph();

  // Calculate the y-value of the graph for each x-bin
  for(Int_t iBin = 0; iBin < combinedLednicky->GetBinsX(); iBin++)
  {
    Double_t combinedBinContent = 1.;
    for(Int_t iLed = 0; iLed < fLednickyEqns.size(); iLed++)
    {
      TGraph *ledEqn = fLednickyEqns[iLed]->GetLednickyGraph();
      Double_t graphValue = ledEqn->GetY()[iBin];
      combinedBinContent += (graphValue - 1.) * fLambdaParameters(iLed);
      delete ledEqn;
    }
    combinedBinContent /= fNorm;
    combinedLednicky->GetY()[iBin] = combinedBinContent;
  }
  return combinedLednicky;
}

void PairSystem::SetLednickyParameters(vector<Double_t> pars)
{
  // Pass the new fit parameters to each LednickyEqn
  for(Int_t iLed = 0; iLed < fLednickyEqns.size(); iLed++){
    fLednickyEqns[iLed]->SetParameters(pars);
  }
  // Normalization parameter should be stored in PairSystem, not
  // passed on to Lednicky Eqn
  fNorm = pars[4];


  
  // Lambda parameters should not be passed on to Lednicky Eqn.
  // If they are coming from Minuit instead of being read from file
  // they should be updated in the PairSystem here.
  if(pars.size() > 5) 
  {
    // Do lambda param stuff
  }
  
}

PairSystem::PairSystem(TH1D *cfData, const vector<LednickyInfo> &ledInfo, TString pairTypeName)
{
  assert(cfData);
  fCF = cfData;
  fNBins = cfData->GetNBinsX();
  fBinWidth = cfData->GetBinWidth;
  fNorm = 1.;
  fPairTypeName = pairTypeName;
  // fCentralityName = centrality;
  // Create a Lednicky eqn for primary-primary correlation function
  // with no transform matrix (i.e. NULL ptr).

  for(Int_t iSys = 0; iSys < ledInfo.size(); iSys++)
  {
    fLambdaParameters.push_back(ledInfo[iSys]->GetLambdaParam);
    LednickyEqn *lednicky = new LednickyEqn(ledInfo[iSys], fNBins, fBinWidth);
    fLednickyEqns.push_back(lednicky);
  }
}



PairSystem::~PairSystem()
{
  if(fCf){
    delete fCF;
    fCF = NULL;
  }
  for(Int_t i = 0; i < fLednickyEqns.size(); i++)
  {
    if(!fLednickyEqns[i]) continue;
    delete fLednickyEqns[i];
    fLednickyEqns[i] = NULL;
  }
}

// void PairSystem::ReadInLambdaParams()
// {

// }
