#include "LednickyInfo.h"
#include "LednickyEqn.h"
#include "PairSystem.h"

#include "TStopwatch.h"
#include <iostream>
#include <assert.h>

using namespace std;


Double_t PairSystem::CalculateFitChisquare()
{
  // cout<<"PairSystem::CalculateFitChisquare"<<endl;
  // Calculate the chisquare difference between the
  // correlation function data and the combined LednickyEqns
  
  // Update the fit parameters and get the full Lednicky Eqn graph
  TGraph *combinedGraph = GetCombinedTGraph();

  //Calculate the total chisquare difference over the fit range.  
  Double_t chi2 = 0.;

  // Find the difference between the fit and the data.  Note that
  // TGraph bins start counting at 0, while histograms bins start
  // counting at 1 (the 0 bin is underflow)
  // cout<<"fLowFitBin:\t"<<fLowFitBin<<".\tfHighFitBin:\t"<<fHighFitBin<<endl;
  for(Int_t iBin = fLowFitBin; iBin < fHighFitBin; iBin++)
  {
    // cout<<"Graph:\t"<<combinedGraph->GetY()[iBin]
    // 	<<".\tCF:\t"<<fCF->GetBinContent(iBin+1)<<endl;
    Double_t diff = combinedGraph->GetY()[iBin] - fCF->GetBinContent(iBin+1);
    Double_t err = fCF->GetBinError(iBin+1);
    // cout<<"Diff:\t"<<diff<<".\tError:\t"<<err<<endl;
    chi2 += pow(diff,2)/pow(err,2);
  }
  delete combinedGraph;
  // cout<<"Chi2 for these parameters:\t"<<chi2<<endl;
  return chi2;
}


TGraph* PairSystem::GetCombinedTGraph()
{
  // TStopwatch timer;
  // cout<<"PairSystem::GetCombinedTGraph"<<endl;
  // Make a TGraph object that sums all the LednickyEqns with
  // lambda parameters.  

  // Set up the graph's dimensions
  // timer.Print(); timer.Continue(); cout<<"GetGraph"<<endl;
  TGraph *combinedLednicky = fLednickyEqns[0]->GetLednickyGraph();
  // timer.Print(); timer.Continue(); cout<<"Got TGraph"<<endl;

  TString graphName = fCF->GetName();
  graphName += "Fit";
  combinedLednicky->SetName(graphName);
  // Calculate the y-value of the graph for each x-bin



  for(Int_t iLed = 0; iLed < fLednickyEqns.size(); iLed++)
  {
    // TStopwatch subTimer;
    TGraph *ledEqn = fLednickyEqns[iLed]->GetLednickyGraph();
    Double_t lambdaParam = fLambdaParameters[iLed];
    // cout<<"Graph value:\t"<<graphValue
    // <<".\tLambdaParam:\t"<<lambdaParam<<endl;

    for(Int_t iBin = 0; iBin < combinedLednicky->GetN(); iBin++)
    {
      Double_t graphValue = ledEqn->GetY()[iBin];
      Double_t partialBinContent = (graphValue - 1.) * lambdaParam;
      // partialBinContent /= fNorm;
      combinedLednicky->GetY()[iBin] += partialBinContent;

    }
    delete ledEqn;
    // subTimer.Print(); 
    // timer.Print(); timer.Continue();
  }

  for(Int_t iBin = 0; iBin < combinedLednicky->GetN(); iBin++)
  {
    combinedLednicky->GetY()[iBin] += 1.;
    combinedLednicky->GetY()[iBin] /= fNorm;

  }










  // for(Int_t iBin = 0; iBin < combinedLednicky->GetN(); iBin++)
  // {
  //   Double_t combinedBinContent = 1.;
  //   for(Int_t iLed = 0; iLed < fLednickyEqns.size(); iLed++)
  //   {
  //     TStopwatch subTimer;
  //     TGraph *ledEqn = fLednickyEqns[iLed]->GetLednickyGraph();
  //     Double_t graphValue = ledEqn->GetY()[iBin];
  //     Double_t lambdaParam = fLambdaParameters[iLed];
  //     // cout<<"Graph value:\t"<<graphValue
  // 	  // <<".\tLambdaParam:\t"<<lambdaParam<<endl;
  //     combinedBinContent += (graphValue - 1.) * lambdaParam;
  //     delete ledEqn;
  //     subTimer.Print(); 
  //   }
  //   combinedBinContent /= fNorm;
  //   combinedLednicky->GetY()[iBin] = combinedBinContent;
  // // timer.Print(); timer.Continue();
  // }
  // timer.Print(); cout<<"End of GetCombinedTGraph"<<endl;
  return combinedLednicky;
}

void PairSystem::SetLednickyParameters(vector<Double_t> pars)
{
  // cout<<"PairSystem::SetLednickyParameters"<<endl;
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

PairSystem::PairSystem(TH1D *cfData, const vector<LednickyInfo> &ledInfo, TString pairTypeName, Int_t sysType):
fLowFitBin(0),
fHighFitBin(50),
fSystemType(sysType)
{
  assert(cfData);
  fCF = cfData;
  fNBins = cfData->GetNbinsX();
  fBinWidth = cfData->GetBinWidth(1);
  fNorm = 1.;
  fPairTypeName = pairTypeName;
  // fCentralityName = centrality;
  // Create a Lednicky eqn for primary-primary correlation function
  // with no transform matrix (i.e. NULL ptr).

  for(Int_t iSys = 0; iSys < ledInfo.size(); iSys++)
  {
    fLambdaParameters.push_back(ledInfo[iSys].GetLambdaParam());
    LednickyEqn *lednicky = new LednickyEqn(ledInfo[iSys], fNBins, fBinWidth);
    fLednickyEqns.push_back(lednicky);
  }
}



PairSystem::~PairSystem()
{
  if(fCF){
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
