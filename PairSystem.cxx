#include "LednickyInfo.h"
#include "LednickyEqn.h"
#include "PairSystem.h"

#include "TStopwatch.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>
#include <assert.h>

using namespace std;

void PairSystem::AddToGraph(TGraph *graph, Double_t d)
{
  for(Int_t iBin = 0; iBin < graph->GetN(); iBin++){
    graph->GetY()[iBin] += d;
  }
}

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
  for(Int_t iBin = fLowFitBin; iBin < fHighFitBin; iBin++)
  {
    Double_t diff = combinedGraph->GetY()[iBin] - fCF->GetBinContent(iBin+1);
    Double_t err = fCF->GetBinError(iBin+1);
    chi2 += pow(diff,2)/pow(err,2);
  }
  delete combinedGraph;
  return chi2;
}


TGraph* PairSystem::GetCombinedTGraph()
{
  // Make a TGraph object that sums all the LednickyEqns with
  // lambda parameters.  

  // Set up the graph's dimensions
  TGraph *combinedLednicky = fLednickyEqns[0]->GetLednickyGraph();

  TString graphName = fCF->GetName();
  graphName += "Fit";
  combinedLednicky->SetName(graphName);
  
  // Initialize the y-values to zero
  for(Int_t iBin = 0; iBin < combinedLednicky->GetN(); iBin++){
    combinedLednicky->GetY()[iBin] = 0.0;
  }


  // Calculate the y-value of the graph for each x-bin.
  // Add the lambda-weighted bin content from each lednicky eqn
  // to the combinedLednicky eqn TGraph
  for(UInt_t iLed = 0; iLed < fLednickyEqns.size(); iLed++)
  {
    TGraph *ledEqn = fLednickyEqns[iLed]->GetLednickyGraph();
    Double_t lambdaParam = fLambdaParameters[iLed];

    for(Int_t iBin = 0; iBin < combinedLednicky->GetN(); iBin++)
    {
      Double_t graphValue = ledEqn->GetY()[iBin];
      Double_t partialBinContent = (graphValue - 1.) * lambdaParam;
      combinedLednicky->GetY()[iBin] += partialBinContent;
    }
    delete ledEqn;
  }

  for(Int_t iBin = 0; iBin < combinedLednicky->GetN(); iBin++)
  {
    combinedLednicky->GetY()[iBin] += 1.;
    combinedLednicky->GetY()[iBin] /= fNorm;

  }
  return combinedLednicky;
}

TCanvas* PairSystem::GetResidualComponentCanvas()
{
  // Create a canvas with all the primary and residual
  // correlation components drawn onto it

  //Prepare the canvas
  TCanvas *can = new TCanvas(fPairTypeName, fPairTypeName);
  TStyle style;
  style.SetPalette(55,0);

  // Prepare legend
  TLegend *leg = new TLegend(0.5, 0.3, 0.9, 0.6);
  leg->SetHeader("Primary and residual correlations");
  
  // Draw the primary-primary correlation function with lambda
  // weight
  TGraph *primaryGraph = fLednickyEqns[0]->GetLednickyGraph();
  Double_t lambda = fLambdaParameters[0];
  AddToGraph(primaryGraph, -1);
  ScaleGraph(primaryGraph, lambda);
  AddToGraph(primaryGraph, 1);
  ScaleGraph(primaryGraph, 1./fNorm);
  primaryGraph->Draw();
  // Make nice axes
  primaryGraph->GetXaxis()->SetRangeUser(0., 1.);
  primaryGraph->GetYaxis()->SetRangeUser(0.85, 1.02);
  primaryGraph->SetLineWidth(2);
  primaryGraph->Draw();

  leg->AddEntry(primaryGraph, "", "L");

  
  // Draw all the residual correlation functions;
  for(UInt_t iLed = 1; iLed < fLednickyEqns.size(); iLed++)
  {
    TGraph *residual = fLednickyEqns[iLed]->GetLednickyGraph();
    Double_t lambdaRes = fLambdaParameters[iLed];
    AddToGraph(residual, -1);
    ScaleGraph(residual, lambdaRes);
    AddToGraph(residual, 1);
    ScaleGraph(residual, 1./fNorm);
    // change symbol/color
    residual->SetLineWidth(2);
    residual->SetLineColor((iLed+1) % 50);
    residual->SetLineStyle((iLed+1) % 11);
    residual->Draw("same");
    leg->AddEntry(residual, "", "L");
  }
  leg->Draw();
  
  return can;
}

void PairSystem::ScaleGraph(TGraph *graph, Double_t d)
{
  for(Int_t iBin = 0; iBin < graph->GetN(); iBin++){
    graph->GetY()[iBin] *= d;
  }
}

void PairSystem::SetLednickyParameters(vector<Double_t> pars)
{
  // Pass the new fit parameters to each LednickyEqn
  for(UInt_t iLed = 0; iLed < fLednickyEqns.size(); iLed++){
    fLednickyEqns[iLed]->SetParameters(pars);
  }
  // Normalization parameter should be stored in PairSystem, not
  // passed on to Lednicky Eqn
  fNorm = pars[4]; 
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
  
  // Create a LednickyEqn from each LednickyInfo
  for(UInt_t iSys = 0; iSys < ledInfo.size(); iSys++)
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
  for(UInt_t i = 0; i < fLednickyEqns.size(); i++)
  {
    if(!fLednickyEqns[i]) continue;
    delete fLednickyEqns[i];
    fLednickyEqns[i] = NULL;
  }
}
