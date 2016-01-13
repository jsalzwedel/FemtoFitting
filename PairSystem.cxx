#include "LednickyInfo.h"
#include "LednickyEqn.h"
#include "PairSystem.h"
#include "ParamType.h"

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
  assert(combinedGraph);
  //Calculate the total chisquare difference over the fit range.  
  Double_t chi2 = 0.;

  // Find the difference between the fit and the data.  Note that
  // TGraph bins start counting at 0, while histograms bins start
  // counting at 1 (the 0 bin is underflow)
  if(fUseLogLikelihood) {
    // Use loglikelihood fitting of numerators and denominators
    for(UInt_t iHist = 0; iHist < fNHists; iHist++) {
      for(Int_t iBin = fLowFitBin; iBin < fHighFitBin; iBin++) {
	Double_t numCount = fNumHists[iHist]->GetBinContent(iBin+1);
	if(numCount == 0) continue;
	Double_t denCount = fDenHists[iHist]->GetBinContent(iBin+1);
	if(denCount == 0) continue;
	Double_t cfVal = combinedGraph->GetY()[iBin];
	cfVal /= fNorms[iHist];

	Double_t log1st = log(cfVal * (numCount + denCount) /
			      (numCount * (cfVal + 1)));
	Double_t log2nd = log((numCount + denCount) /
			      (denCount * (cfVal + 1)));
	chi2 -= 2. * (numCount * log1st + denCount * log2nd);
	// cout<<"Bin:\t"<<iBin<<"\t\tCF:\t"<<cfVal<<endl;
	// cout<<"Bin:\t"<<iBin<<"\tlog1st:\t"<<log1st
	//     <<"\t\tlog2nd:\t"<<log2nd
	//     <<"\t\tChi2:\t"<<chi2
	//     <<"\tNum:\t"<<numCount
	//     <<"\tDen:\t"<<denCount<<endl;
      }
    }
  } else {
    // Do regular chisquare minimization on the correlation function
    for(Int_t iBin = fLowFitBin; iBin < fHighFitBin; iBin++)
    {
      Double_t cfVal = combinedGraph->GetY()[iBin];
      cfVal /= fNorms[0];
      Double_t diff = cfVal - fCF->GetBinContent(iBin+1);
      Double_t err = fCF->GetBinError(iBin+1);
      if(err < 0.0000001) continue;
      chi2 += pow(diff,2)/pow(err,2);
    }
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

  TString graphName;
  if(fCF) graphName = fCF->GetName();
  else graphName = fPairTypeName;
  graphName += "Fit";
  combinedLednicky->SetName(graphName);
  
  // Initialize the y-values to a unity baseline
  for(Int_t iBin = 0; iBin < combinedLednicky->GetN(); iBin++){
    combinedLednicky->GetY()[iBin] = 1.;
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
  if(fUseQuadraticBkg) {
    // Add a quadratic parameter to account for the background.
    for(Int_t iBin = 0; iBin < combinedLednicky->GetN(); iBin++)
    {
      Double_t kStar = fBinWidth * (1.*iBin + 0.5);
      combinedLednicky->GetY()[iBin] += pow(kStar, 2) * fQuadraticBkg;
    }
  }
  return combinedLednicky;
}

TH1D* PairSystem::GetLednickyHist()
{
  // Turn the combined lednicky TGraph into a histogram
  TGraph *g = GetCombinedTGraph();
  assert(g->GetN() == fNBins);

  //
  Double_t minValue = 0.;
  Double_t maxValue = fNBins * fBinWidth;
  TString histName = g->GetName();
  histName+= "LedHist";
  cout<<"Making a lednicky hist with name:\t"<<histName<<endl;
  TH1D *hist = new TH1D(histName, histName, fNBins, minValue, maxValue);

  //For each point in the TGraph, add it to the histogram
  //No bin should be added to more than once
  Double_t *xVals = g->GetX();
  Double_t *yVals = g->GetY();
  for(Int_t iBin = 0; iBin < fNBins; iBin++)
  {
    Int_t binNumber = hist->FindBin(xVals[iBin]);
    assert(binNumber == iBin +1);
    hist->SetBinContent(binNumber, yVals[iBin]);
  }
  return hist;
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
  // ScaleGraph(primaryGraph, 1./fNorm);
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
    // ScaleGraph(residual, 1./fNorm);
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
  if(!fUseLogLikelihood) fNorms[0] = pars[kNorm];
  else {
    for(UInt_t iHists = 0; iHists < fNHists; iHists++) {
      fNorms[iHists] = pars[iHists + kNorm];
    }
  }
}

PairSystem::PairSystem(TH1D *cfData, const vector<LednickyInfo> &ledInfo, TString pairTypeName, Int_t sysType):
fPairTypeName(pairTypeName),
fSystemType(sysType),
fCF(cfData),
fNHists(0),
fQuadraticBkg(0.),
fUseQuadraticBkg(kFALSE),
fUseLogLikelihood(kFALSE),
fLowFitBin(0),
fHighFitBin(50)
{
  // Constructor for use with chisquare fitting
  assert(cfData);
  fNBins = cfData->GetNbinsX();
  fBinWidth = cfData->GetBinWidth(1);
  fNorms.resize(1, 1.);
  // Create a LednickyEqn from each LednickyInfo
  for(UInt_t iSys = 0; iSys < ledInfo.size(); iSys++)
  {
    fLambdaParameters.push_back(ledInfo[iSys].GetLambdaParam());
    LednickyEqn *lednicky = new LednickyEqn(ledInfo[iSys], fNBins, fBinWidth);
    fLednickyEqns.push_back(lednicky);
  }
}

PairSystem::PairSystem(vector<TH1D*> numHists, vector<TH1D*> denHists, const vector<LednickyInfo> &ledInfo, TString pairTypeName, Int_t sysType):
fPairTypeName(pairTypeName),
fSystemType(sysType),
fCF(NULL),
fNumHists(numHists),
fDenHists(denHists),
fQuadraticBkg(0.),
fUseQuadraticBkg(kFALSE),
fUseLogLikelihood(kTRUE),
fLowFitBin(0),
fHighFitBin(50)
{
  // Constructor for log-likelihood fitting using numerator and
  // denominator distributions
  fNHists = fNumHists.size();
  assert(fDenHists.size() == fNHists);
  assert(fNumHists[0]);
  fNBins = fNumHists[0]->GetNbinsX();
  fBinWidth = fNumHists[0]->GetBinWidth(1);
  fNorms.resize(fNHists, 1.);
  for(UInt_t i = 0; i < fNHists; i++) {
    assert(fNumHists[i]);
    assert(fDenHists[i]);
    assert(fNumHists[i]->GetNbinsX() == fNBins);
    assert(fDenHists[i]->GetNbinsX() == fNBins);
    assert(fNumHists[i]->GetBinWidth(1) == fBinWidth);
    assert(fDenHists[i]->GetBinWidth(1) == fBinWidth);
  }
  fCF = NULL; // Not using the correlation function in this mode
  
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
  for(UInt_t i = 0; i < fLednickyEqns.size(); i++) {
    if(!fLednickyEqns[i]) continue;
    delete fLednickyEqns[i];
    fLednickyEqns[i] = NULL;
  }
  for(UInt_t i = 0; i < fNumHists.size(); i++) {
    if(fNumHists[i]) {
      delete fNumHists[i];
      fNumHists[i] = NULL;
    }
  }
  fNumHists.clear();
  for(UInt_t i = 0; i < fDenHists.size(); i++) {
    if(fDenHists[i]) {
      delete fDenHists[i];
      fDenHists[i] = NULL;
    }
  }
  fDenHists.clear();
}
