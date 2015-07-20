
#include "LednickyEqn.h"
// #include "PairSystem.h"
// #include "ParameterConstraint.h"
// #include "Fitter.h"
#include <iostream>
#include "TH2D.h"
#include "TFile.h"
#include <vector>



int main()
{

  TFile transFile("TestTransforms.root");
  TH2D *transform = (TH2D*) transFile.Get("TestTrans1");
  // TH2D *transform = NULL;
  Bool_t isIdentical = false;
  TString name = "TestRes1";
  Int_t nBins = 200;
  Double_t binWidth = 0.01;
  LednickyEqn ledEqn(name, isIdentical, transform, nBins, binWidth);
  Double_t radius = 1.;
  Double_t rf0 = -1.;
  Double_t if0 = 1.;
  Double_t d0 = 3.;
  Double_t parsArr[4] = {radius, rf0, if0, d0};
  vector<Double_t> pars(parsArr, parsArr+4);
  ledEqn.SetParameters(pars);
  TGraph *myGraph = ledEqn.GetLednickyGraph();

  TFile outfile("testFile.root", "update");
  myGraph->Write(myGraph->GetName(), TObject::kOverwrite);
  
  return 0;
}


void UserSetupSystems(Fitter &fitter)
{
  // The user should modify this to suit their fitting needs
  
  // Int_t systems = Fitter::kLL010 
  //   + Fitter::kLL1030
  //   + Fitter::kLL3050
  //   + Fitter::kLA010
  //   + Fitter::kLA1030
  //   + Fitter::kLA3050;

  // Add systems to the analysis

  TString fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root";
  TString histName = "CombinedLLAA0-10KstarMomCorrected";
  TString simpleName = "LL010";
  Bool_t isPrimaryIdentical = kTRUE;
  Bool_t allowImaginaryF0 = kFALSE;
  Double_t initParamsArr[5] = {3., -1., 0., 3., 1.}; 
  vector<Double_t> initParams(initParamsArr);
  // Make initial parameters: Radius, ReF0, ImF0, D0, Normalization 
  Bool_t fixParamsArr[5] = {kFALSE, kFALSE, allowImaginaryF0, kFALSE, kFALSE};
  vector<Bool_t> fixParams(fixParamsArr);

  fitter.CreatePairSystem(simpleName, fileName, histName, isPrimaryIdentical, initParams, fixParams);

  // Add more as needed

}


/*
// main for full program
int main(int argc, char **argv)
{
  
  // Check that the first argument (configuration number)
  // is an integer.
  Int_t config;
  if (sscanf (argv[1], "%i", &config)!=1) { printf ("error - not an integer"); }

  // Setup
  Fitter myFitter();
  myFitter.CreateAllPairSystems(config);
  myFitter.SetFitOptions();

  // Fitting
  myFitter.CreateMinuit();
  myFitter.DoFitting();
  myFitter.SaveOutputPlots();
  
  
}
 */
