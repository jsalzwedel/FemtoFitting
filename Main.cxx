
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
