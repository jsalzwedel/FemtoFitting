
// #include "LednickyEqn.h"
#include "LednickyInfo.h"
// #include "PairSystem.h"
// #include "ParameterConstraint.h"
#include "Fitter.h"
#include <iostream>
#include "TMinuit.h"
#include "TH2D.h"
#include "TFile.h"
#include <vector>
#include <assert.h>

using namespace std;

Fitter *myFitter = NULL;


// int main()
// {

//   TFile transFile("TestTransforms.root");
//   TH2D *transform = (TH2D*) transFile.Get("TestTrans1");
//   // TH2D *transform = NULL;
//   Bool_t isIdentical = false;
//   TString name = "TestRes1";
//   Int_t nBins = 200;
//   Double_t binWidth = 0.01;
//   LednickyEqn ledEqn(name, isIdentical, transform, nBins, binWidth);
//   Double_t radius = 1.;
//   Double_t rf0 = -1.;
//   Double_t if0 = 1.;
//   Double_t d0 = 3.;
//   Double_t parsArr[4] = {radius, rf0, if0, d0};
//   vector<Double_t> pars(parsArr, parsArr+4);
//   ledEqn.SetParameters(pars);
//   TGraph *myGraph = ledEqn.GetLednickyGraph();

//   TFile outfile("testFile.root", "update");
//   myGraph->Write(myGraph->GetName(), TObject::kOverwrite);
  
//   return 0;
// }


TH2D *GetTransformMatrix(TString rootFileName, TString histName)
{
  TFile f(rootFileName,"read");
  TH2D *h = (TH2D*) f.Get(histName);
  assert(h);
  h->SetDirectory(0);
  f.Close();
  return h;
}

void UserSetupSystems(Fitter *fitter)
{
  // The user should modify this to suit their fitting needs

  //************* Add systems to the analysis ****************
  TString fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root";
  TString histName = "CombinedLLAA0-10KstarMomCorrected";
  TString simpleName = "LL010";
  // Make initial parameters: Radius, ReF0, ImF0, D0, Normalization 
  Double_t initParamsArr[5] = {3., -1., 0., 3., 1.}; 
  vector<Double_t> initParams(initParamsArr, initParamsArr+5);
  Double_t minParamsArr[5] = {0., 0., 0., 0., 0.};
  vector<Double_t> minParams(minParamsArr, minParamsArr+5);
  Double_t maxParamsArr[5] = {0., 0., 0., 0., 0.};
  vector<Double_t> maxParams(maxParamsArr, maxParamsArr+5);  

  // Determine which parameters should be fixed in the fitter.
  Bool_t allowImaginaryF0 = kFALSE;
  Bool_t fixParamsArr[5] = {kFALSE, kFALSE, allowImaginaryF0, kFALSE, kFALSE};
  vector<Bool_t> fixParams(fixParamsArr, fixParamsArr+5);


  // Does the primary-primary correlation function consist of
  // identical particles (i.e. both baryons or both antibaryons?)
  
  // Create LednickyInfo for each primary and residual correlation
  // Args: TString name, Double_t lambdaParamter, TH2D *transformMatrix, Bool_t isIdenticalPair
  vector<LednickyInfo> ledInfo;
  LednickyInfo infoLL("LambdaLambda", 0.28, NULL, kTRUE);
  ledInfo.push_back(infoLL);
  
  TString fileNameMatrix = "~/Analysis/lambda/AliAnalysisLambda/Fitting/FemtoFitting/PreparedTransformMatrices.root";
  LednickyInfo infoLS("LambdaSigma", 0.21, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaLambda"), kFALSE);
  ledInfo.push_back(infoLS);

  LednickyInfo infoLX0("LambdaXi0", 0.14, GetTransformMatrix(fileNameMatrix, "TransformMatrixXi0Lambda"), kFALSE);
  ledInfo.push_back(infoLX0);

  LednickyInfo infoLXC("LambdaXiC", 0.14, GetTransformMatrix(fileNameMatrix, "TransformMatrixXiCLambda"), kFALSE);
  ledInfo.push_back(infoLXC);

  LednickyInfo infoSS("SigmaSigma", 0.04, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaSigma"), kTRUE);
  ledInfo.push_back(infoSS);

  LednickyInfo infoSX0("SigmaXi0", 0.05, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaXi0"), kFALSE);
  ledInfo.push_back(infoSX0);

  LednickyInfo infoSXC("SigmaXiC", 0.05, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaXiC"), kFALSE);
  ledInfo.push_back(infoSXC);

  LednickyInfo infoX0XC("Xi0XiC", 0.04, GetTransformMatrix(fileNameMatrix, "TransformMatrixXiCXi0"), kFALSE);
  ledInfo.push_back(infoX0XC);


  // Bool_t isPrimaryIdentical = kTRUE;

  fitter->CreatePairSystem(simpleName, fileName, histName, ledInfo, initParams, minParams, maxParams, fixParams);
  

  //************* Add more systems as needed ******************

}

void UserSetFitOptions(Fitter *myFitter)
{
  myFitter->SetHighFitBin(50);

}


void MinuitFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag)
{
  myFitter->SetParametersAndFit(i, totalChisquare, par);
  return;
}

// main for full program
int main(int argc, char **argv)
{
  


  // Setup
  myFitter = new Fitter();
  UserSetupSystems(myFitter);
  // Add more systems as needed, either here or in UserSetupSystems
  if(myFitter->GetNSystems() < 1) {
    cout<<"No systems to fit."<<endl;
    return 1;
  }
  UserSetFitOptions(myFitter);
  // myFitter->SetFitOptions();


  // Create a TMinuit object.
  // Define, set and fix parameters.
  const Int_t nParams = myFitter->GetNMinuitParams();
  TMinuit *myMinuit = new TMinuit(nParams);
  myMinuit->SetFCN(MinuitFit);
  myFitter->InitializeMinuitParameters(myMinuit);
  // myFitter->CreateMinuit();
  myFitter->DoFitting(myMinuit);
  myFitter->SaveOutputPlots();
  
  return 0;
}
 
