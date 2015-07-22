
#include "LednickyEqn.h"
#include "SystemType.h"
// #include "PairSystem.h"
// #include "ParameterConstraint.h"
// #include "Fitter.h"
#include <iostream>
#include "TH2D.h"
#include "TFile.h"
#include <vector>




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


void UserSetupSystems(Fitter &fitter)
{
  // The user should modify this to suit their fitting needs

  //************* Add systems to the analysis ****************
  TString fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root";
  TString histName = "CombinedLLAA0-10KstarMomCorrected";
  TString simpleName = "LL010";
  // Make initial parameters: Radius, ReF0, ImF0, D0, Normalization 
  Double_t initParamsArr[5] = {3., -1., 0., 3., 1.}; 
  vector<Double_t> initParams(initParamsArr);

  // Determine which parameters should be fixed in the fitter.
  Bool_t allowImaginaryF0 = kFALSE;
  Bool_t fixParamsArr[5] = {kFALSE, kFALSE, allowImaginaryF0, kFALSE, kFALSE};
  vector<Bool_t> fixParams(fixParamsArr);

  // Does the primary-primary correlation function consist of
  // identical particles (i.e. both baryons or both antibaryons?)
  
  // Create SystemInfo for each primary and residual correlation
  // Args: TString name, Double_t lambdaParamter, TH2D *transformMatrix, Bool_t isIdenticalPair
  vector<SystemInfo*> systems;
  SystemInfo *sysLL = new System("LambdaLambda", 0.28, NULL, kTRUE);
  systems.push_back(systems);
  
  TString fileNameMatrix = "~/Analysis/lambda/AliAnalysisLambda/Fitting/FemtoFitting/PreparedTransformMatrices.root";
  SystemInfo *sysLS = new System("LambdaSigma", 0.21, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaLambda"), kFALSE);
  systems.push_back(systems);

  SystemInfo *sysLX0 = new System("LambdaXi0", 0.14, GetTransformMatrix(fileNameMatrix, "TransformMatrixXi0Lambda"), kFALSE);
  systems.push_back(systems);

  SystemInfo *sysLXC = new System("LambdaXiC", 0.14, GetTransformMatrix(fileNameMatrix, "TransformMatrixXiCLambda"), kFALSE);
  systems.push_back(systems);

  SystemInfo *sysSS = new System("SigmaSigma", 0.04, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaSigma"), kTRUE);
  systems.push_back(systems);

  SystemInfo *sysSX0 = new System("SigmaXi0", 0.05, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaXi0"), kFALSE);
  systems.push_back(systems);

  SystemInfo *sysSXC= new System("SigmaXiC", 0.05, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaXiC"), kFALSE);
  systems.push_back(systems);

  SystemInfo *sysX0XC = new System("Xi0XiC", 0.04, GetTransformMatrix(fileNameMatrix, "TransformMatrixXiCXi0"), kFALSE);
  systems.push_back(systems);


  // Bool_t isPrimaryIdentical = kTRUE;

  fitter.CreatePairSystem(simpleName, fileName, histName, isPrimaryIdentical, initParams, fixParams);
  

  //************* Add more systems as needed ******************

}

TH2D *GetTransformMatrix(TString rootFileName, TString histName)
{
  TFile f(rootFileName,"read");
  TH2D *h = (TH2D*) f.Get(histName);
  assert(h);
  histName->SetDirectory(0);
  f.Close();
  return h;
}

// main for full program
int main(int argc, char **argv)
{
  


  // Setup
  Fitter myFitter();
  UserSetupSystems(myFitter);
  // Add more systems as needed, either here or in UserSetupSystems
  if(myFitter.GetNSystems() < 1) {
    cout<<"No systems to fit."<<endl;
    return 1;
  }
  myFitter.SetFitOptions();

  // Fitting
  myFitter.CreateMinuit();
  myFitter.DoFitting();
  myFitter.SaveOutputPlots();
  
  return 0;
}
 
