
#include "LednickyInfo.h"
#include "Fitter.h"
#include <iostream>
#include "TMinuit.h"
#include "TH2D.h"
#include "TFile.h"
#include <vector>
#include <assert.h>

using namespace std;

Fitter *myFitter = NULL;


enum ParamType  {kRad    = 0,
		 kF0Real = 1,
		 kF0Imag = 2,
		 kD0     = 3,
		 kNorm   = 4};

enum SystemType {kLL010, kLL1030, kLL3050,
		 kAA010, kAA1030, kAA3050,
		 kLA010, kLA1030, kLA3050,
		 kLLAA010, kLLAA1030, kLLAA3050};


TH2D *GetTransformMatrix(TString rootFileName, TString histName)
{
  TFile f(rootFileName,"read");
  TH2D *h = (TH2D*) f.Get(histName);
  assert(h);
  h->SetDirectory(0);
  f.Close();
  return h;
}

vector<LednickyInfo> PrepareLednickyInfo(Bool_t isIdentical)
{
  // Create LednickyInfo for each primary and residual correlation
  vector<LednickyInfo> ledInfo;

  // Args: TString name, Double_t lambdaParamter, TH2D *transformMatrix, Bool_t isIdenticalPair
  LednickyInfo infoLL("LambdaLambda", 0.25, NULL, isIdentical); 
  ledInfo.push_back(infoLL);
  
  TString fileNameMatrix = "~/Analysis/lambda/AliAnalysisLambda/Fitting/FemtoFitting/PreparedTransformMatrices.root";
  LednickyInfo infoLS("LambdaSigma", 0.194, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaLambda"), kFALSE);
  ledInfo.push_back(infoLS);

  LednickyInfo infoLX0("LambdaXi0", 0.130, GetTransformMatrix(fileNameMatrix, "TransformMatrixXi0Lambda"), kFALSE);
  ledInfo.push_back(infoLX0);

  LednickyInfo infoLXC("LambdaXiC", 0.117, GetTransformMatrix(fileNameMatrix, "TransformMatrixXiCLambda"), kFALSE);
  ledInfo.push_back(infoLXC);

  LednickyInfo infoSS("SigmaSigma", 0.034, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaSigma"), isIdentical);
  ledInfo.push_back(infoSS);

  LednickyInfo infoSX0("SigmaXi0", 0.049, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaXi0"), kFALSE);
  ledInfo.push_back(infoSX0);

  LednickyInfo infoSXC("SigmaXiC", 0.043, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaXiC"), kFALSE);
  ledInfo.push_back(infoSXC);

  LednickyInfo infoX0XC("Xi0XiC", 0.029, GetTransformMatrix(fileNameMatrix, "TransformMatrixXiCXi0"), kFALSE);
  ledInfo.push_back(infoX0XC);

  return ledInfo;
}

void UserSetupSystems(Fitter *fitter)
{
  // Add systems to the analysis. The user should modify this 
  // to suit their fitting needs


  /////////// Setting up Lambda-Lambda + Antilambda-Antilambda //////////////////
  // 0-10%
  TString fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root";
  TString histName = "CombinedLLAA0-10KstarMomCorrected";
  TString simpleName = "LLAA010";
  // Make initial parameters: Radius, ReF0, ImF0, D0, Normalization 
  Double_t radiiParams[3] = {2.95, 2.67, 2.02};
  Double_t initParamsArr[5] = {radiiParams[0], -.63, 0., 1.36, 1.}; 
  vector<Double_t> initParams(initParamsArr, initParamsArr+5);
  Double_t minParamsArr[5] = {0., 0., 0., 0., 0.};
  vector<Double_t> minParams(minParamsArr, minParamsArr+5);
  Double_t maxParamsArr[5] = {0., 0., 0., 0., 0.};
  vector<Double_t> maxParams(maxParamsArr, maxParamsArr+5);  
  // Determine which parameters should be fixed in the fitter.
  Bool_t fixParamsArr[5] = {kFALSE, kFALSE, kTRUE, kFALSE, kFALSE};
  vector<Bool_t> fixParams(fixParamsArr, fixParamsArr+5);
  // Prepare the lednicky eqn info (lambda parameters, transform matrix locations, whether or not particles are identical)
  vector<LednickyInfo> ledInfoLL = PrepareLednickyInfo(kTRUE);
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLLAA010, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 10-30
  histName = "CombinedLLAA10-30KstarMomCorrected";
  simpleName = "LLAA1030";
  initParams[0] = radiiParams[1];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLLAA1030, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 30-50
  histName = "CombinedLLAA30-50KstarMomCorrected";
  simpleName = "LLAA3050";
  initParams[0] = radiiParams[2];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLLAA3050, ledInfoLL, initParams, minParams, maxParams, fixParams);


  ////////////////// Setting up Lambda-Lambda ////////////////////
  // 0-10%
  fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsLamLam.root";
  histName = "LamLam0-10";
  simpleName = "LL010";
  initParams[0] = radiiParams[0];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLL010, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 10-30
  histName = "LamLam10-30";
  simpleName = "LL1030";
  initParams[0] = radiiParams[1];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLL1030, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 30-50
  histName = "LamLam30-50";
  simpleName = "LL3050";
  initParams[0] = radiiParams[2];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLL3050, ledInfoLL, initParams, minParams, maxParams, fixParams);


  ////////////////// Setting up Antilambda-Antilambda ////////////////////
  // 0-10%
  fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsALamALam.root";
  histName = "ALamALam0-10";
  simpleName = "AA010";
  initParams[0] = radiiParams[0];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kAA010, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 10-30
  histName = "ALamALam10-30";
  simpleName = "AA1030";
  initParams[0] = radiiParams[1];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kAA1030, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 30-50
  histName = "ALamALam30-50";
  simpleName = "AA3050";
  initParams[0] = radiiParams[2];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kAA3050, ledInfoLL, initParams, minParams, maxParams, fixParams);




  ////////////////// Setting up Lambda-Antilambda ////////////////////
  // 0-10%
  fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsLamALamKstarMomCorrected.root";
  histName = "LamALam0-10centrality_varBin5BothFieldsKstarMomCorrected";
  simpleName = "LA010";
  Double_t initParamsArrLA[5] = {radiiParams[0], -1., 1., 3., 1.}; 
  vector<Double_t> initParamsLA(initParamsArrLA, initParamsArrLA + 5);
  Bool_t fixParamsArrLA[5] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
  vector<Bool_t> fixParamsLA(fixParamsArrLA, fixParamsArrLA + 5);
  vector<LednickyInfo> ledInfoLA = PrepareLednickyInfo(kFALSE);
  fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLA010, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);

  // 10-30
  histName = "LamALam10-30centrality_varBin5BothFieldsKstarMomCorrected";
  simpleName = "LA1030";
  initParamsLA[0] = radiiParams[1];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLA1030, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);

  // 30-50
  histName = "LamALam30-50centrality_varBin5BothFieldsKstarMomCorrected";
  simpleName = "LA3050";
  initParamsLA[0] = radiiParams[2];
  // fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLA3050, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);
  //************* Add more systems as needed ******************

}

void UserSetConstraints(Fitter *myFitter)
{
  // Set up any constraints between parameters of different systems.
  // Note: Don't constrain any parameters that are fixed.  It
  // screws up the ndf calculation (so chi^2 / ndf will be
  // off a little bit).

  // Share real f0, imaginary f0, and d0 among partical-antiparticle
  Int_t systemsArrLA[3] = {kLA010, kLA1030, kLA3050};
  vector<Int_t> systemsLA(systemsArrLA, systemsArrLA + 3);
  myFitter->SetupConstraint(kF0Real, systemsLA);
  myFitter->SetupConstraint(kF0Imag, systemsLA);
  myFitter->SetupConstraint(kD0, systemsLA);

  // Share real f0, imaginary f0, and d0 among LambdaLambda + AntilambdaAntilambda
  Int_t systemsArrLLAA[3] = {kLLAA010, kLLAA1030, kLLAA3050};
  vector<Int_t> systemsLLAA(systemsArrLLAA, systemsArrLLAA + 3);
  myFitter->SetupConstraint(kF0Real, systemsLLAA);
  // myFitter->SetupConstraint(kF0Imag, systemsLLAA);
  myFitter->SetupConstraint(kD0, systemsLLAA);

  // Share real f0, imaginary f0, and d0 among LambdaLambda
  Int_t systemsArrLL[3] = {kLL010, kLL1030, kLL3050};
  vector<Int_t> systemsLL(systemsArrLL, systemsArrLL + 3);
  // myFitter->SetupConstraint(kF0Real, systemsLL);
  // myFitter->SetupConstraint(kF0Imag, systemsLL);
  // myFitter->SetupConstraint(kD0, systemsLL);

  // Share real f0, imaginary f0, and d0 among LambdaLambda
  Int_t systemsArrAA[3] = {kAA010, kAA1030, kAA3050};
  vector<Int_t> systemsAA(systemsArrAA, systemsArrAA + 3);
  // myFitter->SetupConstraint(kF0Real, systemsAA);
  // myFitter->SetupConstraint(kF0Imag, systemsAA);
  // myFitter->SetupConstraint(kD0, systemsAA);


  // Constrain LL and AA together
  Int_t systemsArrLLandAA[6] = {kLL010, kLL1030,kLL3050,
				kAA010, kAA1030,kAA3050};
  vector<Int_t> systemsLLandAA(systemsArrLLandAA, systemsArrLLandAA + 6);
  // myFitter->SetupConstraint(kF0Real, systemsLLandAA);
  // myFitter->SetupConstraint(kF0Imag, systemsLLandAA);
  // myFitter->SetupConstraint(kD0, systemsLLandAA);

  // Constrain Averaged LLAA with LA
  Int_t systemsArrLLAALA[6] = {kLLAA010, kLLAA1030, kLLAA3050,
			     kLA010, kLA1030, kLA3050};
  vector<Int_t> systemsLLAALA(systemsArrLLAALA, systemsArrLLAALA + 6);
  // myFitter->SetupConstraint(kF0Real, systemsLLAALA);
  // myFitter->SetupConstraint(kF0Imag, systemsLLAALA);
  // myFitter->SetupConstraint(kD0, systemsLLAALA);


  // share radii among like centralities
  vector<Int_t> systems010;
  systems010.push_back(kLLAA010);
  systems010.push_back(kLA010);
  // systems010.push_back(kLL010);
  // systems010.push_back(kAA010);
  // myFitter->SetupConstraint(kRad, systems010);

  vector<Int_t> systems1030;
  systems1030.push_back(kLLAA1030);
  systems1030.push_back(kLA1030);
  // systems1030.push_back(kLL1030);
  // systems1030.push_back(kAA1030);
  // myFitter->SetupConstraint(kRad, systems1030);

  vector<Int_t> systems3050;
  systems3050.push_back(kLLAA3050);
  systems3050.push_back(kLA3050);
  // systems3050.push_back(kLL3050);
  // systems3050.push_back(kAA3050);
  // myFitter->SetupConstraint(kRad, systems3050);
}

void UserSetFitOptions(Fitter *myFitter)
{
  // Set the upper bin of the fit range
   myFitter->SetHighFitBin(40);
  // How big should the initial starting value fit uncertainty be?
  // This is a catchall for *all* parameter step sizes.
   myFitter->SetStartingStepSize(5.);
  // Output extra plots showing all residual correlation components?
   myFitter->SetDisplayResidualComponents(kTRUE);

   // Use MINOS to find error bars?
   myFitter->SetUseMINOS(kTRUE);
   
  // optional suffix for saved plots and objects
  TString outString = "";
  myFitter->SetOutputString(outString); 
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
  UserSetConstraints(myFitter);
  UserSetFitOptions(myFitter);
  // myFitter->SetFitOptions();


  // Create a TMinuit object.
  // Define, set and fix parameters.
  const Int_t nParams = myFitter->GetNMinuitParams();
  TMinuit *myMinuit = new TMinuit(nParams);
  myMinuit->SetFCN(MinuitFit);
  myFitter->InitializeMinuitParameters(myMinuit);
  // myFitter->CreateMinuit();
  myFitter->DoFitting();
  myFitter->SaveOutputPlots();

  delete myMinuit;
  return 0;
}
 
