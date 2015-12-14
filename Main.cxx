
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

vector<LednickyInfo> PrepareLednickyInfo(Bool_t isIdentical, Bool_t useRootSScaling)
{
  // Create LednickyInfo for each primary and residual correlation
  vector<LednickyInfo> ledInfo;

  // Masses for use with sqrt(s) scaling of scattering amplitude.
  // If not using sqrt(s) scaling of interactions, can optionally set these to 0
  Double_t mLambda = 1.1157;
  Double_t mSigma = 1.1193;
  Double_t mXi0 = 1.3149;
  Double_t mXiC = 1.3217;

  // Args: TString name, Double_t lambdaParamter, TH2D *transformMatrix, Bool_t isIdenticalPair, Bool_t useRootSScaling, Double_t baseMass1, Double_t baseMass2, Double_t actualMass1, Double_t actualMass2
  LednickyInfo infoLL("LambdaLambda", 0.25, NULL, isIdentical, useRootSScaling, mLambda, mLambda, mLambda, mLambda); 
  ledInfo.push_back(infoLL);
  
  TString fileNameMatrix = "~/Analysis/lambda/AliAnalysisLambda/Fitting/FemtoFitting/PreparedTransformMatrices.root";
  LednickyInfo infoLS("LambdaSigma", 0.194, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaLambda"), kFALSE, useRootSScaling, mLambda, mLambda, mLambda, mSigma);
  ledInfo.push_back(infoLS);

  LednickyInfo infoLX0("LambdaXi0", 0.130, GetTransformMatrix(fileNameMatrix, "TransformMatrixXi0Lambda"), kFALSE, useRootSScaling, mLambda, mLambda, mLambda, mXi0);
  ledInfo.push_back(infoLX0);

  LednickyInfo infoLXC("LambdaXiC", 0.117, GetTransformMatrix(fileNameMatrix, "TransformMatrixXiCLambda"), kFALSE, useRootSScaling, mLambda, mLambda, mLambda, mXiC);
  ledInfo.push_back(infoLXC);

  LednickyInfo infoSS("SigmaSigma", 0.034, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaSigma"), isIdentical, useRootSScaling, mLambda, mLambda, mSigma, mSigma);
  ledInfo.push_back(infoSS);

  LednickyInfo infoSX0("SigmaXi0", 0.049, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaXi0"), kFALSE, useRootSScaling, mLambda, mLambda, mSigma, mXi0);
  ledInfo.push_back(infoSX0);

  LednickyInfo infoSXC("SigmaXiC", 0.043, GetTransformMatrix(fileNameMatrix, "TransformMatrixSigmaXiC"), kFALSE, useRootSScaling, mLambda, mLambda, mSigma, mXiC);
  ledInfo.push_back(infoSXC);

  LednickyInfo infoX0XC("Xi0XiC", 0.029, GetTransformMatrix(fileNameMatrix, "TransformMatrixXiCXi0"), kFALSE, useRootSScaling, mLambda, mLambda, mXi0, mXiC);
  ledInfo.push_back(infoX0XC);

  return ledInfo;
}

void UserSetupSystems(Fitter *fitter)
{
  // Add systems to the analysis. The user should modify this 
  // function to suit their fitting needs
  Bool_t useLLAA010Chi2 = kFALSE, useLLAA1030Chi2 = kFALSE, useLLAA3050Chi2 = kFALSE,
         useLL010Chi2 = kFALSE, useLL1030Chi2 = kFALSE, useLL3050Chi2 = kFALSE,
         useAA010Chi2 = kFALSE, useAA1030Chi2 = kFALSE, useAA3050Chi2 = kFALSE,
         useLA010Chi2 = kFALSE, useLA1030Chi2 = kFALSE, useLA3050Chi2 = kFALSE;
  // useLLAA010Chi2  = kTRUE;
  // useLLAA1030Chi2 = kTRUE;
  // useLLAA3050Chi2 = kTRUE;
  useLL010Chi2    = kTRUE;
  useLL1030Chi2   = kTRUE;
  useLL3050Chi2   = kTRUE;
  useAA010Chi2    = kTRUE;
  useAA1030Chi2   = kTRUE;
  useAA3050Chi2   = kTRUE;
  useLA010Chi2    = kTRUE;
  useLA1030Chi2   = kTRUE;
  useLA3050Chi2   = kTRUE;
  Bool_t useRootSScalingLL = kFALSE;
  Bool_t useRootSScalingLA = kFALSE;
  
  /////////// Setting up Lambda-Lambda + Antilambda-Antilambda //////////////////
  // 0-10%
  TString fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/CFs.root";
  TString histName = "Merged/CFLLAA010";
  TString simpleName = "LLAA010";
  // Make initial parameters: Radius, ReF0, ImF0, D0, Normalization 
  Double_t radiiParams[3] = {4., 2.67, 2.02};
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
  vector<LednickyInfo> ledInfoLL = PrepareLednickyInfo(kTRUE, useRootSScalingLL);
  if(useLLAA010Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLLAA010, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 10-30
  histName = "Merged/CFLLAA1030";
  simpleName = "LLAA1030";
  initParams[0] = radiiParams[1];
  if(useLLAA1030Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLLAA1030, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 30-50
  histName = "Merged/CFLLAA3050";
  simpleName = "LLAA3050";
  initParams[0] = radiiParams[2];
  if(useLLAA3050Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLLAA3050, ledInfoLL, initParams, minParams, maxParams, fixParams);


  ////////////////// Setting up Lambda-Lambda ////////////////////
  // 0-10%
  fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/CFs.root";
  histName = "Merged/CFLamLam010";
  simpleName = "LL010";
  initParams[0] = radiiParams[0];
  if(useLL010Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLL010, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 10-30
  histName = "Merged/CFLamLam1030";
  simpleName = "LL1030";
  initParams[0] = radiiParams[1];
  if(useLL1030Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLL1030, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 30-50
  histName = "Merged/CFLamLam3050";
  simpleName = "LL3050";
  initParams[0] = radiiParams[2];
  if(useLL3050Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLL3050, ledInfoLL, initParams, minParams, maxParams, fixParams);


  ////////////////// Setting up Antilambda-Antilambda ////////////////////
  // 0-10%
  fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/CFs.root";
  histName = "Merged/CFALamALam010";
  simpleName = "AA010";
  initParams[0] = radiiParams[0];
  if(useAA010Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kAA010, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 10-30
  histName = "Merged/CFALamALam1030";
  simpleName = "AA1030";
  initParams[0] = radiiParams[1];
  if(useAA1030Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kAA1030, ledInfoLL, initParams, minParams, maxParams, fixParams);

  // 30-50
  histName = "Merged/CFALamALam3050";
  simpleName = "AA3050";
  initParams[0] = radiiParams[2];
  if(useAA3050Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kAA3050, ledInfoLL, initParams, minParams, maxParams, fixParams);




  ////////////////// Setting up Lambda-Antilambda ////////////////////
  // 0-10%
  fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/CFs.root";
  histName = "Merged/CFLamALam010";
  // fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/OneTimeUseTestCF.root";
  // histName = "mm12CF20";

  simpleName = "LA010";
  Double_t initParamsArrLA[5] = {radiiParams[0], -1., 1., 3., 1.}; 
  vector<Double_t> initParamsLA(initParamsArrLA, initParamsArrLA + 5);
  Bool_t fixParamsArrLA[5] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
  vector<Bool_t> fixParamsLA(fixParamsArrLA, fixParamsArrLA + 5);

  vector<LednickyInfo> ledInfoLA = PrepareLednickyInfo(kFALSE, useRootSScalingLA);
  if(useLA010Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLA010, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);

  // 10-30
  histName = "Merged/CFLamALam1030";
  simpleName = "LA1030";
  initParamsLA[0] = radiiParams[1];
  if(useLA1030Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLA1030, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);

  // 30-50
  histName = "Merged/CFLamALam3050";
  simpleName = "LA3050";
  initParamsLA[0] = radiiParams[2];
  if(useLA3050Chi2) fitter->AddPairAnalysisChisquareFit(simpleName, fileName, histName, kLA3050, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);
  //************* Add more systems as needed ******************




  // Log-likelihood fitting.  Do not use in conjunction with
  // chisquare fitting!
  TString fileNumDen = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/CFs.root";
  vector<TString> numNamesLA010;
  numNamesLA010.push_back("mm1/Num/NumLamALam05");
  numNamesLA010.push_back("mm1/Num/NumLamALam510");
  numNamesLA010.push_back("mm2/Num/NumLamALam05");
  numNamesLA010.push_back("mm2/Num/NumLamALam510");
  numNamesLA010.push_back("mm3/Num/NumLamALam05");
  numNamesLA010.push_back("mm3/Num/NumLamALam510");
  numNamesLA010.push_back("pp1/Num/NumLamALam05");
  numNamesLA010.push_back("pp1/Num/NumLamALam510");
  numNamesLA010.push_back("pp2/Num/NumLamALam05");
  numNamesLA010.push_back("pp2/Num/NumLamALam510");
  
  vector<TString> denNamesLA010;
  denNamesLA010.push_back("mm1/Den/DenLamALam05");
  denNamesLA010.push_back("mm1/Den/DenLamALam510");
  denNamesLA010.push_back("mm2/Den/DenLamALam05");
  denNamesLA010.push_back("mm2/Den/DenLamALam510");
  denNamesLA010.push_back("mm3/Den/DenLamALam05");
  denNamesLA010.push_back("mm3/Den/DenLamALam510");
  denNamesLA010.push_back("pp1/Den/DenLamALam05");
  denNamesLA010.push_back("pp1/Den/DenLamALam510");
  denNamesLA010.push_back("pp2/Den/DenLamALam05");
  denNamesLA010.push_back("pp2/Den/DenLamALam510");
  
  initParamsLA[4] = 10.; // Normalization factor accounts for ~1/10 ratio of num vs den pairs
  initParamsLA[0] = radiiParams[0];
  // fitter->AddPairAnalysisLogFit("LA010", fileNumDen, numNamesLA010, denNamesLA010, kLA010, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);

  
  vector<TString> numNamesLA1030;
  numNamesLA1030.push_back("mm1/Num/NumLamALam1015");
  numNamesLA1030.push_back("mm1/Num/NumLamALam1520");
  numNamesLA1030.push_back("mm1/Num/NumLamALam2025");
  numNamesLA1030.push_back("mm1/Num/NumLamALam2530");
  numNamesLA1030.push_back("mm2/Num/NumLamALam1015");
  numNamesLA1030.push_back("mm2/Num/NumLamALam1520");
  numNamesLA1030.push_back("mm2/Num/NumLamALam2025");
  numNamesLA1030.push_back("mm2/Num/NumLamALam2530");
  numNamesLA1030.push_back("mm3/Num/NumLamALam1015");
  numNamesLA1030.push_back("mm3/Num/NumLamALam1520");
  numNamesLA1030.push_back("mm3/Num/NumLamALam2025");
  numNamesLA1030.push_back("mm3/Num/NumLamALam2530");
  numNamesLA1030.push_back("pp1/Num/NumLamALam1015");
  numNamesLA1030.push_back("pp1/Num/NumLamALam1520");
  numNamesLA1030.push_back("pp1/Num/NumLamALam2025");
  numNamesLA1030.push_back("pp1/Num/NumLamALam2530");
  numNamesLA1030.push_back("pp2/Num/NumLamALam1015");
  numNamesLA1030.push_back("pp2/Num/NumLamALam1520");
  numNamesLA1030.push_back("pp2/Num/NumLamALam2025");
  numNamesLA1030.push_back("pp2/Num/NumLamALam2530");
  
  vector<TString> denNamesLA1030;
  denNamesLA1030.push_back("mm1/Den/DenLamALam1015");
  denNamesLA1030.push_back("mm1/Den/DenLamALam1520");
  denNamesLA1030.push_back("mm1/Den/DenLamALam2025");
  denNamesLA1030.push_back("mm1/Den/DenLamALam2530");
  denNamesLA1030.push_back("mm2/Den/DenLamALam1015");
  denNamesLA1030.push_back("mm2/Den/DenLamALam1520");
  denNamesLA1030.push_back("mm2/Den/DenLamALam2025");
  denNamesLA1030.push_back("mm2/Den/DenLamALam2530");
  denNamesLA1030.push_back("mm3/Den/DenLamALam1015");
  denNamesLA1030.push_back("mm3/Den/DenLamALam1520");
  denNamesLA1030.push_back("mm3/Den/DenLamALam2025");
  denNamesLA1030.push_back("mm3/Den/DenLamALam2530");
  denNamesLA1030.push_back("pp1/Den/DenLamALam1015");
  denNamesLA1030.push_back("pp1/Den/DenLamALam1520");
  denNamesLA1030.push_back("pp1/Den/DenLamALam2025");
  denNamesLA1030.push_back("pp1/Den/DenLamALam2530");
  denNamesLA1030.push_back("pp2/Den/DenLamALam1015");
  denNamesLA1030.push_back("pp2/Den/DenLamALam1520");
  denNamesLA1030.push_back("pp2/Den/DenLamALam2025");
  denNamesLA1030.push_back("pp2/Den/DenLamALam2530");
  initParamsLA[0] = radiiParams[1];
  // fitter->AddPairAnalysisLogFit("LA1030", fileNumDen, numNamesLA1030, denNamesLA1030, kLA1030, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);
  
  vector<TString> numNamesLA3050;
  numNamesLA3050.push_back("mm1/Num/NumLamALam3035");
  numNamesLA3050.push_back("mm1/Num/NumLamALam3540");
  numNamesLA3050.push_back("mm1/Num/NumLamALam4045");
  numNamesLA3050.push_back("mm1/Num/NumLamALam4550");
  numNamesLA3050.push_back("mm2/Num/NumLamALam3035");
  numNamesLA3050.push_back("mm2/Num/NumLamALam3540");
  numNamesLA3050.push_back("mm2/Num/NumLamALam4045");
  numNamesLA3050.push_back("mm2/Num/NumLamALam4550");
  numNamesLA3050.push_back("mm3/Num/NumLamALam3035");
  numNamesLA3050.push_back("mm3/Num/NumLamALam3540");
  numNamesLA3050.push_back("mm3/Num/NumLamALam4045");
  numNamesLA3050.push_back("mm3/Num/NumLamALam4550");
  numNamesLA3050.push_back("pp1/Num/NumLamALam3035");
  numNamesLA3050.push_back("pp1/Num/NumLamALam3540");
  numNamesLA3050.push_back("pp1/Num/NumLamALam4045");
  numNamesLA3050.push_back("pp1/Num/NumLamALam4550");
  numNamesLA3050.push_back("pp2/Num/NumLamALam3035");
  numNamesLA3050.push_back("pp2/Num/NumLamALam3540");
  numNamesLA3050.push_back("pp2/Num/NumLamALam4045");
  numNamesLA3050.push_back("pp2/Num/NumLamALam4550");
  
  vector<TString> denNamesLA3050;
  denNamesLA3050.push_back("mm1/Den/DenLamALam3035");
  denNamesLA3050.push_back("mm1/Den/DenLamALam3540");
  denNamesLA3050.push_back("mm1/Den/DenLamALam4045");
  denNamesLA3050.push_back("mm1/Den/DenLamALam4550");
  denNamesLA3050.push_back("mm2/Den/DenLamALam3035");
  denNamesLA3050.push_back("mm2/Den/DenLamALam3540");
  denNamesLA3050.push_back("mm2/Den/DenLamALam4045");
  denNamesLA3050.push_back("mm2/Den/DenLamALam4550");
  denNamesLA3050.push_back("mm3/Den/DenLamALam3035");
  denNamesLA3050.push_back("mm3/Den/DenLamALam3540");
  denNamesLA3050.push_back("mm3/Den/DenLamALam4045");
  denNamesLA3050.push_back("mm3/Den/DenLamALam4550");
  denNamesLA3050.push_back("pp1/Den/DenLamALam3035");
  denNamesLA3050.push_back("pp1/Den/DenLamALam3540");
  denNamesLA3050.push_back("pp1/Den/DenLamALam4045");
  denNamesLA3050.push_back("pp1/Den/DenLamALam4550");
  denNamesLA3050.push_back("pp2/Den/DenLamALam3035");
  denNamesLA3050.push_back("pp2/Den/DenLamALam3540");
  denNamesLA3050.push_back("pp2/Den/DenLamALam4045");
  denNamesLA3050.push_back("pp2/Den/DenLamALam4550");
  initParamsLA[0] = radiiParams[2];
  // fitter->AddPairAnalysisLogFit("LA3050", fileNumDen, numNamesLA3050, denNamesLA3050, kLA3050, ledInfoLA, initParamsLA, minParams, maxParams, fixParamsLA);
















  vector<TString> numNamesLLAA010;
  numNamesLLAA010.push_back("mm1/Num/NumLamLam05");
  numNamesLLAA010.push_back("mm1/Num/NumLamLam510");
  numNamesLLAA010.push_back("mm2/Num/NumLamLam05");
  numNamesLLAA010.push_back("mm2/Num/NumLamLam510");
  numNamesLLAA010.push_back("mm3/Num/NumLamLam05");
  numNamesLLAA010.push_back("mm3/Num/NumLamLam510");
  numNamesLLAA010.push_back("pp1/Num/NumLamLam05");
  numNamesLLAA010.push_back("pp1/Num/NumLamLam510");
  numNamesLLAA010.push_back("pp2/Num/NumLamLam05");
  numNamesLLAA010.push_back("pp2/Num/NumLamLam510");

  numNamesLLAA010.push_back("mm1/Num/NumALamALam05");
  numNamesLLAA010.push_back("mm1/Num/NumALamALam510");
  numNamesLLAA010.push_back("mm2/Num/NumALamALam05");
  numNamesLLAA010.push_back("mm2/Num/NumALamALam510");
  numNamesLLAA010.push_back("mm3/Num/NumALamALam05");
  numNamesLLAA010.push_back("mm3/Num/NumALamALam510");
  numNamesLLAA010.push_back("pp1/Num/NumALamALam05");
  numNamesLLAA010.push_back("pp1/Num/NumALamALam510");
  numNamesLLAA010.push_back("pp2/Num/NumALamALam05");
  numNamesLLAA010.push_back("pp2/Num/NumALamALam510");

  vector<TString> denNamesLLAA010;
  denNamesLLAA010.push_back("mm1/Den/DenLamLam05");
  denNamesLLAA010.push_back("mm1/Den/DenLamLam510");
  denNamesLLAA010.push_back("mm2/Den/DenLamLam05");
  denNamesLLAA010.push_back("mm2/Den/DenLamLam510");
  denNamesLLAA010.push_back("mm3/Den/DenLamLam05");
  denNamesLLAA010.push_back("mm3/Den/DenLamLam510");
  denNamesLLAA010.push_back("pp1/Den/DenLamLam05");
  denNamesLLAA010.push_back("pp1/Den/DenLamLam510");
  denNamesLLAA010.push_back("pp2/Den/DenLamLam05");
  denNamesLLAA010.push_back("pp2/Den/DenLamLam510");

  denNamesLLAA010.push_back("mm1/Den/DenALamALam05");
  denNamesLLAA010.push_back("mm1/Den/DenALamALam510");
  denNamesLLAA010.push_back("mm2/Den/DenALamALam05");
  denNamesLLAA010.push_back("mm2/Den/DenALamALam510");
  denNamesLLAA010.push_back("mm3/Den/DenALamALam05");
  denNamesLLAA010.push_back("mm3/Den/DenALamALam510");
  denNamesLLAA010.push_back("pp1/Den/DenALamALam05");
  denNamesLLAA010.push_back("pp1/Den/DenALamALam510");
  denNamesLLAA010.push_back("pp2/Den/DenALamALam05");
  denNamesLLAA010.push_back("pp2/Den/DenALamALam510");





  
  initParams[4] = 10.; // Normalization factor accounts for ~1/10 ratio of num vs den pairs
  initParams[0] = radiiParams[0];
  // fitter->AddPairAnalysisLogFit("LLAA010", fileNumDen, numNamesLLAA010, denNamesLLAA010, kLLAA010, ledInfoLL, initParams, minParams, maxParams, fixParams);




  
  vector<TString> numNamesLLAA1030;
  numNamesLLAA1030.push_back("mm1/Num/NumLamLam1015");
  numNamesLLAA1030.push_back("mm1/Num/NumLamLam1520");
  numNamesLLAA1030.push_back("mm1/Num/NumLamLam2025");
  numNamesLLAA1030.push_back("mm1/Num/NumLamLam2530");
  numNamesLLAA1030.push_back("mm2/Num/NumLamLam1015");
  numNamesLLAA1030.push_back("mm2/Num/NumLamLam1520");
  numNamesLLAA1030.push_back("mm2/Num/NumLamLam2025");
  numNamesLLAA1030.push_back("mm2/Num/NumLamLam2530");
  numNamesLLAA1030.push_back("mm3/Num/NumLamLam1015");
  numNamesLLAA1030.push_back("mm3/Num/NumLamLam1520");
  numNamesLLAA1030.push_back("mm3/Num/NumLamLam2025");
  numNamesLLAA1030.push_back("mm3/Num/NumLamLam2530");
  numNamesLLAA1030.push_back("pp1/Num/NumLamLam1015");
  numNamesLLAA1030.push_back("pp1/Num/NumLamLam1520");
  numNamesLLAA1030.push_back("pp1/Num/NumLamLam2025");
  numNamesLLAA1030.push_back("pp1/Num/NumLamLam2530");
  numNamesLLAA1030.push_back("pp2/Num/NumLamLam1015");
  numNamesLLAA1030.push_back("pp2/Num/NumLamLam1520");
  numNamesLLAA1030.push_back("pp2/Num/NumLamLam2025");
  numNamesLLAA1030.push_back("pp2/Num/NumLamLam2530");
  
  numNamesLLAA1030.push_back("mm1/Num/NumALamALam1015");
  numNamesLLAA1030.push_back("mm1/Num/NumALamALam1520");
  numNamesLLAA1030.push_back("mm1/Num/NumALamALam2025");
  numNamesLLAA1030.push_back("mm1/Num/NumALamALam2530");
  numNamesLLAA1030.push_back("mm2/Num/NumALamALam1015");
  numNamesLLAA1030.push_back("mm2/Num/NumALamALam1520");
  numNamesLLAA1030.push_back("mm2/Num/NumALamALam2025");
  numNamesLLAA1030.push_back("mm2/Num/NumALamALam2530");
  numNamesLLAA1030.push_back("mm3/Num/NumALamALam1015");
  numNamesLLAA1030.push_back("mm3/Num/NumALamALam1520");
  numNamesLLAA1030.push_back("mm3/Num/NumALamALam2025");
  numNamesLLAA1030.push_back("mm3/Num/NumALamALam2530");
  numNamesLLAA1030.push_back("pp1/Num/NumALamALam1015");
  numNamesLLAA1030.push_back("pp1/Num/NumALamALam1520");
  numNamesLLAA1030.push_back("pp1/Num/NumALamALam2025");
  numNamesLLAA1030.push_back("pp1/Num/NumALamALam2530");
  numNamesLLAA1030.push_back("pp2/Num/NumALamALam1015");
  numNamesLLAA1030.push_back("pp2/Num/NumALamALam1520");
  numNamesLLAA1030.push_back("pp2/Num/NumALamALam2025");
  numNamesLLAA1030.push_back("pp2/Num/NumALamALam2530");
  
  vector<TString> denNamesLLAA1030;
  denNamesLLAA1030.push_back("mm1/Den/DenLamLam1015");
  denNamesLLAA1030.push_back("mm1/Den/DenLamLam1520");
  denNamesLLAA1030.push_back("mm1/Den/DenLamLam2025");
  denNamesLLAA1030.push_back("mm1/Den/DenLamLam2530");
  denNamesLLAA1030.push_back("mm2/Den/DenLamLam1015");
  denNamesLLAA1030.push_back("mm2/Den/DenLamLam1520");
  denNamesLLAA1030.push_back("mm2/Den/DenLamLam2025");
  denNamesLLAA1030.push_back("mm2/Den/DenLamLam2530");
  denNamesLLAA1030.push_back("mm3/Den/DenLamLam1015");
  denNamesLLAA1030.push_back("mm3/Den/DenLamLam1520");
  denNamesLLAA1030.push_back("mm3/Den/DenLamLam2025");
  denNamesLLAA1030.push_back("mm3/Den/DenLamLam2530");
  denNamesLLAA1030.push_back("pp1/Den/DenLamLam1015");
  denNamesLLAA1030.push_back("pp1/Den/DenLamLam1520");
  denNamesLLAA1030.push_back("pp1/Den/DenLamLam2025");
  denNamesLLAA1030.push_back("pp1/Den/DenLamLam2530");
  denNamesLLAA1030.push_back("pp2/Den/DenLamLam1015");
  denNamesLLAA1030.push_back("pp2/Den/DenLamLam1520");
  denNamesLLAA1030.push_back("pp2/Den/DenLamLam2025");
  denNamesLLAA1030.push_back("pp2/Den/DenLamLam2530");

  denNamesLLAA1030.push_back("mm1/Den/DenALamALam1015");
  denNamesLLAA1030.push_back("mm1/Den/DenALamALam1520");
  denNamesLLAA1030.push_back("mm1/Den/DenALamALam2025");
  denNamesLLAA1030.push_back("mm1/Den/DenALamALam2530");
  denNamesLLAA1030.push_back("mm2/Den/DenALamALam1015");
  denNamesLLAA1030.push_back("mm2/Den/DenALamALam1520");
  denNamesLLAA1030.push_back("mm2/Den/DenALamALam2025");
  denNamesLLAA1030.push_back("mm2/Den/DenALamALam2530");
  denNamesLLAA1030.push_back("mm3/Den/DenALamALam1015");
  denNamesLLAA1030.push_back("mm3/Den/DenALamALam1520");
  denNamesLLAA1030.push_back("mm3/Den/DenALamALam2025");
  denNamesLLAA1030.push_back("mm3/Den/DenALamALam2530");
  denNamesLLAA1030.push_back("pp1/Den/DenALamALam1015");
  denNamesLLAA1030.push_back("pp1/Den/DenALamALam1520");
  denNamesLLAA1030.push_back("pp1/Den/DenALamALam2025");
  denNamesLLAA1030.push_back("pp1/Den/DenALamALam2530");
  denNamesLLAA1030.push_back("pp2/Den/DenALamALam1015");
  denNamesLLAA1030.push_back("pp2/Den/DenALamALam1520");
  denNamesLLAA1030.push_back("pp2/Den/DenALamALam2025");
  denNamesLLAA1030.push_back("pp2/Den/DenALamALam2530");


  
  initParams[0] = radiiParams[1];
  // fitter->AddPairAnalysisLogFit("LLAA1030", fileNumDen, numNamesLLAA1030, denNamesLLAA1030, kLLAA1030, ledInfoLL, initParams, minParams, maxParams, fixParams);




  vector<TString> numNamesLLAA3050;
  numNamesLLAA3050.push_back("mm1/Num/NumLamLam3035");
  numNamesLLAA3050.push_back("mm1/Num/NumLamLam3540");
  numNamesLLAA3050.push_back("mm1/Num/NumLamLam4045");
  numNamesLLAA3050.push_back("mm1/Num/NumLamLam4550");
  numNamesLLAA3050.push_back("mm2/Num/NumLamLam3035");
  numNamesLLAA3050.push_back("mm2/Num/NumLamLam3540");
  numNamesLLAA3050.push_back("mm2/Num/NumLamLam4045");
  numNamesLLAA3050.push_back("mm2/Num/NumLamLam4550");
  numNamesLLAA3050.push_back("mm3/Num/NumLamLam3035");
  numNamesLLAA3050.push_back("mm3/Num/NumLamLam3540");
  numNamesLLAA3050.push_back("mm3/Num/NumLamLam4045");
  numNamesLLAA3050.push_back("mm3/Num/NumLamLam4550");
  numNamesLLAA3050.push_back("pp1/Num/NumLamLam3035");
  numNamesLLAA3050.push_back("pp1/Num/NumLamLam3540");
  numNamesLLAA3050.push_back("pp1/Num/NumLamLam4045");
  numNamesLLAA3050.push_back("pp1/Num/NumLamLam4550");
  numNamesLLAA3050.push_back("pp2/Num/NumLamLam3035");
  numNamesLLAA3050.push_back("pp2/Num/NumLamLam3540");
  numNamesLLAA3050.push_back("pp2/Num/NumLamLam4045");
  numNamesLLAA3050.push_back("pp2/Num/NumLamLam4550");
  
  numNamesLLAA3050.push_back("mm1/Num/NumALamALam3035");
  numNamesLLAA3050.push_back("mm1/Num/NumALamALam3540");
  numNamesLLAA3050.push_back("mm1/Num/NumALamALam4045");
  numNamesLLAA3050.push_back("mm1/Num/NumALamALam4550");
  numNamesLLAA3050.push_back("mm2/Num/NumALamALam3035");
  numNamesLLAA3050.push_back("mm2/Num/NumALamALam3540");
  numNamesLLAA3050.push_back("mm2/Num/NumALamALam4045");
  numNamesLLAA3050.push_back("mm2/Num/NumALamALam4550");
  numNamesLLAA3050.push_back("mm3/Num/NumALamALam3035");
  numNamesLLAA3050.push_back("mm3/Num/NumALamALam3540");
  numNamesLLAA3050.push_back("mm3/Num/NumALamALam4045");
  numNamesLLAA3050.push_back("mm3/Num/NumALamALam4550");
  numNamesLLAA3050.push_back("pp1/Num/NumALamALam3035");
  numNamesLLAA3050.push_back("pp1/Num/NumALamALam3540");
  numNamesLLAA3050.push_back("pp1/Num/NumALamALam4045");
  numNamesLLAA3050.push_back("pp1/Num/NumALamALam4550");
  numNamesLLAA3050.push_back("pp2/Num/NumALamALam3035");
  numNamesLLAA3050.push_back("pp2/Num/NumALamALam3540");
  numNamesLLAA3050.push_back("pp2/Num/NumALamALam4045");
  numNamesLLAA3050.push_back("pp2/Num/NumALamALam4550");
  
  vector<TString> denNamesLLAA3050;
  denNamesLLAA3050.push_back("mm1/Den/DenLamLam3035");
  denNamesLLAA3050.push_back("mm1/Den/DenLamLam3540");
  denNamesLLAA3050.push_back("mm1/Den/DenLamLam4045");
  denNamesLLAA3050.push_back("mm1/Den/DenLamLam4550");
  denNamesLLAA3050.push_back("mm2/Den/DenLamLam3035");
  denNamesLLAA3050.push_back("mm2/Den/DenLamLam3540");
  denNamesLLAA3050.push_back("mm2/Den/DenLamLam4045");
  denNamesLLAA3050.push_back("mm2/Den/DenLamLam4550");
  denNamesLLAA3050.push_back("mm3/Den/DenLamLam3035");
  denNamesLLAA3050.push_back("mm3/Den/DenLamLam3540");
  denNamesLLAA3050.push_back("mm3/Den/DenLamLam4045");
  denNamesLLAA3050.push_back("mm3/Den/DenLamLam4550");
  denNamesLLAA3050.push_back("pp1/Den/DenLamLam3035");
  denNamesLLAA3050.push_back("pp1/Den/DenLamLam3540");
  denNamesLLAA3050.push_back("pp1/Den/DenLamLam4045");
  denNamesLLAA3050.push_back("pp1/Den/DenLamLam4550");
  denNamesLLAA3050.push_back("pp2/Den/DenLamLam3035");
  denNamesLLAA3050.push_back("pp2/Den/DenLamLam3540");
  denNamesLLAA3050.push_back("pp2/Den/DenLamLam4045");
  denNamesLLAA3050.push_back("pp2/Den/DenLamLam4550");
  
  denNamesLLAA3050.push_back("mm1/Den/DenALamALam3035");
  denNamesLLAA3050.push_back("mm1/Den/DenALamALam3540");
  denNamesLLAA3050.push_back("mm1/Den/DenALamALam4045");
  denNamesLLAA3050.push_back("mm1/Den/DenALamALam4550");
  denNamesLLAA3050.push_back("mm2/Den/DenALamALam3035");
  denNamesLLAA3050.push_back("mm2/Den/DenALamALam3540");
  denNamesLLAA3050.push_back("mm2/Den/DenALamALam4045");
  denNamesLLAA3050.push_back("mm2/Den/DenALamALam4550");
  denNamesLLAA3050.push_back("mm3/Den/DenALamALam3035");
  denNamesLLAA3050.push_back("mm3/Den/DenALamALam3540");
  denNamesLLAA3050.push_back("mm3/Den/DenALamALam4045");
  denNamesLLAA3050.push_back("mm3/Den/DenALamALam4550");
  denNamesLLAA3050.push_back("pp1/Den/DenALamALam3035");
  denNamesLLAA3050.push_back("pp1/Den/DenALamALam3540");
  denNamesLLAA3050.push_back("pp1/Den/DenALamALam4045");
  denNamesLLAA3050.push_back("pp1/Den/DenALamALam4550");
  denNamesLLAA3050.push_back("pp2/Den/DenALamALam3035");
  denNamesLLAA3050.push_back("pp2/Den/DenALamALam3540");
  denNamesLLAA3050.push_back("pp2/Den/DenALamALam4045");
  denNamesLLAA3050.push_back("pp2/Den/DenALamALam4550");
  



  
  initParams[0] = radiiParams[2];
  // fitter->AddPairAnalysisLogFit("LLAA3050", fileNumDen, numNamesLLAA3050, denNamesLLAA3050, kLLAA3050, ledInfoLL, initParams, minParams, maxParams, fixParams);

  

  
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
  // myFitter->SetupConstraint(kF0Real, systemsLLAA);
  // myFitter->SetupConstraint(kD0, systemsLLAA);

  // Share real f0, imaginary f0, and d0 among LambdaLambda
  Int_t systemsArrLL[3] = {kLL010, kLL1030, kLL3050};
  vector<Int_t> systemsLL(systemsArrLL, systemsArrLL + 3);
  // myFitter->SetupConstraint(kF0Real, systemsLL);
  // myFitter->SetupConstraint(kD0, systemsLL);

  // Share real f0, imaginary f0, and d0 among LambdaLambda
  Int_t systemsArrAA[3] = {kAA010, kAA1030, kAA3050};
  vector<Int_t> systemsAA(systemsArrAA, systemsArrAA + 3);
  // myFitter->SetupConstraint(kF0Real, systemsAA);
  // myFitter->SetupConstraint(kD0, systemsAA);


  // Constrain LL and AA together
  Int_t systemsArrLLandAA[6] = {kLL010, kLL1030,kLL3050,
				kAA010, kAA1030,kAA3050};
  vector<Int_t> systemsLLandAA(systemsArrLLandAA, systemsArrLLandAA + 6);
  myFitter->SetupConstraint(kF0Real, systemsLLandAA);
  myFitter->SetupConstraint(kD0, systemsLLandAA);

  // Constrain Averaged LLAA with LA
  Int_t systemsArrLLAALA[6] = {kLLAA010, kLLAA1030, kLLAA3050,
			     kLA010, kLA1030, kLA3050};
  vector<Int_t> systemsLLAALA(systemsArrLLAALA, systemsArrLLAALA + 6);
  // myFitter->SetupConstraint(kF0Real, systemsLLAALA);
  // myFitter->SetupConstraint(kD0, systemsLLAALA);


  // share radii among like centralities
  vector<Int_t> systems010;
  // systems010.push_back(kLLAA010);
  systems010.push_back(kLA010);
  systems010.push_back(kLL010);
  systems010.push_back(kAA010);
  myFitter->SetupConstraint(kRad, systems010);

  vector<Int_t> systems1030;
  // systems1030.push_back(kLLAA1030);
  systems1030.push_back(kLA1030);
  systems1030.push_back(kLL1030);
  systems1030.push_back(kAA1030);
  myFitter->SetupConstraint(kRad, systems1030);

  vector<Int_t> systems3050;
  // systems3050.push_back(kLLAA3050);
  systems3050.push_back(kLA3050);
  systems3050.push_back(kLL3050);
  systems3050.push_back(kAA3050);
  myFitter->SetupConstraint(kRad, systems3050);
}

void UserSetFitOptions(Fitter *myFitter)
{
  // Set the upper bin of the fit range
   myFitter->SetHighFitBin(40);
   
  // How big should the initial starting value fit uncertainty be?
  // This is a catchall for *all* parameter step sizes.
   myFitter->SetStartingStepSize(.1);


   // Use MINOS to find error bars?
   myFitter->SetUseMINOS(kFALSE);
   
  // optional suffix for saved plots and objects
  TString outString = "";
  myFitter->SetOutputString(outString);

  // Output extra plots showing all residual correlation components?
  myFitter->SetDisplayResidualComponents(kTRUE);
  
  // Output extra histograms of lednicky Eqns?
  myFitter->SetOutputLednickyHists(kTRUE);
  
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

  // Get the parameters set up in Fitter
  myFitter->SetupInitialParameterVectors();
  
  // Create a TMinuit object.
  // Define, set and fix parameters.
  const Int_t nParams = myFitter->GetNMinuitParams();
  cout<<"Number of Minuit Parameters:\t"<<nParams<<endl;
  TMinuit *myMinuit = new TMinuit(nParams);
  myMinuit->SetFCN(MinuitFit);
  myFitter->InitializeMinuitParameters(myMinuit);
  // myFitter->CreateMinuit();
  myFitter->DoFitting();
  myFitter->SaveOutputPlots();

  delete myMinuit;
  return 0;
}
 
