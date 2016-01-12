//**********************************************************
// Fitting correlation functions
//**********************************************************

#include "Fitter.h"
#include "TH1D.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <iostream>
#include <assert.h>
#include <cstdlib>



using namespace std;


Fitter::Fitter():
  fNParamsTotal(0),
  fFixedParams(0),
  fFitBins(40),
  fHighFitBin(40),
  fLowFitBin(0),
  fNSystems(0),
  fNLogLikelihoodDataSets(0),
  fMaxMinuitCalls(20000),
  fStepSize(0.1),
  fUseEstimatedLambdaParams(kTRUE),
  fUseMINOS(kFALSE),
  fDisplayResidualComponents(kFALSE),
  fUseChisquareFitting(kFALSE),
  fUseLogLikelihoodFitting(kFALSE),
  fOutputString(""),
  fFitCalls(0),
  fChisquare(0.),
  fMinuitVerbosity(0)
{
  fTimer = new TStopwatch();
  fParamNames.push_back("Radius");
  fParamNames.push_back("ReF0");
  fParamNames.push_back("ImF0");
  fParamNames.push_back("D0");
  fParamNames.push_back("Norm");
  
}

Fitter::~Fitter()
{
  for(UInt_t i = 0; i < fPairSystems.size(); i++)
  {
    if(!fPairSystems[i]) continue;
    delete fPairSystems[i];
    fPairSystems[i] = NULL;
  }
  for(UInt_t i = 0; i < fParamConstraints.size(); i++)
  {
    if(!fParamConstraints[i]) continue;
    delete fParamConstraints[i];
    fParamConstraints[i] = NULL;
  }
  // if(fMinuit) {delete fMinuit; fMinuit = NULL;}
  if(fTimer) {delete fTimer; fTimer = NULL;}
}

void Fitter::AddPairAnalysisChisquareFit(TString simpleName, TString fileName, TString cfName, Int_t sysType, const vector<LednickyInfo> &ledInfo, vector<Double_t> initParams, vector<Double_t> minParams, vector<Double_t> maxParams, vector<Bool_t> fixParams)
{
  // Add a fit analysis that will do a chisquare fit of a
  // correlation function.
  if(fUseLogLikelihoodFitting) {
    cerr<<"ERROR: Attempt to add chisquare analysis after a "
	<<"log-likelihood fit analysis was already included."
	<<endl;
    std::exit(EXIT_FAILURE);
  }
  fUseChisquareFitting = kTRUE;
  CreatePairSystemChisquare(simpleName, fileName, cfName, sysType, ledInfo);
  UInt_t numberOfNormParams = 1;
  PushBackParams(simpleName, initParams, minParams, maxParams, fixParams, numberOfNormParams);
}

void Fitter::AddPairAnalysisLogFit(TString simpleName, TString fileName, vector<TString> numNames, vector<TString> denNames, Int_t sysType, const vector<LednickyInfo> &ledInfo, vector<Double_t> initParams, vector<Double_t> minParams, vector<Double_t> maxParams, vector<Bool_t> fixParams)
{
  // Add a fit analysis that will do a log-likelihood fit
  // directly to numerator and denominator distributions
  if(fUseChisquareFitting) {
    cerr<<"ERROR: Attempt to add log-likelihood fit analysis after "
	<<"a chisquare fit analysis was already included."
	<<endl;
    std::exit(EXIT_FAILURE);
  }
  fUseLogLikelihoodFitting = kTRUE;
  CreatePairSystemLog(simpleName, fileName, numNames, denNames, sysType, ledInfo);
  UInt_t numberOfNormParams = numNames.size();
  fNLogLikelihoodDataSets += numberOfNormParams;
  PushBackParams(simpleName, initParams, minParams, maxParams, fixParams, numberOfNormParams);
}


void Fitter::CreatePairSystemChisquare(TString simpleName, TString fileName, TString cfName, Int_t sysType, const vector<LednickyInfo> &ledInfo)
{
  // Create a pair system for doing a chisquare fit to a
  // correlation function.
  TFile inFile(fileName, "read");
  TH1D *cf = (TH1D*) inFile.Get(cfName);
  assert(cf);
  cf->SetDirectory(0);
  PairSystem *system = new PairSystem(cf, ledInfo, simpleName, sysType);
  fPairSystems.push_back(system);
}

void Fitter::CreatePairSystemLog(TString simpleName, TString fileName, vector<TString> numNames, vector<TString> denNames, Int_t sysType, const vector<LednickyInfo> &ledInfo)
{
  // Create a pair system for a log-likelihood fit
  // of numerator and denominator distributions
  TFile inFile(fileName, "read");
  assert(numNames.size());
  assert(numNames.size() == denNames.size());
  vector<TH1D*> numHists;
  vector<TH1D*> denHists;
  for(UInt_t i = 0; i < numNames.size(); i++) {
    TH1D *num = (TH1D*) inFile.Get(numNames[i]);
    TH1D *den = (TH1D*) inFile.Get(denNames[i]);
    assert(num && den);
    num->SetDirectory(0);
    den->SetDirectory(0);
    numHists.push_back(num);
    denHists.push_back(den);
  }
  PairSystem *system = new PairSystem(numHists, denHists, ledInfo, simpleName, sysType);
  fPairSystems.push_back(system);
  cout<<"Added pair system with "<<numNames.size()<<" num and den pairs"<<endl;
}






void Fitter::DoFitting()
{
  cout<<"DoFitting"<<endl;
  assert(fMinuit); // fMinuit needs to have been created by now
  Timer();
  // Run the fit procedure

  
  Double_t arglist[5] = {0,0,0,0,0}; //Arguments that can be passed with Minuit commands
  Int_t errFlag = 0;
  arglist[0] = 1;
  fMinuit->mnexcm("CALL FCN", arglist, 1, errFlag);

  
  // Set how verbose the output is (from no output at -1, to max at 3)
  arglist[0] = fMinuitVerbosity;
  fMinuit->mnexcm("SET PRINT", arglist, 1, errFlag);

  Timer();
  // Run Migrad with a specified max number of calls
  arglist[0] = fMaxMinuitCalls;
  arglist[1] = 0.1;
  fMinuit->mnexcm("MIGRAD", arglist, 1, errFlag);
  Timer();
  
  fMinuit->mnexcm("SHOw CORrelations", arglist, 1, errFlag);

  if(fUseMINOS) {
    fMinuit->mnexcm("MINOS", arglist, 1, errFlag);
    fMinuit->mnexcm("SHOw CORrelations", arglist, 1, errFlag);
  }

  cout<<"Finalchi2:\t"<<fChisquare<<endl
      <<"Fit bins:\t"<<fFitBins<<endl
      <<"Actual Minuit Pars:\t"<<fMinuitParNames.size()<<endl
      <<"Fixed params\t"<<fFixedParams<<endl
      <<"Free Params:\t"<<fMinuitParNames.size() - fFixedParams<<endl
      <<"Chi2/ndf:\t"<<GetChisquarePerNDF()<<endl
      <<"Pvalue:\t"<<GetPvalue()<<endl;

}


Int_t Fitter::GetConstrainedParamIndex(const Int_t currentSys, const Int_t currentPar)
{
  // Find index of the earlier matching constrained parameter
  
  Int_t sysType = fPairSystems[currentSys]->GetSystemType();

  // Loop through the constraints until we find the relevant one
  for(UInt_t iCon = 0; iCon < fParamConstraints.size(); iCon++)
  {
    // Check for constraints on this type of parameter
    ParameterConstraint *constraint = fParamConstraints[iCon];
    const Int_t parIndex = constraint->GetConstrainedParam();
    if(parIndex != currentPar) continue;
    
    const vector<Int_t> &consSystems = fParamConstraints[iCon]->GetConstrainedSystems();

    for(UInt_t iSys = 1; iSys < consSystems.size(); iSys++)
    {
      
      // Does this system have this constraint?
      if(consSystems[iSys] == sysType) 
      {
	// Get the type of the first system in the constraint
	const Int_t sysTypePrior = consSystems[0];
	// Now find the index of the system with that type.
	Int_t sysIndexPrior = -1;
	for(Int_t iPairSys = 0; iPairSys < fNSystems; iPairSys++)
	{
	  if(sysTypePrior == fPairSystems[iPairSys]->GetSystemType())
	  {
	    sysIndexPrior = iPairSys;
	  }
	}
	assert(sysIndexPrior != -1);
	
	// Return the Minuit index of the parameter.
	Int_t absoluteIndex = 0;
	for(Int_t iSysInner = 0; iSysInner < sysIndexPrior; iSysInner++) {
	  absoluteIndex += fNParamsPerSystem[iSysInner];
	}
	
        absoluteIndex += parIndex;
	return absoluteIndex;
      }
    } 
  }
  // Should not get to this point
  return 100000;
}

Double_t Fitter::GetChisquarePerNDF()
{
  Double_t fitBins = fFitBins;
  if(fUseLogLikelihoodFitting) {
    fitBins *= fNLogLikelihoodDataSets;
  }
  Double_t ndf = 1.*fitBins - (1.*fMinuitParNames.size() - 1.*fFixedParams);
  return fChisquare/ndf;
}

Double_t Fitter::GetPvalue()
{
  Double_t ndf = 1.*fFitBins - (1.*fMinuitParNames.size() - 1.*fFixedParams);
  return TMath::Prob(fChisquare,ndf); 
}

void Fitter::InitializeMinuitParameters(TMinuit *minuit)
{
  // Define the fit parameters
  fMinuit = minuit;
  Double_t startingStepSize = fStepSize;
  for(UInt_t iPar = 0; iPar < fMinuitParNames.size(); iPar++){
    fMinuit->DefineParameter(iPar, fMinuitParNames[iPar], fMinuitParInitial[iPar],  startingStepSize, fMinuitParMinimum[iPar], fMinuitParMaximum[iPar]);
    if(fMinuitParIsFixed[iPar]) fMinuit->FixParameter(iPar);
  }
  
}

Bool_t Fitter::IsParameterConstrained(const Int_t currentSys, const Int_t currentPar)
{
  // Check to see if this parameter has been constrained 
  // to be the same as the parameter from an earlier system
  // cout<<"IsParameterConstrained\n";
  Int_t sysType = fPairSystems[currentSys]->GetSystemType();
  for(UInt_t iCon = 0; iCon < fParamConstraints.size(); iCon++)
  {
    // Check for constraints on this type of parameter
    ParameterConstraint *constraint = fParamConstraints[iCon];
    if(constraint->GetConstrainedParam() != currentPar) continue;
    const vector<Int_t> &consSystems = fParamConstraints[iCon]->GetConstrainedSystems();
    for(UInt_t iSys = 1; iSys < consSystems.size(); iSys++)
    {
      if(consSystems[iSys] == sysType) {
	return true;
      }
    } 
  }
  return kFALSE;
}

void Fitter::PushBackParams(TString simpleName, vector<Double_t> initParams, vector<Double_t> minParams, vector<Double_t> maxParams, vector<Bool_t> fixParams, UInt_t nNormParams)
{
  // Add all the parameters to the fitter

  // If we need extra normalization parameters (for log fitting),
  // add them here. Add extra copies of the supplied norm param.
  for(UInt_t i = 1; i < nNormParams; i++) {
    initParams.push_back(initParams[4]);
    minParams.push_back(minParams[4]);
    maxParams.push_back(maxParams[4]);
    fixParams.push_back(fixParams[4]);
  }
  // Add these parameters to the fitter's list
  fSystemNames.push_back(simpleName);
  fInitParams.push_back(initParams);
  fMinParams.push_back(minParams);
  fMaxParams.push_back(maxParams);
  fFixParams.push_back(fixParams);
  Int_t nParams = initParams.size();
  fNParamsPerSystem.push_back(nParams);
  fNParamsTotal += nParams;
  fNSystems++;
  for(UInt_t i = 0; i < fixParams.size(); i++)
  {
    if(fixParams[i]) fFixedParams++;
  }
}

void Fitter::SaveOutputPlots()
{
  // Save output plots for each fit
  TFile outFile("FitResults.root","update");
  TDirectory *dirLedHists = outFile.GetDirectory("LednickyHists");
  if(fOutputLednickyHist) {
    if(!dirLedHists) {
      outFile.mkdir("LednickyHists");
      dirLedHists = outFile.GetDirectory("LednickyHists");
    }
  }
  // TStyle myStyle;
  gStyle->SetOptStat(0);
  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    outFile.cd();
    TH1D *cf = fPairSystems[iSys]->GetCF();
    if(cf) {
      cf->SetAxisRange(0.7,1.05,"Y");
      cf->SetAxisRange(0.,1.,"X");
      cf->SetMarkerStyle(20);
      cf->SetMarkerColor(1);
      cf->SetLineColor(1);
      cf->GetXaxis()->SetLabelSize(0.06);
      cf->GetYaxis()->SetLabelSize(0.06);
      cf->GetXaxis()->SetTitleSize(0.06);
      cf->GetYaxis()->SetTitleSize(0.06);
      cf->GetXaxis()->SetNdivisions(505);
      cf->GetYaxis()->SetNdivisions(505);
      cf->GetXaxis()->SetTitleOffset(0.8);
      cf->GetXaxis()->CenterTitle();

      TString cfName = cf->GetName();
      // cfName += fOutputString;
      cf->Write(cfName, TObject::kOverwrite);
    }
    TGraph *g = fPairSystems[iSys]->GetCombinedTGraph();
    TString gName = g->GetName();
    gName += fOutputString;
    g->Write(gName, TObject::kOverwrite);
    TCanvas c1("c1","c1");
    if(cf) {
      cf->DrawCopy();
      g->Draw("same");
    }
    else g->Draw();
    TString plotName = "Plot";
    plotName += gName;
    c1.SaveAs(plotName + ".pdf");
    c1.SaveAs(plotName + ".png");
    c1.SaveAs(plotName + ".eps");
    cout<<"Saved file "<<plotName<<endl;
    if(fDisplayResidualComponents) SaveResidualComponentPlot(iSys);
    if(fOutputLednickyHist) {
      dirLedHists->cd();
      TH1D *ledHist = fPairSystems[iSys]->GetLednickyHist();
      // TString ledHistName = ledHist->GetName();
      ledHist->Write(ledHist->GetName()+fOutputString, TObject::kOverwrite);
    }
  }

  
}

void Fitter::SaveResidualComponentPlot(Int_t sys)
{
  // Save output plot showing all primary and residual correlation
  // components of the correlation function.
  TCanvas *components = fPairSystems[sys]->GetResidualComponentCanvas();
  if(!components) {
    cout<<"No Residual Components plot found\n";
    return;
  }

  TH1D *cf = fPairSystems[sys]->GetCF();
  TGraph *g = fPairSystems[sys]->GetCombinedTGraph();
  if(cf) cf->DrawCopy("same");
  if(g) g->Draw("same");
  
  TString plotName = "Plot";
  plotName += components->GetName();
  plotName += fOutputString;
  components->SaveAs(plotName + "Residuals.pdf");
  components->SaveAs(plotName + "Residuals.png");
  components->SaveAs(plotName + "Residuals.eps");
  cout<<"Saved file "<<plotName<<endl;
  delete components; components = NULL;
}

void Fitter::SetHighFitBin(Int_t bin)
{
  for(UInt_t iSys = 0; iSys < fNSystems; iSys++)
  {
    fPairSystems[iSys]->SetHighFitBin(bin);
  }
  fHighFitBin = bin;
  fFitBins = (fHighFitBin - fLowFitBin) * fNSystems;
}

void Fitter::SetLowFitBin(Int_t bin)
{
  for(UInt_t iSys = 0; iSys < fNSystems; iSys++)
  {
    fPairSystems[iSys]->SetLowFitBin(bin);
  }
  fLowFitBin = bin;
  fFitBins = (fHighFitBin - fLowFitBin) * fNSystems;
}

void Fitter::SetMinuitVerbosity(Int_t verbosity)
{
  // Set the verbosity level of minuit fitting
  // -1, 0, 1, 2, 3 (-1 for minimal output, 3 for max output)
  fMinuitVerbosity = verbosity;
  if(fMinuit) {
    Double_t arglist[1] = {(Double_t)fMinuitVerbosity};
    Int_t errFlag = 0;
    fMinuit->mnexcm("SET PRINT", arglist, 1, errFlag);
  }
  
}

void Fitter::SetParametersAndFit(Int_t& i, Double_t &totalChisquare, Double_t *par)
{
  // cout<<"SetParametersAndFitBegin:\t"<<++fFitCalls<<endl;
  Int_t nCalls = fMinuit->fNfcn;
  if(nCalls % 100 == 0) {
    cout<<"Begin fit iteration:\t"<<nCalls<<".\t"
	<<"Last Chi^2:\t"<<fChisquare<<".\t";
    Timer();
  }
  // Take the TMinuit parameters, set the parameters for each
  // pair system, and get the resulting chisquare of the fit.

  // We'll copy par into the parameters vector, and insert 
  // any constrained parameters into their appropriate
  // positions
  vector<Double_t> parameters(fNParamsTotal);
  Int_t constrainedParams = 0;
  Int_t currentParamIndex = -1; //First param will increment to 0
  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    Int_t nParams = fNParamsPerSystem[iSys];
    for(Int_t iPar = 0; iPar < nParams; iPar++) 
    {
      currentParamIndex++;
      // Int_t thisParamIndex = iSys * fNParams + iPar;
      // If the parameter is constrained, minuit won't have a value
      // for it.  We'll need to copy the value from the 
      // corresponding constrained parameter.
      if(IsParameterConstrained(iSys, iPar)){
	Int_t priorIndex = GetConstrainedParamIndex(iSys, iPar);
	// cout<<"ThisIndex:\t"<<thisParamIndex<<".\t"<<parameters[thisParamIndex]<<endl;
	// cout<<"PriorIndex:\t"<<priorIndex<<".\t"<<parameters[priorIndex]<<endl;
	assert(currentParamIndex > priorIndex);
	parameters[currentParamIndex] = parameters[priorIndex];
	constrainedParams++;
  	continue;
      }
      parameters[currentParamIndex] = par[currentParamIndex - constrainedParams];
    }
  }

  // Break up the total parameter vector into a mini-vector for
  // each PairSystem. Then pass the vector to the system and 
  // get chisquare
  totalChisquare = 0;
  Int_t nCopiedParams = 0;
  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    Int_t nParams = fNParamsPerSystem[iSys];
    vector<Double_t> pairSystemPars(&parameters[nCopiedParams],
				    &parameters[nCopiedParams + nParams]);
    nCopiedParams += nParams;
    fPairSystems[iSys]->SetLednickyParameters(pairSystemPars);
    totalChisquare += fPairSystems[iSys]->CalculateFitChisquare();
  }
  assert(nCopiedParams == fNParamsTotal);
  // cout<<"SetParametersAndFitEnd:\t"<<endl;
  // Timer();
  fChisquare = totalChisquare;
}

void Fitter::SetupInitialParameterVectors()
{
  // Get all the parameters put into arrays to be passed to Minuit.
  // Exclude constrained parameters from the list.

  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    for(Int_t iPar = 0; iPar < fNParamsPerSystem[iSys]; iPar++)
    {
      // Check to see if this parameter has been constrained 
      // to be the same as the parameter from an earlier system
      if(IsParameterConstrained(iSys, iPar)) {
	// This parameter is constrained. It has already been set 
	// up for a previous system.  Ignore it here.
	continue;
      }
      // New parameter.  Now set it up.
      TString parName = "";
      if(iPar > 3) {
	parName += "Norm";
	Int_t normIndex = iPar - 3;
	parName += normIndex;
	parName += fSystemNames[iSys];
      }
      else parName = fParamNames[iPar] + fSystemNames[iSys];
      fMinuitParNames.push_back(parName);
      fMinuitParInitial.push_back(fInitParams[iSys][iPar]);
      fMinuitParMinimum.push_back(fMinParams[iSys][iPar]);
      fMinuitParMaximum.push_back(fMaxParams[iSys][iPar]);
      fMinuitParIsFixed.push_back(fFixParams[iSys][iPar]);
    } // end parameter loop
  } // end system loop
}


void Fitter::SetupConstraint(Int_t param, vector<Int_t> systems)
{
  ParameterConstraint *constraint = new ParameterConstraint(param, systems, fPairSystems);
  fParamConstraints.push_back(constraint);
}


void Fitter::Timer()
{
  if(!fTimer) return;
  fTimer->Print();
  fTimer->Continue();
}
