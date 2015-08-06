//**********************************************************
// Fitting correlation functions
//**********************************************************

#include "Fitter.h"
#include <iostream>
#include "TH1D.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TMath.h"

#include <assert.h>
#include "TCanvas.h"

using namespace std;


Fitter::Fitter():
  fNParams(5),
  fFixedParams(0),
  fNSystems(0),
  fMaxMinuitCalls(10000),
  fUseEstimatedLambdaParams(kTRUE),
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



void Fitter::CreatePairSystem(TString simpleName, TString fileName, TString histName, Int_t sysType, const vector<LednickyInfo> &ledInfo, vector<Double_t> initParams, vector<Double_t> minParams, vector<Double_t> maxParams, vector<Bool_t> fixParams)
{
  TFile inFile(fileName, "read");
  TH1D *cf = (TH1D*) inFile.Get(histName);
  assert(cf);
  cf->SetDirectory(0);
  PairSystem *system = new PairSystem(cf, ledInfo, simpleName, sysType);
  fPairSystems.push_back(system);
  fSystemNames.push_back(simpleName);
  fInitParams.push_back(initParams);
  fMinParams.push_back(minParams);
  fMaxParams.push_back(maxParams);
  fFixParams.push_back(fixParams);
  fNSystems++;
  for(UInt_t i = 0; i < fixParams.size(); i++)
  {
    if(fixParams[i]) fFixedParams++;
  }
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

  // Maybe have all output results go to a file?

  Timer();
  // Run Migrad with a specified max number of calls
  arglist[0] = fMaxMinuitCalls;
  arglist[1] = 0.1;
  fMinuit->mnexcm("MIGRAD", arglist, 1, errFlag);
  Timer();
  
  fMinuit->mnexcm("SHOw CORrelations", arglist, 1, errFlag);

  cout<<"Finalchi2:\t"<<fChisquare<<endl
      <<"Fit bins:\t"<<fFitBins<<endl
      <<"Actual Minuit Pars:\t"<<fMinuitParNames.size()<<endl
      <<"Fixed params\t"<<fFixedParams<<endl
      <<"Free Params:\t"<<fMinuitParNames.size() - fFixedParams<<endl
      <<"Chi2/ndf:\t"<<GetChisquarePerNDF()<<endl
      <<"Pvalue:\t"<<GetPvalue()<<endl;

  // If outputting to file, return output to terminal now
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
	Int_t sysIndexPrior = 100;
	for(Int_t iPairSys = 0; iPairSys < fNSystems; iPairSys++)
	{
	  if(sysTypePrior == fPairSystems[iPairSys]->GetSystemType())
	  {
	    sysIndexPrior = iPairSys;
	  }
	}
	assert(sysIndexPrior != 100);
	// Return the absolute index of the first parameter in the 
	// constraint
        Int_t absoluteIndex = sysIndexPrior * fNParams + parIndex;
	return absoluteIndex;
      }
    } 
  }
  // Should not get to this point
  return 100000;
}

Double_t Fitter::GetChisquarePerNDF()
{
  Double_t ndf = 1.*fFitBins - (1.*fMinuitParNames.size() - 1.*fFixedParams);
  return fChisquare/ndf;
}

Double_t Fitter::GetPvalue()
{
  Double_t ndf = 1.*fFitBins - (1.*fMinuitParNames.size() - 1.*fFixedParams);
  return TMath::Prob(fChisquare,ndf); 
}

void Fitter::InitializeMinuitParameters(TMinuit *minuit)
{

  SetupInitialParameters();
  // Define the fit parameters
  fMinuit = minuit;
  Double_t startingStepSize = 0.1;
  
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

void Fitter::SaveOutputPlots()
{
  // Save output plots for each fit
  TFile outFile("FitResults.root","update");
  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    TH1D *cf = fPairSystems[iSys]->GetCF();
    cf->SetAxisRange(0.7,1.05,"Y");
    cf->SetAxisRange(0.,1.,"X");
    cf->Write(cf->GetName(), TObject::kOverwrite);
    TGraph *g = fPairSystems[iSys]->GetCombinedTGraph();
    g->Write(g->GetName(), TObject::kOverwrite);
    TCanvas c1("c1","c1");
    cf->DrawCopy();
    g->Draw("same");
    TString plotName = "Plot";
    plotName += g->GetName();
    c1.SaveAs(plotName + ".pdf");
    c1.SaveAs(plotName + ".png");
    cout<<"Saved file "<<plotName<<endl;
  }
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
  
  // Take the TMinuit parameters, set the parameters for each
  // pair system, and get the resulting chisquare of the fit.

  // We'll copy par into the parameters vector, and insert 
  // any constrained parameters into their appropriate
  // positions
  vector<Double_t> parameters(fNSystems * fNParams);
  Int_t constrainedParams = 0;
  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    for(Int_t iPar = 0; iPar < fNParams; iPar++) 
    {
      Int_t thisParamIndex = iSys * fNParams + iPar;
      // If the parameter is constrained, minuit won't have a value
      // for it.  We'll need to copy the value from the 
      // corresponding constrained parameter.
      if(IsParameterConstrained(iSys, iPar)){
	Int_t priorIndex = GetConstrainedParamIndex(iSys, iPar);
	// cout<<"ThisIndex:\t"<<thisParamIndex<<".\t"<<parameters[thisParamIndex]<<endl;
	// cout<<"PriorIndex:\t"<<priorIndex<<".\t"<<parameters[priorIndex]<<endl;
	assert(thisParamIndex > priorIndex);
	parameters[thisParamIndex] = parameters[priorIndex];
	constrainedParams++;
  	continue;
      }
      parameters[thisParamIndex] = par[thisParamIndex - constrainedParams];
    }
  }

  // Break up the total parameter vector into a mini-vector for
  // each PairSystem. Then pass the vector to the system and 
  // get chisquare
  totalChisquare = 0;
  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    
    vector<Double_t> pairSystemPars(&parameters[iSys * fNParams],
				    &parameters[(iSys +1) * fNParams]);
    fPairSystems[iSys]->SetLednickyParameters(pairSystemPars);
    totalChisquare += fPairSystems[iSys]->CalculateFitChisquare();
  }
  // cout<<"SetParametersAndFitEnd:\t"<<endl;
  // Timer();
  fChisquare = totalChisquare;
}


void Fitter::SetupInitialParameters()
{
  // Get all the parameters put into arrays to be passed to Minuit.
  // Exclude constrained parameters from the list.

  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    for(Int_t iPar = 0; iPar < fNParams; iPar++)
    {
      // Check to see if this parameter has been constrained 
      // to be the same as the parameter from an earlier system
      if(IsParameterConstrained(iSys, iPar)) {
	// This parameter is constrained. It has already been set 
	// up for a previous system.  Ignore it here.
	continue;
      }
      // New parameter.  Now set it up.
      TString parName = fParamNames[iPar] + fSystemNames[iSys];
      fMinuitParNames.push_back(parName);
      fMinuitParInitial.push_back(fInitParams[iSys][iPar]);
      fMinuitParMinimum.push_back(fMinParams[iSys][iPar]);
      fMinuitParMaximum.push_back(fMaxParams[iSys][iPar]);
      fMinuitParIsFixed.push_back(fFixParams[iSys][iPar]);
    } // end parameter loop
  } // end system loop
}

void Fitter::Timer()
{
  if(!fTimer) return;
  fTimer->Print();
  fTimer->Continue();
}

void Fitter::SetupConstraint(Int_t param, vector<Int_t> systems)
{
  ParameterConstraint *constraint = new ParameterConstraint(param, systems, fPairSystems);
  fParamConstraints.push_back(constraint);
}
