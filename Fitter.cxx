//**********************************************************
// Fitting correlation functions
//**********************************************************

#include "Fitter.h"
#include "LednickyInfo.h"
#include "PairSystem.h"
#include "ParameterConstraint.h"



Fitter::Fitter():
  fNParams(5),
  fNSystems(0),
  fMaxMinuitCalls(10000),
  fUseEstimatedLambdaParams(kTRUE)
{
  fParamNames.push_back("Radius");
  fParamNames.push_back("ReF0");
  fParamNames.push_back("ImF0");
  fParamNames.push_back("D0");
  fParamNames.push_back("Norm");
  
}

Fitter::~Fitter()
{
  for(Int_t i = 0; i < fPairSystem.size(); i++)
  {
    if(!fPairSystems[i]) continue;
    delete fPairSystems[i];
    fPairSystems[i] = NULL;
  }
  for(Int_t i = 0; i < fParamConstraints.size(); i++)
  {
    if(!fParamConstraints[i]) continue;
    delete fParamConstraints[i];
    fParamConstraints[i] = NULL;
  }
  if(fMinuit) {delete fMinuit; fMinuit = NULL;}
}

Bool_t Fitter::IsParameterConstrained(const Int_t currentSys, const Int_t currentPar)
{
  // Check to see if this parameter has been constrained 
  // to be the same as the parameter from an earlier system
  
  for(Int_t iCon = 0; iCon < fParamConstraints.size(); iCon++)
  {
    // Check for constraints on this type of parameter
    ParameterConstraint *constraint = fParamConstraints[iCon];
    if(constraint->GetConstrainedParam() != par) continue;
    
    vector<Int_t> consSystems = fParamConstraints[iCon]->GetConstrainedSystems();

    for(Int_t iSys = 1; iSys < consSystems; iSys++)
    {
      if(consSystems[iSys] == sys) return true;
    } 
  }
  return kFALSE;
}

// void Fitter::CreateAllPairSystems(Int_t configuration)
// {
//   // Create all the pair systems objects that will be used
//   // in the fitting.  Call CreateSinglePairSystem for each
//   // combination of pair type and centrality that will be 
//   // fit.

//   vector<TString> fileNames;
//   vector<TString> histNames;
//   vector<Bool_t>  isIdenticalPrimary;

//   GetHistConfiguration(configuration, fileNames, histNames, isIdenticalPrimary);
//   fNSystems = fileNames.size();  
//   assert(fNSystems == histNames.size());
//   assert(fNSystems > 0);

//   for(int iSystem = 0; iSystem < fNSystems; iSystem++)
//   {
//     TFile inFile((fileNames[iSystem]), "read");
//     TH1D *cf = inFile.Get(histNames[iSystem]);
//     assert(cf);
//     cf->SetDirectory(0);
//     PairSystem *system = new PairSystem(cf, isIdenticalPrimary[iSystem]);
//     fPairSystems.push_back(system);
//   }
// }


void Fitter::CreatePairSystem(TString simpleName, TString fileName, TString histName, const vector<LednickyInfo> &ledInfo, vector<Double_t> initParams, vector<Bool_t> fixParams)
{
  TFile inFile(fileName, "read");
  TH1D *cf = inFile.Get(histName);
  assert(cf);
  cf->SetDirectory(0);
  PairSystem *system = new PairSystem(cf, ledInfo, simpleName);
  fPairSystems.push_back(system);
  fSystemNames.push_back(simpleName);
  fInitParams.push_back(initParams);
  fFixParams.push_back(fixParams);
}

void Fitter::CreateMinuit(/* */)
{
  // Create a TMinuit object.
  // Define, set and fix parameters.


  // Using the Fit Options and the number of LednickyEqns,
  // determine how many fit parameters there will be.

  fMinuit = new TMinuit(fMinuitParNames.size());
  fMinuit->SetFCN(SetParametersAndFit);
  InitializeParameters(fMinuit);


}

void Fitter::DoFitting()
{
  // Run the fit procedure
  Double_t arglist[5] = {0,0,0,0,0}; //Arguments that can be passed with Minuit commands
  Int_t errFlag = 0;
  // arglist[0] = 1;
  // fMinuit->mnexcm("CALL FCN", arglist, 1, errFlag);

  // Set how verbose the output is (from no output at -1, to max at 3)
  arglist[0] = 0;
  fMinuit->mnexcm("SET PRINT", arglist, 1, errFlag);

  // Maybe have all output results go to a file?


  // Run Migrad with a specified max number of calls
  arglist[0] = fMaxMinuitCalls;
  myMinuit->mnexcm("MIGRAD", arglist, 1, errFlag);

  // If outputting to file, return output to terminal now
}

void Fitter::GetHistConfiguration(Int_t config, vector<TString> &fileNames, vector<TString> &histNames, vector<Bool_t> &isIdenticalPrimary)
{
  // Returns a vector of TFile names and corresponding hist names.
  // Can put many different correlation function histograms in 
  // here.  Call any combination of them by adding their numbers
  // to the config variable.

  if(kLL010 & config) {
    fileNames.push_back("/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root");
    histNames.push_back("CombinedLLAA0-10KstarMomCorrected");
    isIdenticalPrimary.push_back(kTRUE);
    // Add more as needed
  }
  if(kLL1030 & config) {
    fileNames.push_back("/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root");
    histNames.push_back("CombinedLLAA10-30KstarMomCorrected");
    isIdenticalPrimary.push_back(kTRUE);
    // Add more as needed
  }
  if(kLL3050 & config) {
    fileNames.push_back("/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsLamALamKstarMomCorrected.root");
    histNames.push_back("CombinedLLAA30-50KstarMomCorrected");
    isIdenticalPrimary.push_back(kTRUE);
    // Add more as needed
  }
  if(kLA010 & config) {
    fileNames.push_back("/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root");
    histNames.push_back("LamALam0-10centrality_varBin5BothFieldsKstarMomCorrected");
    isIdenticalPrimary.push_back(kFALSE);
    // Add more as needed
  }
  if(kLA1030 & config) {
    fileNames.push_back("/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root");
    histNames.push_back("LamALam10-30centrality_varBin5BothFieldsKstarMomCorrected");
    isIdenticalPrimary.push_back(kFALSE);
    // Add more as needed
  }
  if(kLA3050 & config) {
    fileNames.push_back("/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root");
    histNames.push_back("LamALam30-50centrality_varBin5BothFieldsKstarMomCorrected");
    isIdenticalPrimary.push_back(kFALSE);
    // Add more as needed
  }

}

void Fitter::InitializeParameters(TMinuit *minuit)
{
  // Define the fit parameters

  Double_t startingStepSize = 0.001;
  
  // Allow the parameters to have no bounds
  Double_t parameterMinimum = 0;
  Double_t parameterMaximum = 0;

  for(int iPar = 0; iPar < fMinuitParNames.size(); iPar++){
    minuit->DefineParameter(iPar, fParNames[iPar], fParInitial[iPar],  startingStepSize, parameterMinimum, parameterMaximum);
    if(fParIsFixed[iPar]) myMinut->FixParameter(iPar);
  }
  
}


void Fitter::SaveOutputPlots()
{
  // Save output plots for each fit
}

void Fitter::SetFitOptions(const Int_t constraintConfig)
{
  // Set things like:
  // -How many correlations to fit
  // -What parameters to fix
  // -What parameters to share


  SetupParameterVectors();
  SetupParameterConstraints(constraintConfig);
  SetupInitialParameters();
 
  
}

void Fitter::SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag)
{
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
      if(IsParamConstrained(iSys, iPar, true)){
	Int_t priorIndex = GetConstrainedParamIndex(iSys, iPar);
	assert(thisParamIndex > priorIndex);
	parameters[thisParamIndex] = parameters[priorIndex];
	constrainedParams++;
  	continue;
      }
      parameters[thisParamIndex] = par[thisParamIndex - constrainedParams];
    }
  }

  assert(sizeof(par)/sizeof(Double_t) == fNParams * fNSystems - constrainedParams);

  // Break up the total parameter vector into a mini-vector for
  // each PairSystem. Then pass the vector to the system and 
  // get chisquare
  totalChisquare = 0;
  for(Int_t iSys = 0; iSys < fNSystems; iSys++)
  {
    
    vector<Double_t> pairSystemPars(&parameters[iSys * fNParams],
				    &parameters[(iSys +1) * fNParams]);
    fPairSystems[iSys]->SetLednickyParams(pairSystemPars);
    totalChisquare += fPairSystems[iSys]->CalculateFitChisquare();
  }
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
      if(IsParameterConstrained(thisSys, iPar)) {
	// This parameter is constrained. It has already been set 
	// up for a previous system.  Ignore it here.
	continue;
      }
      // New parameter.  Now set it up.
      TString parName = fParamNames[iPar] + fSystemNames[iSys];
      fMinuitParNames.push_back(parName);
      fMinuitParInitial.push_back(fInitParams[iSys][iPar]);
      fMinuitParIsFixed.push_back(fInitParams[iSys][iPar]);
    } // end parameter loop
  } // end system loop
}







void Fitter::SetupParameterConstraints(const Int_t config)
{
  
  if(config /*  & ... */){
    ConstrainF0D0();
  }
  if(config /*  & ... */){
    ConstrainRadii();
  }

}

void Fitter::SetupParameterVectors()
{
  // Initialize parameter values for Minuit
  // Set min and max values

 
  
  // Set the size of the parameter vectors.
  fParNames.resize(fNParams*fNSystems);
  fParInitial.resize(fNParams*fNSystems);
  fParMinimum.resize(fNParams*fNSystems);
  fParMaximum.resize(fNParams*fNSystems);
  fParCurrent.resize(fNParams*fNSystems);
  fParIsFixed.resize(fNParams*fNSystems);

  // Basic param setup

}

// void Fitter::SetUseEstimatedLambdaParams(Bool_t useParam)
// {
//   fUseEstimatedLambdaParams = useParam;
//   if(!useParam) fNParams = 6; // Or 5 + however many lambda params to fit
// }

void Fitter::ConstrainRadii()
{
  // Constrain radii between particle-particle and particle-antiparticle systems
  ParamType param = kRadii;  
  
  Int_t systemsArr010[2] = {kLL010, kLa010};
  vector<Int_t> systems010(systemsArr010);
  ParameterConstraint *constraint = new ParameterConstraint(param, systems010);
  fParamConstraints.push_back(constraint);

  Int_t systemsArr1030[2] = {kLL1030, kLa1030};
  vector<Int_t> systems1030(systemsArr1030);
  constraint = new ParameterConstraint(param, systems1030);
  fParamConstraints.push_back(constraint);

  Int_t systemsArr3050[2] = {kLL3050, kLa3050};
  vector<Int_t> systems3050(systemsArr3050);
  constraint = new ParameterConstraint(param, systems1030);
  fParamConstraints.push_back(constraint);
}

void Fitter::ConstrainF0D0()
{
  // Constrain all LL systems to use the same scattering
  // length and d0.  Constrain the LA systems the same way.

  Int_t systemsArrLL[3] = {kLL010, kLL1030, kLL3050};
  vector<Int_t> systemsLL(systemsArrLL);

  Int_t systemsArrLA[3] = {kLA010, kLA1030, kLA3050};
  vector<Int_t> systemsLA(systemsArrLA);

  ParamType param = kF0Real;  
  ParameterConstraint *constraint = new ParameterConstraint(param, systemsLL);
  fParamConstraints.push_back(constraint);

  constraint = new ParameterConstraint(param, systemsLA);
  fParamConstraints.push_back(constraint);

  param = kF0Imag;
  constraint = new ParameterConstraint(param, systemsLL);
  fParamConstraints.push_back(constraint);

  constraint = new ParameterConstraint(param, systemsLA);
  fParamConstraints.push_back(constraint);

  param = kD0;
  constraint = new ParameterConstraint(param, systemsLL);
  fParamConstraints.push_back(constraint);

  constraint = new ParameterConstraint(param, systemsLA);
  fParamConstraints.push_back(constraint);
}
