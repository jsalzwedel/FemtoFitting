//**********************************************************
// Fitting correlation functions
//**********************************************************

#include "Fitter.h"
#include "PairSystem.h"
#include "ParameterConstraint.h"

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
  myFitter.CreateMinuit(/* */);
  myFitter.DoFitting();
  myFitter.SaveOutputPlots();
  
  
}

Fitter::Fitter()
{
  fNParams(6); 
  fNSystems(0);
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
}


void Fitter::CreateAllPairSystems(Int_t configuration)
{
  // Create all the pair systems objects that will be used
  // in the fitting.  Call CreateSinglePairSystem for each
  // combination of pair type and centrality that will be 
  // fit.

  vector<TString> fileNames;
  vector<TString> histNames;
  vector<Bool_t>  isIdenticalPrimary;

  GetHistConfiguration(configuration, fileNames, histNames, isIdenticalPrimary);
  fNSystems = fileNames.size();  
  assert(fNSystems == histNames.size());
  assert(fNSystems > 0);

  for(int iSystem = 0; iSystem < fNSystems; iSystem++)
  {
    TFile inFile((fileNames[iSystem]), "read");
    TH1D *cf = inFile.Get(histNames[iSystem]);
    assert(cf);
    cf->SetDirectory(0);
    PairSystem *system = new PairSystem(cf, isIdenticalPrimary[iSystem]);
    fPairSystems.push_back(system);
  }
}

void Fitter::CreateMinuit(/* */)
{
  // Create a TMinuit object.
  // Define, set and fix parameters.


  // Using the Fit Options and the number of LednickyEqns,
  // determine how many fit parameters there will be.

  fMinuit = new TMinuit(fNParams);
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

  for(int iPar = 0; iPar < fNParams; iPar++){
    minuit->DefineParameter(iPar, fParNames[iPar], fParInitial[iPar],  startingStepSize, fParMinimum[iPar], fParMaximum[iPar]);
    if(fParIsFixed[iPar]) myMinut->FixParameter(iPar);
  }
  
}


void Fitter::SaveOutputPlots()
{
  // Save output plots for each fit
}

void Fitter::SetFitOptions()
{
  // Set things like:
  // -How many correlations to fit
  // -What parameters to fix
  // -What parameters to share


  SetupParameters();
  SetupParameterConstraints();

 
  
}

void Fitter::SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag)
{
  // Take the TMinuit parameters, set them for each pair 
  // system, and get the resulting Chisquare of the fit.


  totalChisquare = 0;

  for(Int_t iSys = 0; iSys < fPairSystems.size(); iSys++)
  {
    vector<Double_t> newPars;
    for(Int_t iPar = 0; iPar < fNParams; iPar++) {
      newPars[iPar] = par[fNParams * iSys + iPar];
    }
    fPairSystems[iSys]->SetLednickyParams(newPars);
    totalChisquare += fPairSystems[iSys]->CalculateFitChisquare();
  }
}

void Fitter::SetupParameters()
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
