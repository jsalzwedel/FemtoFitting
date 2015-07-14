

Double_t PairSystem::CalculateFitChisquare(vector<Double_t> pars)
{
  // Calculate the chisquare difference between the
  // correlation function data and the combined LednickyEqns
  
  // Update the fit parameters
  SetLednickyParameters(pars);


  // Make a TGraph object that sums all the LednickyEqns with
  // lambda parameters.  

  

  //Take the chisquare difference over the fit range.  
  Double_t chisquare;
  
  //...

  return chisquare;
}

void PairSystem::CreateNewLednickyEqn(TString name, Bool_t isIdentical, TH2D *transformMatrix/*,  other params? */)
{
  // Called by Fitter
  LednickyEqn *lednicky = new LednickyEqn(name, isIdentical, transformMatrix);
  fLednickyEqns.push_back(lednicky);
}

void PairSystem::SetLednickyParameters(vector<Double_t> pars)
{
  // Pass the new fit parameters to each LednickyEqn
  for(Int_t iLed = 0; iLed < fLednickyEqns.size(); iLed++){
    fLednickyEqns[iLed]->SetParameters(pars);
  }
}

PairSystem::PairSystem(TH1D *cfData, Bool_t isIdenticalPrimary)
{
  fCF = cfData;
  
  // Create a Lednicky eqn for primary-primary correlation function
  // and for each of the residual correlations.
  
  // Make the primary-primary lednicky eqn.  
  TH2D *transformLamLam; // = ...
  CreateNewLednickyEqn("LamLam", isIdenticalPrimary, NULL);
  
  // Make all the residual correlation lednicky eqns
  TH2D *transformLamSig; // = ... 
  CreateNewLednickyEqn("LamSig", kFALSE, transformLamSig);

  TH2D *transformLamXiC; // = ... 
  CreateNewLednickyEqn("LamXiC", kFALSE, transformLamXiC);

  TH2D *transformLamXi0; // = ... 
  CreateNewLednickyEqn("LamXi0", kFALSE, transformLamXi0);

  TH2D *transformSigSig; // = ... 
  CreateNewLednickyEqn("SigSig", isIdenticalPrimary, transformSigSig);  // This is identical unless we have particle-antiparticle system

  TH2D *transformSigXiC; // = ... 
  CreateNewLednickyEqn("SigXiC", kFALSE, transformSigXiC);
 
  TH2D *transformSigXi0; // = ... 
  CreateNewLednickyEqn("SigXi0", kFALSE, transformSigXi0);
}

PairSystem::~PairSystem()
{
  if(fCf){
    delete fCF;
    fCF = NULL;
  }
  for(Int_t i = 0; i < fLednickyEqns.size(); i++)
  {
    if(!fLednickyEqns[i]) continue;
    delete fLednickyEqns[i];
    fLednickyEqns[i] = NULL;
  }
}

void PairSystem::ReadInLambdaParams()
{

}
