

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

void PairSystem::CreateNewLednickyEqn(Bool_t isIdentical, TH2D *transformMatrix/*,  other params? */)
{
  // Called by Fitter
  LednickyEqn *lednicky = new LednickyEqn(isIdentical, transformMatrix);
  fLednickyEqns.push_back(lednicky);
}

void PairSystem::SetLednickyParameters(vector<Double_t> pars)
{
  // Pass the new fit parameters to each LednickyEqn
  for(Int_t iLed = 0; iLed < fLednickyEqns.size(); iLed++){
    fLednickyEqns[iLed]->SetParameters(pars);
  }
}

PairSystem::PairSystem(TH1D *cfData)
{
  fCF = cfData;
  
  // Create a Lednicky eqn for primary-primary correlation function
  // and for each of the residual correlations.


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
