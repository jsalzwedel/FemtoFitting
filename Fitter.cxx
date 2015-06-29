//**********************************************************
// Fitting correlation functions
//**********************************************************

int main(int argc, char **argv)
{
  //Determine how many 

}

void SetFitOptions()
{
  // Set things like:
  // -How many correlations to fit
  // -What parameters to fix
  // -What parameters to share
  // -location of data
  // -location of transform matrices
  
}

void CreateLednickyEqns()
{
  // Create all the lednicky eqn objects that will be used
  // in the fitting
  // Maybe rendered obsolete by AddResidual?  Will still need
  // to create the primary-primary correlation function.
}

TMinuit *CreateMinuit()
{
  // Create a TMinuit object.
  // Define, set and fix parameters.

}

void AddResidual(TString residualFileName, TString residualHistName)
{
  // Maybe makes a ResidualCorrelation object?  Object would have
  // its own Lednicky Eqn object?
  // Also, a pointer to a residual correlation histogram

}

void SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag)
{
  // Take the TMinuit parameters par, set them in the Lednicky eqns,
  // and get the resulting Chisquare of the fit.

}
