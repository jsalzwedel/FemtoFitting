//**********************************************************
// Fitting correlation functions
//**********************************************************

#include "Fitter.h"
#include "PairSystem.h"

int main(int argc, char **argv)
{
  
  Fitter myFitter();
  myFitter.CreatePairSystems();
  myFitter.SetFitOptions();

  myFitter.CreateMinuit(/* */);
  myFitter.DoFitting();
  myFitter.SaveOutputPlots();
  
  
}

Fitter::Fitter()
{

}

Fitter::~Fitter()
{

}


void Fitter::CreatePairSystems()
{
  // Create all the pair systems objects that will be used
  // in the fitting


  //Example
  


}

void Fitter::CreateMinuit(/* */)
{
  // Create a TMinuit object.
  // Define, set and fix parameters.

  Int_t nParameters = 0;

  // Using the Fit Options and the number of LednickyEqns,
  // determine how many fit parameters there will be.

  fMinuit = new TMinuit(nParameters);
  fMinuit->SetFCN(SetParametersAndFit);

  // Define the fit parameters
  // for(int iPar = 0; iPar < nParameters; iPar++){
  //   fMinuit->DefineParameter(iPar, parNames[iPar], parInit[iPar],  0.001, parMin[iPar],  parMax[iPar]);
  //   if(parIsFixed[iPar]) myMinut->FixParameter(iPar);
  // }
  

}

void Fitter::DoFitting()
{
  // Run the fit procedure
  Double_t arglist[5] = {0,0,0,0,0}; //Arguments that can be passed with Minuit commands
  Int_t errFlag = 0;
  // arglist[0] = 1;
  // fMinuit->mnexcm("CALL FCN", arglist, 1, errFlag);

  // Set how verbose the output is (from no output at -1, to max at 3)
  arglist[0] =0;
  fMinuit->mnexcm("SET PRINT", arglist, 1, errFlag);

  // Maybe have all output results go to a file?


  // Run Migrad with a specified max number of calls
  arglist[0] = fMaxMinuitCalls;
  myMinuit->mnexcm("MIGRAD", arglist, 1, errFlag);

  // If outputting to file, turn return output to terminal now
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
  // -location of data
  // -location of transform matrices
  
}

void Fitter::SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag)
{
  // Take the TMinuit parameters par, set them in the Lednicky eqns,
  // and get the resulting Chisquare of the fit.

}
