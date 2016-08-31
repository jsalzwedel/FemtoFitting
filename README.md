# FemtoFitting
Generalized Lednicky fit code for 1D femtoscopic correlation functions with residual correlations (no Coulomb interactions).

## Setup

Before running, set up user-specific options in Main.cxx via the following functions: PrepareLednickyInfo, UserSetupSystems, UserSetConstraints, and UserSetFitOptions.

* PrepareLednickyInfo
  1. Create a LednickyInfo object (contains info that is used to set up the form of the Lednicky Eqn)
  2. In the LednickyInfo constructor, set the lambda parameter (must be set - cannot be a free parameter in this fitter)
  3. In constructor, set the smear matrix (momentum resolution and/or residual correlation matrix). Can also pass NULL if no transform is needed (i.e. perfect resolution primary particles)
  4. In constructor, set whether particles in the pair are identical.
  5. Can also set mass values if you want to include sqrt(s) scaling of the scattering parameters for residual correlations.
  6. Add LednickyInfo object to ledInfo vector.
  7. Repeat steps 1-6 for each residual correlation.
* UserSetupSystems
  1. Add a new PairSystem to the fitter via fitter->AddPairAnalysisChisquareFit or AddPairAnalysisLogFit
    1. Make TString for the path to the ROOT file containing your correlation functions (or num/den histograms)
    2. Make TString for the path to the histogram *in* the root file.
    3. In initParamsArr, setup starting parameter values for radius, Ref0, Imf0, d0, linear bkg, quadratic bkg, and normalization
    4. Optional: Setup min/max parameter constraints (not recommended - screws with error bars)
    5. In fixParamsArr, setup which parameters should be fixed during fit process.
    6. Create the LednickyInfo vector via PrepareLednickyInfo.
    7. Choose a unique SystemType enumeration (see top of Main.cxx). The user can modify the list of enumerated values as desired.
    7. Add pair system to fitter, passing the above parameters
    8. Tell fitter if you want to use linear and/or quadratic background parameter (via fitter->SetUseLinearBkgPoly and SetUseQuadBkgPoly)
  2. Add more pair systems as needed (e.g. for different pair types, centrality ranges, or kT regions)
  2. If using existing parameter setup, can fix various params by uncommenting lines towards the top of UserSetupSystems
  3. Don't mix chisquare fitting with log-likelihood fitting!
* UserSetConstraints
  1. Setup any desired parameter constraints between different PairSystems. For example, might want to force Lambda-Lambda to use the same Ref0 value for each centrality bin.
  2. The parameters are identified via enums (see ParamType.h).
  3. The constrained PairSystems are identified via the SystemType enum that they were created with.
* UserSetFitOptions
  1. Set fit range
  2. Set upper and lower background fit range (if using non-flat bkg)
  3. Set output options as desired




## How to run

Compile code via make (Requires ROOT v6 for C++11 compatibility).  Run code with ./runMe

Fit results are output in the terminal as well as to ./FitResults/SystematicFitData.csv. Each time the fitter is run, a new line of results is appended

The column names can be found in FitResults/SystematicColumnNames.csv
