//...


class Fitter{
 public:
  Fitter();
  virtual ~Fitter();

 private:
  void SetParametersAndFit(Int_t& i, Double_t *x, Double_t &totalChisquare, Double_t *par, Int_t iflag);

  TMinuit *fMinuit;
  vector<*LednickyEqn> fLedEquations;
}
