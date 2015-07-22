#ifndef LednickyInfo_H
#define LednickyInfo_H


class LednickyInfo {
 public: 
  LednickyInfo(TString systemName, Double_t lambdaParam, TH2D *transformMatrix, Bool_t isIdenticalPair);
  ~LednickyInfo();
  TString GetSystemName() const {return fSystemName;};
  Double_t GetLambdaParam() const {return fLambdaParam;};
  TH2D *GetTransformMatrix() const {return fTransformMatrix;};
  Bool_t GetIsIdenticalPair() const {return fIsIdenticalPair;};
  

 private:
  TString fSystemName;
  Double_t fLambdaParam;
  TH2D* fTransformMatrix;
  Bool_t fIsIdenticalPair;
};



#endif
