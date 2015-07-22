#ifndef SystemInfo_H
#define SystemInfo_H


class SystemInfo {
 public: 
  SystemInfo(TString systemName, Double_t lambdaParam, TH2D *transformMatrix, Bool_t isIdenticalPair);
  ~SystemInfo();
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
