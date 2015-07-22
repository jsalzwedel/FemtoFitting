void MoveTransformMatric(TString inFileName, TString matrixName)
{
  TFile inFile(inFileName);
  TH2D *transform = (TH2D*) inFile.Get(matrixName);
  
  TFile outFile("PreparedTransformMatrices.root","update");
  if(transform) {
    outFile.cd();
    transform->RebixX(4);
    transform->RebinY(4);
    transform->Write(transform->GetName(), TObject::kOverwrite);
  }
  else cout<<"Could not find matrix with that name"<<endl;
  outFile.Close();
  inFile.Close();
}
