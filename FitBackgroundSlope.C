TF1 *FitBackgroundSlope(TH1D *cf)
{
  // Fit the large k* region of a correlation function with a line
  // To extract the background slope.
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  // Setup TF1
  TF1 *line = new TF1("line","[0]*x + [1]", 0.0, 2.);
  line->SetParName(0, "Slope");
  line->SetParName(1, "Intercept");

  // set fit range

  // Fit and draw
  TCanvas *c1 = new TCanvas("Original", "Original");
  cf->Fit(line,"", "", 0.5, 1.);
  cf->DrawCopy();
  line->SetRange(0., 2.);
  line->DrawCopy("same");
  return line;
}

TH1D *CorrectCFSlope(TH1D *cf, TF1 *background)
{
  // Take a correlation function, divide out the background function,
  // and save the CF to a new file.

  TString newName = cf->GetName();
  newName += "BkgCorrected";
  // TH1D *correctedCF = (TH1D*) cf->Clone(newName);
  TH1D *correctedCF = new TH1D("tmp", "tmp", 10, 0, 10);
  cf->Copy(*correctedCF);
  correctedCF->SetName(newName);
  correctedCF->Divide(background);

  TCanvas *c2 = new TCanvas("Corrected","Corrected");
  correctedCF->DrawCopy("phiste");
  TF1 *flatLine = new TF1("flat","1", 0., 2.);
  flatLine->DrawCopy("same");
  return correctedCF;
}


void GenerateBackgroundCorrections(TString fileName, TString histName)
{

  // TString fileName = "/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root";
  // TString histName = "CombinedLLAA30-50KstarMomCorrected";

  TFile * inFile = new TFile(fileName);
  TH1D *cf = (TH1D*) inFile->Get(histName);

  TF1 *fit = FitBackgroundSlope(cf);

  TH1D *correctedCF = CorrectCFSlope(cf, fit);

  SaveCorrectedCF(correctedCF);
}

void SaveCorrectedCF(TH1D *cf)
{
  TString outputFileName = "BkgCorrectedCFs.root";
  TFile outFile(outputFileName,"update");
  outFile.cd();
  cf->Write(cf->GetName(), TObject::kOverwrite);
}

void RunAllCorrections()
{
  vector<TString> fileNames;
  fileNames.push_back("/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsCombinedLLAAMomCorrected.root");
  fileNames.push_back("/home/jai/Analysis/lambda/AliAnalysisLambda/Results/AnalysisResults/cfsLamALamKstarMomCorrected.root");

  vector<TString> histNames;
  histNames.push_back("CombinedLLAA0-10KstarMomCorrected");
  histNames.push_back("CombinedLLAA10-30KstarMomCorrected");
  histNames.push_back("CombinedLLAA30-50KstarMomCorrected");

  histNames.push_back("LamALam0-10centrality_varBin5BothFieldsKstarMomCorrected");
  histNames.push_back("LamALam10-30centrality_varBin5BothFieldsKstarMomCorrected");
  histNames.push_back("LamALam30-50centrality_varBin5BothFieldsKstarMomCorrected");

  for(Int_t i = 0; i < 3; i++)
  {
    GenerateBackgroundCorrections(fileNames[0], histNames[i]);
  }

  for(Int_t i = 3; i < 6; i++)
  {
    GenerateBackgroundCorrections(fileNames[1], histNames[i]);
  }
  

}

