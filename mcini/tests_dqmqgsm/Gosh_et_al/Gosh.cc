{
#include "Riostream.h"
#include "string.h"
#include "TLatex.h"
#include <cmath>

gROOT->Reset();
gROOT->ProcessLine("#include <vector>");
gROOT->LoadMacro("ReadData.cc");
  bool graypalette = 0;
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(62);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);

// Make a canva and a pad
TCanvas* depo = new TCanvas("depo","Gosh et al",200,10,1000,800);
TPad* pad = new TPad("pad","The pad with the historgram",0.02,0.02,0.98,0.98,21);
pad->SetBorderMode(0); pad->SetFillColor(kWhite); pad->SetTickx(); pad->SetTicky(); pad->Draw();
pad->cd();
  
TH2D* ptHistoBack = new TH2D("","",1000,0,1000,1000,0,0.25);
  ptHistoBack->SetXTitle("p_{T}, MeV/c per nucleon");
  ptHistoBack->SetYTitle("P(p_{T})");
  ptHistoBack->GetXaxis()->SetTitleOffset(0.8);
  ptHistoBack->GetYaxis()->SetTitleOffset(1.);
  ptHistoBack->GetXaxis()->SetNdivisions(210);
  ptHistoBack->GetXaxis()->SetLabelOffset(0.007);
  ptHistoBack->SetStats(kFALSE);
ptHistoBack->Draw();

float xVal[100];
float yVal[100];
float xErr[100];
float yErr[100];


for(int i=0; i<100; i++) {
  xErr[i] = 20;
  yErr[i] = sqrt(500*yVal[i])/500.;
}


int nlines;
nlines = ReadData("gosh_pt.dat", xVal, 1., yVal, 0.01);

double integr = 0;
double pZ_init = 0;

for(int i=0; i<100; i++) {
  integr += yVal[i] ;
}
//std::cout<<integr<<"\n";


  TH1D* ptHisto = new TH1D("momentum_histo","", nlines-1,0.,1000);
  TH1D* ptHistoExp = new TH1D("momentum_histo","", nlines-1,0.,1000);

std::vector<double>* px;
std::vector<double>* py;
std::vector<double>* pz;
std::vector<double>* Z;
std::vector<double>* A;


TFile* ReadFile = new TFile("OOPt.root");
TTree *tree = (TTree*) ReadFile->Get("Glauber");
tree->SetBranchAddress("Z_on_A", &Z);
tree->SetBranchAddress("A_on_A", &A);
tree->SetBranchAddress("pX_on_A", &px);
tree->SetBranchAddress("pY_on_A", &py);
tree->SetBranchAddress("pZ_on_A", &pz);

TTree* init_par = (TTree*) ReadFile->Get("Conditions");
init_par->SetBranchAddress("pZ_in_MeV_on_A", &pZ_init);
init_par->GetEntry(0);

TGraphErrors* gosh1994 = new TGraphErrors(nlines-1, xVal, yVal, xErr, yErr);

  gosh1994->SetMarkerStyle(8);
  gosh1994->SetMarkerColor(kBlue);
  gosh1994->SetLineColor(kBlue);
  gosh1994->SetLineStyle(kSolid);
  gosh1994->SetMarkerSize(1.);
  gosh1994->SetLineWidth(1);

for(int k = 0; k < nlines; k++){ptHistoExp->Fill(xVal[k],yVal[k]);}

double norm = tree->GetEntries();

	for(int k = 0; k < tree->GetEntries(); k++){
		tree->GetEntry(k);
		
		double pX = 0;
		double pY = 0;
		double pZ = 0;
		double pT = 0;
		double tan_theta = 0;
		double theta = 0;
		double eta = 0;			           

		for(int k = 0; k< px->size(); k++ ){pX=px->at(k);pY=py->at(k);pZ=pz->at(k);
		                                   tan_theta = std::sqrt(pX*pX + pY*pY)/pZ;
						    theta = atan(tan_theta);
						    //pT = std::sqrt(pX*pX + pY*pY)/A->at(k);
						    pT = pZ_init*tan_theta/A->at(k);
						    ptHisto->Fill(pT);
						   
						   }

	}

delete ReadFile;

//gosh1994->Draw("P");
ptHisto->Scale(1./ptHisto->Integral());
ptHisto->SetLineWidth(2.);
ptHistoExp->SetLineWidth(2.);
ptHisto->SetLineStyle(kSolid);
ptHistoExp->SetLineStyle(kDashed);
ptHisto->SetLineColor(kBlue);
ptHistoExp->SetLineColor(kBlack);
//ptHisto->Scale(25.*0.666/norm);
ptHisto->Draw("HIST SAME");
ptHistoExp->Draw("HIST SAME");

auto leg = new TLegend(0.23,0.58,0.84,0.79,"2.1A GeV O on light targets ");
leg->AddEntry(ptHisto,"AAMCC, O-O ");
//leg->AddEntry(gosh1994, "Gosh et al (1994), O-CNO");
leg->AddEntry(ptHistoExp, "Gosh et al (1994), O-CNO");
leg->Draw();

pad->Update();

pad->Print("gosh_pt.png");

}
