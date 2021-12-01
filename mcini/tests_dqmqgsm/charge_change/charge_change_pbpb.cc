{
#include "Riostream.h"
#include "string.h"
#include "TLatex.h"

gROOT->Reset();
gROOT->ProcessLine("#include <vector>");
gROOT->LoadMacro("ReadData.cc");
// Styles "Plain", "Bold", "Video", "Pub" are available
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


TCanvas* depo = new TCanvas("depo","NPP",200,10,1000,800);
TPad* pad = new TPad("pad","The pad with the historgram",0.02,0.02,0.98,0.98,21);
pad->SetBorderMode(0); pad->SetFillColor(kWhite); pad->SetTickx(); pad->SetTicky(); pad->Draw();


pad->cd(); pad->SetLogy(1);

  
  TH2D* np_histo = new TH2D("np_histo ","", 2000 ,-0.5 , 82.5, 36 ,2,1e4 );
  np_histo->SetStats(kFALSE);
  np_histo->SetXTitle("Z");
  np_histo->SetYTitle("#frac{d #sigma}{dZ}, mb");
  np_histo->SetLabelSize(0.038);
  np_histo->GetXaxis()->SetTitleOffset(0.8);
  np_histo->GetYaxis()->SetTitleOffset(1.);
  np_histo->GetXaxis()->SetNdivisions(210);
  np_histo->GetXaxis()->SetLabelOffset(0.007);
  
TH1D* charge_change = new TH1D("","",83,-0.5,82.5);
TH1D* charge_change_evap = new TH1D("","",83,-0.5,82.5);
TH1D* charge_change_fission = new TH1D("","",83,-0.5,82.5);
TH1D* charge_change_mf = new TH1D("","",83,-0.5,82.5);
charge_change->SetLineColor(kBlack);
charge_change_evap->SetLineColor(4);
//charge_change_evap->SetLineStyle(kDashed);
charge_change_fission->SetLineColor(8);
//charge_change_fission->SetLineStyle(kDotted);
charge_change_mf->SetLineColor(40);
//charge_change_mf->SetLineStyle(kDashDotted);

double sigma = 5.43;
double norm = 100000;
float b;
double kinEn;
std::vector<double>* Z;
std::vector<double>* A;

int nlines;
float xVal[100];
float yVal[100];
float xErr[100];
float yErr[100];


for(int i=0; i<100; i++) {
  xErr[i] = 0.;
  yErr[i] = 0.;
}

nlines = ReadData("../PbPb_158A_GeV.dat", xVal, 1., yVal, 1.);
TGraphErrors* PbPb_GSI_graph = new TGraphErrors(nlines-1, xVal, yVal, xErr, yErr);

  PbPb_GSI_graph->SetMarkerStyle(20);
  PbPb_GSI_graph->SetMarkerColor(kRed+1);
  PbPb_GSI_graph->SetLineColor(kBlue);
  PbPb_GSI_graph->SetLineStyle(kSolid);
  PbPb_GSI_graph->SetMarkerSize(1.);
  PbPb_GSI_graph->SetLineWidth(0);

TFile* ReadedFile1_em = new TFile("../bolpb.root");
TH2D*  ReadedHisto1_em = (TH2D*) ReadedFile1_em->Get("h999");
TH1D*  emHist = ReadedHisto1_em->ProjectionY();
emHist->SetLineColor(46);

TFile* ReadFile = new TFile("../pbpb158AGeV_100.root");
TTree *tree = (TTree*) ReadFile->Get("Glauber");
tree->SetBranchAddress("impact_parameter", &b);
tree->SetBranchAddress("A_on_A", &A);
tree->SetBranchAddress("Z_on_A", &Z);



TTree *init_par = (TTree*) ReadFile->Get("Conditions");

init_par->SetBranchAddress("Xsect_total", &sigma);
init_par->GetEntry(0);
sigma*=1000; 
norm = tree->GetEntries();
	for(int k = 0; k < tree->GetEntries(); k++){
		tree->GetEntry(k);

		double n = 0;
		double p = 0;
		double Fheavy_z = 0;
		double Sheavy_z = 0;
		double mult = 0;
		double totZ = 0;

		for(int z = 0; z< Z->size(); z++ ){
			totZ += Z->at(z);
			if(Z->at(z) > Fheavy_z){Sheavy_z = Fheavy_z; Fheavy_z = Z->at(z);}
			if(A->at(z) > 4){mult +=1;};
			}
		for(int z = 0; z< Z->size(); z++ ){			
			charge_change->Fill(Z->at(z));}

		if((Fheavy_z - Sheavy_z) > Fheavy_z-10 && (mult  == 1)){
		for(int z = 0; z< Z->size(); z++ ){			
			charge_change_evap->Fill(Z->at(z));}}

		if((Fheavy_z - Sheavy_z) < Fheavy_z*0.5 && (mult  == 2)){
		for(int z = 0; z< Z->size(); z++ ){			
			charge_change_fission->Fill(Z->at(z));}}

		if((Fheavy_z - Sheavy_z) < Fheavy_z*0.3333 && (mult  > 2)){
		for(int z = 0; z< Z->size(); z++ ){			
			charge_change_mf->Fill(Z->at(z));}}
			
	}

np_histo->Draw();

PbPb_GSI_graph->Draw("P");
charge_change->Scale(sigma/norm);
charge_change_evap->Scale(sigma/norm);
charge_change_fission->Scale(sigma/norm);
charge_change_mf->Scale(sigma/norm);
emHist->Draw("HIST SAME");
charge_change->Add(emHist);
charge_change->Draw("HIST SAME");
charge_change_evap->Draw("HIST SAME");
charge_change_fission->Draw("HIST SAME");
charge_change_mf->Draw("HIST SAME");

TLatex* sqrtS = new TLatex(4.5,12050, "PbPb, 158A GeV");
sqrtS->SetTextFont(42);
sqrtS->Draw();

TLegend* leg2 = new TLegend(0.265,0.605,0.83,0.86,"AAMCCv2, Hybrid");
leg2->AddEntry(charge_change,"all", "lp");
leg2->AddEntry(charge_change_evap,"Evaporation, (M_{A>4} = 1; Z_{max} - Z_{premax} > Z_{max} - 5)", "lp");
leg2->AddEntry(charge_change_fission,"Fission, (M_{A>4} = 2; Z_{max} - Z_{premax} < 0.5*Z_{max} )", "lp");
leg2->AddEntry(charge_change_mf,"MF, (M_{A>4} > 3; Z_{max} - Z_{premax} < Z_{max}/3 )", "lp");
leg2->AddEntry(emHist, "EM");
leg2->AddEntry(PbPb_GSI_graph,"exp; Dekhisi (2000)","lp");
leg2->Draw();

pad->Update();


gPad->Print("biasedPbPb11MeV.png");
}
