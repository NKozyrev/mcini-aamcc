#include "TGaxis.h"
#include "TRandom.h"
#include "Riostream.h"
#include "string.h"
#include "TLatex.h"
{bool graypalette = 0;
  
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

	gStyle->SetOptTitle(0);
	//gPad->SetLogy();
	double sigma, sigma1;
	//double N2=2e+3;

   //TGaxis::SetMaxDigits(1); 
   TCanvas *c1 = new TCanvas("c1","hists with different scales",800,800);
   //c1 -> Divide(2,1,0.001,0.001);
   //c1 -> cd(1);
   gStyle->SetOptStat(kFALSE);
   std::vector<double>* A;
   std::vector<double>* Z;
   std::vector<double>* A1;
   std::vector<double>* Z1;

   TFile* ReadFile = new TFile("../pbpb158AGeV_20k_Emax8.root");
   TH1D* h1 = new TH1D("h1","", 31, -0.5, 30.5);
   TTree *tree = (TTree*) ReadFile->Get("Glauber");   
   tree->SetBranchAddress("A_on_B", &A);
   tree->SetBranchAddress("Z_on_B", &Z);
   TTree *init_par = (TTree*) ReadFile->Get("Conditions");
   //init_par->SetBranchAddress("Xsect_total", &sigma);
   TFile* ReadFile1 = new TFile("../pbpb158AGeV_20k_Emax8.root");//Pbpnrw_17GeV_Aladin_0_5perc_10k_oldFermitest
   TH1D* h2 = new TH1D("h2","", 31, -0.5, 30.5);
   TTree *tree1 = (TTree*) ReadFile1->Get("Glauber");   
   tree1->SetBranchAddress("A_on_B", &A1);
  tree1->SetBranchAddress("Z_on_B", &Z1);
   TTree *init_par1 = (TTree*) ReadFile1->Get("Conditions");
   // init_par1->SetBranchAddress("Xsect_total", &sigma1);
   init_par->GetEntry(0);
   init_par1->GetEntry(0);
   double p=0,p1=0, nn=0, dd=0;
	double N=1*1e+4;
   for(int k = 0; k < tree->GetEntries(); k++){
      tree->GetEntry(k);

      for (int fragment = 0; fragment < (A->size()); fragment++){
            /*if (A->at(fragment) != 1 || Z->at(fragment)!= 0)*/h1->Fill(A->at(fragment));
            if (A->at(fragment) == 1 && Z->at(fragment) == 1) p+=1;
            if (A->at(fragment) == 1 && Z->at(fragment) == 0) nn+=1;
            if (A->at(fragment) == 2 && Z->at(fragment) == 1) dd+=1;
         }
      }

   for(int k = 0; k < tree1->GetEntries(); k++){
      tree1->GetEntry(k);

      for (int fragment = 0; fragment < (A1->size()); fragment++){
            /*if (A1->at(fragment) != 1 || Z1->at(fragment)!= 0)*/h2->Fill(A1->at(fragment));
         if (A1->at(fragment) == 1 && Z1->at(fragment) == 1) p1+=1;
         }
      }
      p=p/N;  p1=p1/N;
cout << "neutrons: " << nn/N << "; protons: " << p << "; deutrons: " << dd/N << endl;
TH1D* h3 = new TH1D("h3", "h3 title", 9, -0.5,8.5);
   //h3->SetBinContent(0,0);
   TH2D* back = new TH2D("","",100,1,12,100,0,20);
   back->Draw("HIST SAME");
   back->SetYTitle("<N_{all}>");
   back->SetXTitle("A");
   
   h2->SetYTitle("<N_{all}>");
   h2->SetXTitle("A");
   h2->SetTitleOffset(0.9, "y");
   h2->SetTitleSize(0.04);
   h2->GetXaxis()->SetTitleOffset(0.5);
   h2->GetXaxis()->SetLabelSize(0.03);
   h2->GetYaxis()->SetLabelSize(0.03);
   h2->GetXaxis()->SetLabelFont(62);
   h2->GetYaxis()->SetLabelFont(62);


    h1->Scale(Float_t(1/double(tree->GetEntries())));
    h2->Scale(Float_t(1/double(tree1->GetEntries())));

	h1->SetLineStyle(kSolid);
	h1->SetLineWidth(3);
	h2->SetLineWidth(3);
   h1->SetLineColor(1);
   h2->SetMaximum(19.);


   h1->GetXaxis()->SetRange(1,12);
   h2->GetXaxis()->SetRange(1,12);
   h2->GetYaxis()->SetRange(0,30);
   h1->SetTitleSize(0.04, "xy");
   h1->SetTitleFont(62, "xy");
   h1->SetTitleOffset(0.9, "y");
   h2->SetLineStyle(7);
   h2->SetLineColor(kRed);
  
   h2->Draw("SAME");
   c1->Update();
   h3->SetBinContent(1,16);h3->SetBinContent(2,0.5);
   h3->GetXaxis()->SetRange(1,12);
   h3->SetLineStyle(9);
   h3->SetLineWidth(3);
   h3->SetLineColor(kBlue);
   h3->Draw("HIST SAME");
   //h1->Draw("HIST SAME");

//std::cout<<h1->GetBinContent(1)<<"\n";

   auto legend = new TLegend(0.35,0.61,0.71,0.74);
   //legend -> AddEntry(h1, "AAMCC, with noneqiulibrium fragmentation");
   legend -> AddEntry(h2, "AAMCC, w/o noneqiulibrium fragmentation");
   legend -> AddEntry(h3, "NA49");
   legend->SetTextSize(0.035);
   legend -> Draw();

   TPaveText *pt = new TPaveText(0.26,0.75,0.86,0.88,"NDC NB");
   pt->SetTextSize(0.032);
   pt->SetFillColor(0);
   pt->SetTextAlign(12);
   pt->AddText("^{208}Pb, 158A GeV, b #approx 2.1 fm ");
   std::string str = "#LTN_{p, h}#GT = " + std::to_string(p) + "; #LTN_{p, e}#GT = " + std::to_string(p1);
   const char *c =str.c_str(); 
   pt->SetTextFont(42);
   pt->AddText(c);
   //pt->Draw(); 



   c1->Update();

}
