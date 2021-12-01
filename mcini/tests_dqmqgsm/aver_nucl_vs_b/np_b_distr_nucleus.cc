#include "Riostream.h"
#include "string.h"
#include "TLatex.h"
#include "TGaxis.h"
#include <fstream>
{
gROOT->Reset();
gROOT->ProcessLine("#include <vector>");
//TGaxis::SetMaxDigits(1); 

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
  gStyle->SetLineWidth(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(62);
  gStyle->SetLineScalePS(0.3);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);

  ofstream out("ALICE.txt");

 
// Make a canva and a pad
TCanvas* depo = new TCanvas("depo","NPP",1800,900);
depo->SetFillColor(0);
depo->SetBorderMode(0);
depo->SetBorderSize(2);
depo->SetFrameBorderMode(0);
depo->SetFrameBorderMode(0);
gStyle->SetOptStat(0);
depo -> Divide(2,1,0.001,0.001);
depo -> cd(1);
  
  TH1D* histo3_p = new TH1D("histo3_p","", 20, 0, 20);
  TH1D* histo3_n = new TH1D("histo3_n","", 20, 0, 20);
  TH1D* histo4_p = new TH1D("histo3_p","", 20, 0, 20);
  TH1D* histo4_n = new TH1D("histo3_n","", 20, 0, 20);
  TH1D* neutron_hist = new TH1D("neutron ","Neutron distribution", 20, 0, 20);
  TH1D* proton_hist = new TH1D("proton ","Proton distribution", 20, 0, 20);
  TH1D* bhist = new TH1D("bhist","",20,0,20);
  
  for (int k=0; k<20; k++){
    double x = double(k) + 0.5;
    if(9.016 - 0.887*x + 1.410*pow(x,2) - 0.118*pow(x,3) + 0.002*pow(x,4) >= 0) neutron_hist->SetBinContent(k+1, 9.016 - 0.887*x + 1.410*pow(x,2) - 0.118*pow(x,3) + 0.002*pow(x,4));
    else neutron_hist->SetBinContent(k+1, 0);
    if(1.24 - 0.756*x + 0.556*pow(x,2) - 0.055*pow(x,3) + 0.0014*pow(x,4) >= 0) proton_hist->SetBinContent(k+1, 1.24 - 0.756*x + 0.556*pow(x,2) - 0.055*pow(x,3) + 0.0014*pow(x,4)); 
    else proton_hist->SetBinContent(k+1, 0);
  }
  
 TFile corr("Pb_Pb_glau_react_fact.root");
 TH1D* factor;
 corr.GetObject("factor",factor);
  

//TH1D* h1 = new TH1D("h1", "h1", 180, 0, 18);

double kinEn, sigma; 
Double_t N,Nact=0;
Int_t Nb;
Double_t Npall = 0, Npiter=0;
Int_t P_bin = 0;
//Double_t bmin, bmax, per1=0.8, per2=1.0;
float b; 
Int_t Nvoid;
std::vector<double>* Z;
std::vector<double>* A;
float num_n_b[20] = {0};float num_p_b[20] = {0}; float num_n_at_b[20] = {0}; float num_p_at_b[20] = {0}; float num_void_at_b[20] = {0}; float num_react_at_b[20] = {0};



TFile* ReadFile3 = new TFile("dcmqgsm_in_aamcc_format_lhc.root"); //"auau12AGeV_20k_Gold.root"); // MST

TTree *tree3 = (TTree*) ReadFile3->Get("Glauber"); 
tree3->SetBranchAddress("A_on_A", &A);
tree3->SetBranchAddress("Z_on_A", &Z);
tree3->SetBranchAddress("impact_parameter", &b);




  for(int k = 0; k < tree3->GetEntries(); k++){
    tree3->GetEntry(k);
    bhist->Fill(b);
    for (int fragment = 0; fragment < (A->size()); fragment++){
      if (A->at(fragment) == 1 && Z->at(fragment) == 0){
        num_n_at_b[int(b)]+= 1;
      }
      if (A->at(fragment) == 1 && Z->at(fragment) == 1){
        num_p_at_b[int(b)]+= 1;
      }
      num_void_at_b[int(b)] += Nvoid;
      num_react_at_b[int(b)] += 1;
    }
     num_n_b[int(b)]+=1;
     num_p_b[int(b)]+=1;

  }
 //cout << Npall/100000<< endl;

  for (int k = 0; k < 20; k++)
  {
     histo3_p->Fill(float(k)+0.5, (num_p_at_b[k]/num_p_b[k])*factor->GetBinContent(factor->FindBin(k+0.5)));
     histo3_n->Fill(float(k)+0.5, (num_n_at_b[k]/num_n_b[k])*factor->GetBinContent(factor->FindBin(k+0.5)));
  }
  
cout<<histo3_n->GetBinContent(20)<<endl;

//num_n_b = {0}; num_p_b = {0}; num_n_at_b = {0}; num_p_at_b = {0}; num_void_at_b = {0}; num_react_at_b = {0};

 
/*TFile* ReadFile4 = new TFile("./check_mst/PbPb_5_02atev_nomst.root"); // NO MST

TTree *tree4 = (TTree*) ReadFile4->Get("Glauber"); 
tree4->SetBranchAddress("A_on_B", &A);
tree4->SetBranchAddress("Z_on_B", &Z);
tree4->SetBranchAddress("impact_parameter", &b);
 
   for(int k = 0; k < tree3->GetEntries(); k++){
    tree4->GetEntry(k);
    bhist->Fill(b);
    int flag_n_b = 0;int flag_p_b = 0;
    for (int fragment = 0; fragment < (A->size()); fragment++){
      if (A->at(fragment) == 1 && Z->at(fragment) == 0){
        num_n_at_b[int(b)]+= 1;
        flag_n_b = 1;
      }
      if (A->at(fragment) == 1 && Z->at(fragment) == 1){
        num_p_at_b[int(b)]+= 1;
        flag_p_b = 1;
      }
      num_void_at_b[int(b)] += Nvoid;
      num_react_at_b[int(b)] += 1;
    }
    if(flag_n_b == 1) num_n_b[int(b)]+=1;
    if(flag_p_b == 1) num_p_b[int(b)]+=1;

  }*/
 


  for (int k = 0; k < 20; k++)
  {
     histo4_p->Fill(float(k)+0.5, (num_p_at_b[k]/num_p_b[k]));//*factor->GetBinContent(factor->FindBin(k+0.5)));
     histo4_n->Fill(float(k)+0.5, (num_n_at_b[k]/num_n_b[k]));//*factor->GetBinContent(factor->FindBin(k+0.5)));
  }

     

  neutron_hist->SetStats(kFALSE);
  neutron_hist->SetYTitle("#LT N_{n} #GT");
  neutron_hist->SetXTitle("b, fm");
  neutron_hist->GetYaxis()->SetTitleFont(42);
  neutron_hist->GetXaxis()->SetTitleFont(42);
  neutron_hist->SetTitleSize(0.06, "xy");  
  neutron_hist->SetMaximum(50);
  neutron_hist->SetLineWidth(1);
  neutron_hist->SetMinimum(0.);
  neutron_hist->GetXaxis()->SetLabelSize(0.045);
  neutron_hist->GetYaxis()->SetLabelSize(0.045);
//np_histo->GetZaxis()->SetLabelSize(0.03);
  neutron_hist->Draw("HIST");


  histo3_n->GetXaxis()->SetLabelSize(0.045);
  histo3_n->GetYaxis()->SetLabelSize(0.045);
  histo3_n->SetYTitle("#LT N_{n} #GT");
  histo3_n->SetXTitle("b, fm");
  histo3_n->GetYaxis()->SetTitleFont(42);
  histo3_n->GetXaxis()->SetTitleFont(42);
  histo3_n->SetTitleSize(0.06, "xy");
histo3_n->SetLineStyle(5);
histo3_n->SetLineWidth(2);
histo3_n->SetLineColor(30);
histo3_n->GetXaxis()->SetRange(1,17);
histo3_n->Draw("SAME HIST");

histo4_n->SetLineStyle(7);
histo4_n->SetLineWidth(1);
histo4_n->SetLineColor(kBlack);
//histo4_n->Draw("SAME HIST");

depo->Update();
   
auto legend = new TLegend(0.617,0.759,0.8759,0.8625);
   legend -> AddEntry(neutron_hist, "ALICE");
   legend -> AddEntry(histo3_n, "DCM-QGSM");
  // legend -> AddEntry(histo4_n, "AAMCC NO MST");
   legend -> SetTextSize(0.045);
   legend -> Draw();

depo -> cd(2);

proton_hist->SetStats(kFALSE);
proton_hist->SetYTitle("#LT N_{p} #GT");
proton_hist->SetXTitle("b, fm");
proton_hist->SetLabelSize(0.045);
proton_hist->SetTitleSize(0.06, "xy"); 
proton_hist->SetMaximum(15);

proton_hist->SetLineWidth(2);
proton_hist->SetMinimum(0.);
proton_hist->GetXaxis()->SetLabelSize(0.045);
proton_hist->GetYaxis()->SetLabelSize(0.045);
proton_hist->GetYaxis()->SetTitleFont(42);
proton_hist->GetXaxis()->SetTitleFont(42);
//np_histo->GetZaxis()->SetLabelSize(0.03);
proton_hist->Draw("HIST");


histo3_p->GetXaxis()->SetLabelSize(0.045);
  histo3_p->GetYaxis()->SetLabelSize(0.045);
  histo3_p->SetYTitle("#LT N_{p} #GT");
  histo3_p->SetXTitle("b, fm");
  histo3_p->GetYaxis()->SetTitleFont(42);
  histo3_p->GetXaxis()->SetTitleFont(42);
  histo3_p->SetTitleSize(0.06, "xy");
histo3_p->SetStats(kFALSE);
histo3_p->SetYTitle("#LT N_{p} #GT");
histo3_p->SetXTitle("b, fm");
histo3_p->SetLabelSize(0.045);
histo3_p->SetTitleSize(0.06, "xy"); 


histo3_p->SetLineStyle(5);
histo3_p->SetLineWidth(2);
histo3_p->SetLineColor(30);
histo3_p->GetXaxis()->SetRange(1,17);
histo3_p->Draw("SAME HIST");

histo4_p->SetLineStyle(7);
histo4_p->SetLineWidth(2);
histo4_p->SetLineColor(kBlack);
//histo4_p->Draw("SAME HIST");

depo->Update();

	TLatex* text = new TLatex(2,32,"PbPb, #sqrt{#it{s}_{NN}} = 5.02 TeV"); 
	text->SetTextFont(42);
  text->SetTextSize(0.045);
	text->Draw();

  TLatex* text2 = new TLatex(2,29.56,"Data: ALICE Collaboration, ALICE-PUBLIC-2020-001"); 
  text2->SetTextFont(42);
  text2->SetTextSize(0.045);
//  text2->Draw();

auto legend2 = new TLegend(0.617,0.759,0.8759,0.8625);
legend2 -> SetTextSize(0.08);
   legend2 -> AddEntry(proton_hist, "ALICE");
   legend2 -> AddEntry(histo3_p, "AAMCC MST");
   //legend2 -> AddEntry(histo4_p, "AAMCC NO MST");
   
   //legend2 -> Draw();

depo->Update();

}
