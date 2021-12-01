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

  ofstream out("ALICE.txt");

 
// Make a canva and a pad
TCanvas* depo = new TCanvas("depo","NPP",1800,900);
gStyle->SetOptStat(0);
depo -> Divide(2,1,0.001,0.001);
depo -> cd(1);
  
  TH1D* histo3_p = new TH1D("histo3_p","", 20, 0, 20);
  TH1D* histo3_n = new TH1D("histo3_n","", 20, 0, 20);
  TH1D* histo4_p = new TH1D("histo4_p","", 20, 0, 20);
  TH1D* histo4_n = new TH1D("histo4_n","", 20, 0, 20);
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
Int_t P_bin = 0;
//Double_t bmin, bmax, per1=0.8, per2=1.0;
float b; 
Int_t Nvoid;
std::vector<double>* Z;
std::vector<double>* A;
float num_n_b[20] = {0};float num_p_b[20] = {0}; float num_n_at_b[20] = {0}; float num_p_at_b[20] = {0}; float num_void_at_b[20] = {0}; float num_react_at_b[20] = {0};


/*
TFile* ReadFile3 = new TFile("../pbpb5020GeV_100k_pT2x.root");

TTree *tree3 = (TTree*) ReadFile3->Get("Glauber"); 
tree3->SetBranchAddress("A_on_B", &A);
tree3->SetBranchAddress("Z_on_B", &Z);
tree3->SetBranchAddress("impact_parameter", &b);*/

std::vector<double> MassOnSideA;
std::vector<double> ChargeOnSideA;
float btree;
TFile* ReadFile3 = new TFile("aamcc.root","RECREATE","Demo ROOT file with histograms & trees");
TTree* tree31 = new TTree("Glauber","Events from glauber modeling");
cout<<"test0"<<endl;
tree31->Branch("A_on_A", "std::vector" ,&MassOnSideA);
tree31->Branch("Z_on_A", "std::vector" ,&ChargeOnSideA);
tree31->Branch("impact_parameter", &btree, "impact_parameter/f");
cout<<"test1"<<endl;

TChain* fChain=new TChain("events");
for (int i = 0; i <= 49; i++) {fChain->Add(Form("%s/dcmqgsm_%d.root", "/media/sf_Ubuntu_share/lhc_dcm", i));}
    //TFile* ReadFile = new TFile("../input/QA_dcmqgsm.root");

cout<<"Numer of entries equal "<<fChain->GetEntries()<<endl;
UEvent* fEvent = new UEvent;
EventInitialState* fIniState = new EventInitialState;
fChain->SetBranchAddress("event", &fEvent); 
fChain->SetBranchAddress("iniState", &fIniState); 
    
    Long64_t lNEvents = fChain->GetEntries();
    Long64_t fNpa;
    UParticle* fParticle;
                    
  for (long i = 0; i < lNEvents; i++)
  {
        fChain->GetEntry(i);
        //cout<<"test2"<<endl;
        fNpa = fEvent->GetNpa();
        btree = fEvent->GetB();
        for (int j=0;j<fNpa;j++)
        {
            //cout<<"test3"<<endl;
            fParticle = fEvent->GetParticle(j);
            if (fParticle -> T() == 1)
            {
              if (fParticle->GetPdg()>1e9) {
                MassOnSideA.push_back(fParticle->GetPdg()/10%1000);
                ChargeOnSideA.push_back(fParticle->GetPdg()/10000%1000);
              }
              else if (fParticle->GetPdg()==2212) {
                MassOnSideA.push_back(1.);
                ChargeOnSideA.push_back(1.);
              }
              else if (fParticle->GetPdg()==2112) {
                MassOnSideA.push_back(1.);
                ChargeOnSideA.push_back(0.);
              }
            }
        }
        //cout<<MassOnSideA<<endl;
        tree31->Fill();
        MassOnSideA.clear();
        ChargeOnSideA.clear();
        //cout<<"test4"<<endl;
    }
  ReadFile3->Write();
  ReadFile3->Close();
cout<<"test"<<endl;
TFile* ReadFile31 = new TFile("aamcc.root");
TTree *tree3 = (TTree*) ReadFile31->Get("Glauber"); 
tree3->SetBranchAddress("A_on_A", &A);
tree3->SetBranchAddress("Z_on_A", &Z);
tree3->SetBranchAddress("impact_parameter", &b);


  for(int k = 0; k < tree3->GetEntries(); k++){
    tree3->GetEntry(k);
    bhist->Fill(b);
    cout<<b<<endl;
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
    /*if(flag_n_b == 1)*/ num_n_b[int(b)]+=1;
    /*if(flag_p_b == 1)*/ num_p_b[int(b)]+=1;

  }
 

  for (int k = 0; k < 20; k++)
  {
     histo3_p->Fill(float(k)+0.5, (num_p_at_b[k]/num_p_b[k])*factor->GetBinContent(factor->FindBin(k+0.5)));
     histo3_n->Fill(float(k)+0.5, (num_n_at_b[k]/num_n_b[k])*factor->GetBinContent(factor->FindBin(k+0.5)));
  }
  
//num_n_b[20] = {0}; num_p_b[20] = {0}; num_n_at_b[20] = {0}; num_p_at_b[20] = {0}; num_void_at_b[20] = {0}; num_react_at_b[20] = {0};

 
/*TFile* ReadFile4 = new TFile("../pbpb5020GeV_100k_pT2x.root");

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
    num_n_b[int(b)]+=1;
    num_p_b[int(b)]+=1;

  }
 


  for (int k = 0; k < 20; k++)
  {
     histo4_p->Fill(float(k)+0.5, (num_p_at_b[k]/num_p_b[k])*factor->GetBinContent(factor->FindBin(k+0.5)));
     histo4_n->Fill(float(k)+0.5, (num_n_at_b[k]/num_n_b[k])*factor->GetBinContent(factor->FindBin(k+0.5)));
  }*/

     

  neutron_hist->SetStats(kFALSE);
  neutron_hist->SetYTitle("#LT N_{n} #GT");
  neutron_hist->SetXTitle("b, fm");
  neutron_hist->SetLabelSize(0.03);
  neutron_hist->SetTitleSize(0.035, "xy");  
  neutron_hist->SetMaximum(60);




//np_histo->SetLineStyle(10);
//np_histo->SetLineStyle(kSolid);
neutron_hist->SetLineWidth(2);
neutron_hist->SetMinimum(0.);
neutron_hist->GetXaxis()->SetLabelSize(0.03);
neutron_hist->GetYaxis()->SetLabelSize(0.03);
//np_histo->GetZaxis()->SetLabelSize(0.03);
neutron_hist->Draw("SAME HIST");

histo3_n->SetLineStyle(5);
histo3_n->SetLineWidth(2);
histo3_n->SetLineColor(30);
histo3_n->Draw("SAME HIST");

histo4_n->SetLineStyle(7);
histo4_n->SetLineWidth(2);
histo4_n->SetLineColor(kBlack);
//histo4_n->Draw("SAME HIST");

depo->Update();
   
auto legend = new TLegend(0.617,0.759,0.8759,0.8625);
   legend -> AddEntry(neutron_hist, "ALICE");
   legend -> AddEntry(histo3_n, "AAMCC, d(E*, A_{pf.})");
   //legend -> AddEntry(histo4_n, "AAMCC, 10 MeV");
   legend -> SetTextSize(0.024);
   legend -> Draw();

depo -> cd(2);

proton_hist->SetStats(kFALSE);
proton_hist->SetYTitle("#LT N_{p} #GT");
proton_hist->SetXTitle("b, fm");
proton_hist->SetLabelSize(0.03);
proton_hist->SetTitleSize(0.035, "xy"); 
proton_hist->SetMaximum(35);

proton_hist->SetLineWidth(2);
proton_hist->SetMinimum(0.);
proton_hist->GetXaxis()->SetLabelSize(0.03);
proton_hist->GetYaxis()->SetLabelSize(0.03);
//np_histo->GetZaxis()->SetLabelSize(0.03);
proton_hist->Draw("HIST");

histo3_p->SetLineStyle(5);
histo3_p->SetLineWidth(2);
histo3_p->SetLineColor(30);
histo3_p->Draw("SAME HIST");

histo4_p->SetLineStyle(7);
histo4_p->SetLineWidth(2);
histo4_p->SetLineColor(kBlack);
//histo4_p->Draw("SAME HIST");

depo->Update();

	TLatex* text = new TLatex(2,32,"PbPb, #sqrt{#it{s}_{NN}} = 5.02 TeV"); 
	text->SetTextFont(42);
	text->Draw();

auto legend2 = new TLegend(0.617,0.759,0.8759,0.8625);
   legend2 -> AddEntry(proton_hist, "ALICE, hot H-2");
   legend2 -> AddEntry(histo3_p, "AAMCC, cold H-2");
   legend2 -> SetTextSize(0.024);
   //legend2 -> Draw();

depo->Update();

}
