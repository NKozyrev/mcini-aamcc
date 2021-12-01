{
#include "Riostream.h"
#include "string.h"
#include "TLatex.h"
#include <stdio.h> 
#include <time.h> 

// https://link.springer.com/content/pdf/10.1007/s002180050403.pdf EMU-01/12
time_t start, end;
time(&start);

gROOT->Reset();
gROOT->ProcessLine("#include <vector>");


// Styles "Plain", "Bold", "Video", "Pub" are available
gROOT->SetStyle("Bold");
gStyle->Reset("Plain");
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
gStyle->SetPalette(1);
gStyle->SetCanvasColor(10);
gStyle->SetCanvasBorderMode(0);
gStyle->SetFrameLineWidth(1);
gStyle->SetFrameFillColor(kWhite);
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gStyle->SetPadBottomMargin(0.15);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetHistLineColor(kRed);
gStyle->SetFuncWidth(2);
gStyle->SetFuncColor(kGreen);
gStyle->SetLabelSize(0.045,"xyz");
gStyle->SetLabelOffset(0.005,"y");
gStyle->SetLabelOffset(0.005,"x");
gStyle->SetLabelColor(kBlack,"xyz");
gStyle->SetTitleSize(0.07,"xyz");
gStyle->SetTitleOffset(1.1,"y");
gStyle->SetTitleOffset(1.0,"x");
gStyle->SetTextSizePixels(26);
gStyle->SetTextFont(42);
gStyle->SetPalette(kRainBow);
gStyle->SetLineWidth(1);
gStyle->SetCanvasColor(0);
gStyle->SetTitleFillColor(0);
gStyle->SetTitleBorderSize(0);
gStyle->SetStatColor(0);
gStyle->SetHistLineWidth(0);
gStyle->SetLineScalePS(0.3);
//  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 
gStyle->SetLegendFillColor(kWhite);
//  gStyle->SetFillColor(kWhite);
gStyle->SetLegendFont(42);

// Make a canva and a pad
TCanvas* depo = new TCanvas("depo","Comparison with EMU-01/12; Au on AgBr",0,0,1280,960);
TPad* pad = new TPad("pad","The pad",0.02,0.02,0.98,0.98,21);
pad->SetBorderMode(0); pad->SetFillColor(kWhite); pad->SetTickx(); pad->SetTicky(); pad->Draw();
//pad->SetGridx(1);
//pad->SetGridy(1);
pad->SetLogz();
pad->Divide(2,2, 1e-12,1e-12);

// First
double sigma = 0;
float b;
std::vector<double>* Z;
std::vector<double>* A;
std::vector<double>* pXA;
std::vector<double>* pYA;
std::vector<double>* pZA;

// Second
double sigma_2 = 0;
float b_2;
std::vector<double>* Z_2;
std::vector<double>* A_2;
std::vector<double>* pXA2;
std::vector<double>* pYA2;
std::vector<double>* pZA2;

TFile* ReadFile = new TFile("../au2ag10_7AGeV_20k_Emax8.root");
TTree *tree = (TTree*) ReadFile->Get("Glauber");
tree->SetBranchAddress("impact_parameter", &b);
tree->SetBranchAddress("A_on_A", &A);
tree->SetBranchAddress("pX_on_A", &pXA);
tree->SetBranchAddress("pY_on_A", &pYA);
tree->SetBranchAddress("pZ_on_A", &pZA);
tree->SetBranchAddress("Z_on_A", &Z);

TFile* ReadFile2 = new TFile("../au2br10_7AGeV_20k_Emax8.root");
TTree *tree2 = (TTree*) ReadFile2->Get("Glauber");
tree2->SetBranchAddress("impact_parameter", &b_2);
tree2->SetBranchAddress("A_on_A", &A_2);
tree2->SetBranchAddress("pX_on_A", &pXA2);
tree2->SetBranchAddress("pY_on_A", &pYA2);
tree2->SetBranchAddress("pZ_on_A", &pZA2);
tree2->SetBranchAddress("Z_on_A", &Z_2);

TTree *init_par = (TTree*) ReadFile->Get("Conditions");
init_par->SetBranchAddress("Xsect_total", &sigma);
init_par->GetEntry(0);

TTree *init_par_2 = (TTree*) ReadFile2->Get("Conditions");
init_par_2->SetBranchAddress("Xsect_total", &sigma_2);
init_par_2->GetEntry(0);

double arr_zb[20];
double arr_zb_nz[40];
double arr_zmax[20];
double arr_zmax_2[20];
double arr_err_zmax[20];
double arr_err_zmax_2[20];
double arr_imf[20];
double arr_imf_2[20];
double arr_err_imf[20];
double arr_err_imf_2[20];
double arr_nz[40];
double arr_nz_2[40];
double arr_err_nz1[40];
double arr_err_nz1_2[40];
double arr_nz2[40];
double arr_nz2_2[40];
double arr_err_nz2[40];
double arr_err_nz2_2[40];
double arr_err_x[20];
double arr_err_x_nz[40];

// Empty arrays
for (int j = 0; j<40; j++){
  arr_zb_nz[j] = 0;
  arr_nz[j] = 0;
  arr_nz2[j] = 0;
  arr_err_nz1[j] = 0;
  arr_err_nz1_2[j] = 0;
  arr_nz_2[j] = 0;
  arr_nz2_2[j] = 0;
  arr_err_nz2[j] = 0;
  arr_err_nz2_2[j] = 0;
  arr_err_x_nz[j] = 0;
  if(j < 20){
    arr_zb[j] = 0;
    arr_zmax[j] = 0;
    arr_zmax_2[j] = 0;
    arr_err_zmax[j] = 0;
    arr_err_zmax_2[j] = 0;
    arr_imf[j] = 0;
    arr_imf_2[j] = 0;
    arr_err_imf[j] = 0;
    arr_err_imf_2[j] = 0;
    arr_err_x[j] = 0;
  }
}

double counter_z_max = 0;
double counter_imf = 0;
double counter_nz = 0;
double counter_errors = 0;
double zmax = 0;
double z_b3 = 0;
double z_b2 = 0;
double m_imf = 0;
double z_imf = 0;
double n_z = 0;
double n_z2 = 0;

std::cout << "Start of the first file" << endl;

for (int j = 0; j<40; j++){
  counter_z_max = 0;
  counter_imf = 0;
  counter_nz = 0;
  for(int k = 0; k < tree->GetEntries(); k++){
    tree->GetEntry(k);
    z_b3 = 0;
    z_b2 = 0;
    zmax = 0;
    z_imf = 0;
    n_z = 0;
    n_z2 = 0;
    for(int i = 0; i<Z->size(); i++){
      // Z_b3
      if (Z->at(i) >= 3){
        z_b3 +=Z->at(i);
      }
      // Z_bound
      if (Z->at(i) >= 2){
        z_b2 +=Z->at(i);
      }
      // Z_MAX
      if (Z->at(i) > zmax ){
        zmax = Z->at(i);
      }
      // M_IMF
      if (Z->at(i) >= 3 && Z->at(i) <= 30){
        z_imf += 1;
      }
      // N_{Z=1}
      if (Z->at(i) == 1 && (TMath::ATan(pow(pow(pXA->at(i), 2) + pow(pYA->at(i), 2), 0.5)/pZA->at(i)) < 17.2/1000)){
        n_z += 1;
      }
      // N_{Z=2}
      if (Z->at(i) == 2){
        n_z2 += 1;
      }
    }
    // Z_MAX
    if (z_b3 > 0 + 4*j && z_b3 <= 4 + 4*j && j < 20){
      arr_zmax[j] += zmax;
      arr_err_zmax[j] += pow(zmax, 2);
      counter_z_max += 1;
    }
    // M_IMF
    if (z_b2 > 0 + 4*j && z_b2 <= 4 + 4*j && j < 20){
      arr_imf[j] += z_imf;
      arr_err_imf[j] += pow(z_imf, 2);
      counter_imf += 1;
    }
    // N_{Z=1} and N_{Z=2} 
    if (z_b3 > 2 + 2*j && z_b3 <= 4 + 2*j){
      arr_nz[j] += n_z;
      arr_err_nz1[j] += pow(n_z, 2);
      arr_nz2[j] += n_z2;
      arr_err_nz2[j] += pow(n_z2, 2);
      counter_nz += 1;
    }
  }
  if (j < 20){
    arr_zb[j] = 2 + 4*j; // For Z_MAX and IMF
  }
  arr_zb_nz[j] = 3 + 2*j; // For NZ_1 and NZ_2 
  // Z_MAX     
  if (j < 20){
    if(counter_z_max != 0){
      arr_zmax[j] = arr_zmax[j]/counter_z_max;
      arr_err_zmax[j] = pow(arr_err_zmax[j]/counter_z_max - pow(arr_zmax[j], 2), 0.5)/pow(counter_z_max - 1, 0.5);
    }
    // M_IMF
    if(counter_imf != 0){
      arr_imf[j] = arr_imf[j]/counter_imf;
      arr_err_imf[j] = pow(arr_err_imf[j]/counter_imf - pow(arr_imf[j], 2), 0.5)/pow(counter_imf - 1, 0.5);
    }
  }
  // N_{Z=1} and N_{Z=2} 
  if(counter_nz != 0){
    arr_nz[j] = arr_nz[j]/counter_nz;
    arr_err_nz1[j] = pow(arr_err_nz1[j]/counter_nz - pow(arr_nz[j], 2), 0.5)/pow(counter_nz - 1, 0.5);
    arr_nz2[j] = arr_nz2[j]/counter_nz;
    arr_err_nz2[j] = pow(arr_err_nz2[j]/counter_nz - pow(arr_nz2[j], 2), 0.5)/pow(counter_nz - 1, 0.5);
  }
}

std::cout << "End of the first file" << endl;

std::cout << "Start of the second file" << endl;

for (int j = 0; j<40; j++){
  counter_z_max = 0;
  counter_imf = 0;
  counter_nz = 0;
  for(int k = 0; k < tree2->GetEntries(); k++){
    tree2->GetEntry(k);
    z_b3 = 0;
    z_b2 = 0;
    zmax = 0;
    z_imf = 0;
    n_z = 0;
    n_z2 = 0;
    for(int i = 0; i<Z_2->size(); i++){
      // Z_b3
      if (Z_2->at(i) >= 3){
        z_b3 +=Z_2->at(i);
      }
      // Z_bound
      if (Z_2->at(i) >= 2){
        z_b2 +=Z_2->at(i);
      }
      // Z_MAX
      if (Z_2->at(i) > zmax ){
        zmax = Z_2->at(i);
      }
      // M_IMF
      if (Z_2->at(i) >= 3 && Z_2->at(i) <= 30){
        z_imf += 1;
      }
      // N_{Z=1}
      if (Z_2->at(i) == 1 && (TMath::ATan(pow(pow(pXA2->at(i), 2) + pow(pYA2->at(i), 2), 0.5)/pZA2->at(i)) < 17.2/1000)){
        n_z += 1;
      }
      // N_{Z=2}
      if (Z_2->at(i) == 2){
        n_z2 += 1;
      }
    }
    // Z_MAX
    if (z_b3 > 0 + 4*j && z_b3 <= 4 + 4*j && j < 20){
      arr_zmax_2[j] += zmax;
      arr_err_zmax_2[j] += pow(zmax, 2);
      counter_z_max += 1;
    }
    // M_IMF
    if (z_b2 > 0 + 4*j && z_b2 <= 4 + 4*j && j < 20){
      arr_imf_2[j] += z_imf;
      arr_err_imf_2[j] += pow(z_imf, 2);
      counter_imf += 1;
    }
    // N_{Z=1} and N_{Z=2} 
    if (z_b3 > 2 + 2*j && z_b3 <= 4 + 2*j){
      arr_nz_2[j] += n_z;
      arr_err_nz1_2[j] += pow(n_z, 2);
      arr_nz2_2[j] += n_z2;
      arr_err_nz2_2[j] += pow(n_z2, 2);
      counter_nz += 1;
    }
  }
  if (j < 20){
    arr_zb[j] = 2 + 4*j; // For Z_MAX and IMF
  }
  arr_zb_nz[j] = 3 + 2*j; // For NZ_1 and NZ_2 
  // Z_MAX     
  if (j < 20){
    if(counter_z_max != 0){
      arr_zmax_2[j] = arr_zmax_2[j]/counter_z_max;
      arr_err_zmax_2[j] = pow(arr_err_zmax_2[j]/counter_z_max - pow(arr_zmax_2[j], 2), 0.5)/pow(counter_z_max - 1, 0.5);
    }
    // M_IMF
    if(counter_imf != 0){
      arr_imf_2[j] = arr_imf_2[j]/counter_imf;
      arr_err_imf_2[j] = pow(arr_err_imf_2[j]/counter_imf - pow(arr_imf_2[j], 2), 0.5)/pow(counter_imf - 1, 0.5);
    }
    else
    {
      arr_imf_2[j] = 0;
    }
  }
  // N_{Z=1} and N_{Z=2} 
  if(counter_nz != 0){
    arr_nz_2[j] = arr_nz_2[j]/counter_nz;
    arr_err_nz1_2[j] = pow(arr_err_nz1_2[j]/counter_nz - pow(arr_nz_2[j], 2), 0.5)/pow(counter_nz - 1, 0.5);
    arr_nz2_2[j] = arr_nz2_2[j]/counter_nz;
    arr_err_nz2_2[j] = pow(arr_err_nz2_2[j]/counter_nz - pow(arr_nz2_2[j], 2), 0.5)/pow(counter_nz - 1, 0.5);
  }
}

std::cout << "End of the second file" << endl;

// Au + AgBr
for (int j = 0; j<40; j++){
  if (j < 20){
    arr_zmax[j] = (arr_zmax[j] + sigma_2/sigma*arr_zmax_2[j])/(1 + sigma_2/sigma);
    arr_imf[j] = (arr_imf[j] + sigma_2/sigma*arr_imf_2[j])/(1 + sigma_2/sigma);
  }
  arr_nz[j] = (arr_nz[j] + sigma_2/sigma*arr_nz_2[j])/(1 + sigma_2/sigma);
  arr_nz2[j] = (arr_nz2[j] + sigma_2/sigma*arr_nz2_2[j])/(1 + sigma_2/sigma);
}
for (int j = 0; j<40; j++){
  if (j < 20){
    arr_err_zmax[j] = pow((pow(arr_err_zmax[j],2) + pow((sigma_2/sigma),2)*pow(arr_err_zmax_2[j],2)),0.5)/(1 + sigma_2/sigma);
    arr_err_imf[j] = pow((pow(arr_err_imf[j],2) + pow((sigma_2/sigma),2)*pow(arr_err_imf_2[j],2)),0.5)/(1 + sigma_2/sigma);
  }
  arr_err_nz1[j] = pow((pow(arr_err_nz1[j],2) + pow((sigma_2/sigma),2)*pow(arr_err_nz1_2[j],2)),0.5)/(1 + sigma_2/sigma);
  arr_err_nz2[j] = pow((pow(arr_err_nz2[j],2) + pow((sigma_2/sigma),2)*pow(arr_err_nz2_2[j],2)),0.5)/(1 + sigma_2/sigma);
}

// Z_MAX ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Experimental data Adamovich
double x[37],y[37];
double err_x[37];
for (int j = 0; j<37; j++){
  err_x[j] = 0;
}
int i=0;
ifstream file("zmax_exp_adamovich.txt");
while(i<37)
{
  char c;
  file>>x[i]>>y[i];
  if(file.eof())
       break;
  //std::cout << x[i] << y[i] <<endl;
  i++;
}

// Experimental data Aladin
double x2[39],y2[39];
double err_x2[39];
for (int j = 0; j<39; j++){
  err_x2[j] = 0;
}
i=0;
ifstream file2("zmax_exp_aladin.txt");
while(i<39)
{
  char c;
  file2>>x2[i]>>y2[i];
  if(file2.eof())
       break;
  //std::cout << x2[i] << y2[i] <<endl;
  i++;
}

pad->cd(1);

TMultiGraph* mg = new TMultiGraph();  
mg->GetXaxis()->SetTitle("#it{z}_{b3}");
mg->GetYaxis()->SetTitle("<#it{Z}_{max}>");
mg->GetYaxis()->SetTitleOffset(1.15);

TGraphErrors* plot_AAMCC_1 = new TGraphErrors(20,arr_zb,arr_zmax, arr_err_x, arr_err_zmax);
TGraphErrors* plot_adamovich = new TGraphErrors(37,x,y,err_x,err_x);
TGraphErrors* plot_botvina = new TGraphErrors(37,x2,y2,err_x2,err_x2);

plot_adamovich->SetTitle("");
plot_adamovich->SetMarkerStyle(22);
plot_adamovich->SetMarkerColor(kCyan+2);
plot_adamovich->SetMarkerSize(2.2);

plot_botvina->SetTitle("");
plot_botvina->SetMarkerStyle(29);
plot_botvina->SetMarkerColor(9);
plot_botvina->SetMarkerSize(2.2);

plot_AAMCC_1->SetTitle("");
plot_AAMCC_1->SetMarkerStyle(kFullCircle);
plot_AAMCC_1->SetMarkerColor(kMagenta+1);
plot_AAMCC_1->GetXaxis()->SetTitle("");
plot_AAMCC_1->GetYaxis()->SetTitle(""); 
plot_AAMCC_1->SetMarkerSize(2.0);

mg->Add(plot_botvina);
mg->Add(plot_adamovich);
mg->Add(plot_AAMCC_1);
mg->Draw("AP");
gPad->Modified();
mg->GetXaxis()->SetLimits(0.0001,80);
mg->SetMinimum(0.);
mg->SetMaximum(80.);

auto legend = new TLegend(0.15069133,0.724623,0.553378,0.90100);
legend->AddEntry(plot_AAMCC_1,"AAMCC","p");
legend->AddEntry(plot_adamovich,"Adamovich et al. 97","p");
legend->AddEntry(plot_botvina,"Botvina et al. 95","p");
legend->SetTextSize(gStyle->GetTextSize()*1);
legend->Draw();

TLatex* sqrtS = new TLatex(4.24,52.7, "10A GeV Au on AgBr");
sqrtS->SetTextSize(0.05);
sqrtS->Draw();
pad->Update();

// IMF ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double x_imf[39],y_imf[39];
double err_x_imf[39], err_y_imf[39];
for (int j = 0; j<39; j++){
  err_x_imf[j] = 0;
  err_y_imf[j] = 0;
}
// Array of exp data
i=0;
ifstream file_imf("z_bound_imf_exp.txt");
while(i<39)
{
  char c;
  file_imf>>x_imf[i]>>y_imf[i];
  if(file_imf.eof())
       break;
  //std::cout << x_imf[i] << y_imf[i]<<endl;
  i++;
}
// errors
i=0;
ifstream file_imf_err("z_bound_imf_exp_errors.txt");
while(i<39)
{
  char c;
  file_imf_err >> err_y_imf[i];
  if(file_imf_err.eof())
       break;
  i++;
}
for (int i = 0; i<39; i++)
{
  err_y_imf[i] = err_y_imf[i] - y_imf[i];
}

pad->cd(2);

TMultiGraph* mg_imf = new TMultiGraph();
mg_imf->GetXaxis()->SetTitle("#it{z}_{bound}");
mg_imf->GetYaxis()->SetTitle("<#it{M}_{IMF}>");
mg_imf->GetYaxis()->SetTitleOffset(1.15);

TGraphErrors* plot_imf_1 = new TGraphErrors(20,arr_zb,arr_imf, arr_err_x, arr_err_imf);
TGraphErrors* plot_imf_exp = new TGraphErrors(39,x_imf,y_imf,err_x_imf,err_y_imf);

plot_imf_1->SetTitle("");
plot_imf_1->SetMarkerStyle(kFullCircle);
plot_imf_1->SetMarkerColor(kMagenta+1);
plot_imf_1->SetMarkerSize(2.1);

plot_imf_exp->SetTitle("");
plot_imf_exp->SetMarkerStyle(22);
plot_imf_exp->SetMarkerColor(kCyan+2);
plot_imf_exp->GetXaxis()->SetTitle("");
plot_imf_exp->GetYaxis()->SetTitle(""); 
plot_imf_exp->SetMarkerSize(2.0);
mg_imf->GetXaxis()->SetLimits(0.1,80);
mg_imf->SetMinimum(0);
//mg_imf->SetMaximum(4.);
mg_imf->Add(plot_imf_exp);
mg_imf->Add(plot_imf_1);
mg_imf->Draw("AP"); 
gPad->Modified();
pad->Update();

// N_{Z=1} ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double x_nz1[19],y_nz1[19];
double err_x_nz1[19], err_y_nz1[19];
for (int j = 0; j<19; j++){
  err_x_nz1[j] = 0;
}
i=0;
ifstream file_nz1("nz1_exp.txt");
while(i<19)
{
  char c;
  file_nz1>>x_nz1[i]>>y_nz1[i] >> err_y_nz1[i];
  if(file_nz1.eof())
       break;
  //std::cout << x_nz1[i] << y_nz1[i] << err_y_nz1[i] <<endl;
  i++;
}
for (int j = 0; j<19; j++){
  err_y_nz1[j] = y_nz1[j] - err_y_nz1[j];
}

pad->cd(3);

TMultiGraph* mg_nz1 = new TMultiGraph();
mg_nz1->GetXaxis()->SetTitle("#it{z}_{b3}");
mg_nz1->GetYaxis()->SetTitle("<#it{N}_{Z = 1}>");
mg_nz1->GetYaxis()->SetTitleOffset(1.15);

TGraphErrors* plot_nz1 = new TGraphErrors(40,arr_zb_nz,arr_nz, arr_err_x_nz, arr_err_nz1);
TGraphErrors* plot_nz1_exp = new TGraphErrors(19,x_nz1,y_nz1,err_x_nz1,err_y_nz1);

plot_nz1_exp->SetTitle("");
plot_nz1_exp->SetMarkerStyle(22);
plot_nz1_exp->SetMarkerColor(kCyan+2);
plot_nz1_exp->SetMarkerSize(2.1);

plot_nz1->SetTitle("");
plot_nz1->SetMarkerStyle(kFullCircle);
plot_nz1->SetMarkerColor(kMagenta+1);
plot_nz1->GetXaxis()->SetTitle("");
plot_nz1->GetYaxis()->SetTitle(""); 
plot_nz1->SetMarkerSize(2.0); 

mg_nz1->GetXaxis()->SetLimits(4,80);
mg_nz1->SetMinimum(0.);
mg_nz1->Add(plot_nz1_exp);
mg_nz1->Add(plot_nz1);
mg_nz1->Draw("AP");
gPad->Modified();
pad->Update();

// N_{Z=2} ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double x_nz2[39],y_nz2[39];
double err_x_nz2[39], err_y_nz2[39];
for (int j = 0; j<39; j++){
  err_x_nz2[j] = 0;
}
i=0;
ifstream file_nz2("nz2_exp.txt");
while(i<39)
{
  char c;
  file_nz2>>x_nz2[i]>>y_nz2[i] >> err_y_nz2[i];
  if(file_nz2.eof())
       break;
  //std::cout << x_nz2[i] << y_nz2[i] << err_y_nz2[i] <<endl;
  i++;
}
for (int j = 0; j<39; j++){
  err_y_nz2[j] = y_nz2[j] - err_y_nz2[j];
}

pad->cd(4);

TMultiGraph* mg_nz2 = new TMultiGraph();
mg_nz2->GetXaxis()->SetTitle("#it{z}_{b3}");
mg_nz2->GetYaxis()->SetTitle("<#it{N}_{Z = 2}>");
mg_nz2->GetYaxis()->SetTitleOffset(1.15);

TGraphErrors* plot_nz2 = new TGraphErrors(40,arr_zb_nz,arr_nz2, arr_err_x_nz, arr_err_nz2);
TGraphErrors* plot_nz2_exp = new TGraphErrors(39,x_nz2,y_nz2,err_x_nz2,err_y_nz2);

plot_nz2_exp->SetTitle("");
plot_nz2_exp->SetMarkerStyle(22);
plot_nz2_exp->SetMarkerColor(kCyan+2);
plot_nz2_exp->SetMarkerSize(2.1);

plot_nz2->SetTitle("");
plot_nz2->SetMarkerStyle(kFullCircle);
plot_nz2->SetMarkerColor(kMagenta+1);
plot_nz2->GetXaxis()->SetTitle("");
plot_nz2->GetYaxis()->SetTitle("");
plot_nz2->SetMarkerSize(2.0); 

mg_nz2->GetXaxis()->SetLimits(4,80);
mg_nz2->SetMinimum(0.);

mg_nz2->Add(plot_nz2_exp);
mg_nz2->Add(plot_nz2);
mg_nz2->Draw("AP");
gPad->Modified();
//mg_nz2->SetMaximum(9.);

pad->Update();

time(&end);
double seconds = difftime(end, start);
std::cout << "Time of calc in minutes: " << seconds/(60) << endl;

}
