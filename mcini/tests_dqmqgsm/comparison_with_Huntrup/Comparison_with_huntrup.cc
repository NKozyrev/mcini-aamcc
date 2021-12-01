{
#include "Riostream.h"
#include "string.h"
#include "TLatex.h"
#include <stdio.h> 
#include <time.h> 

// https://journals.aps.org/prc/pdf/10.1103/PhysRevC.61.034903 - HUNTRUP
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
gStyle->SetPadColor(10);
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gStyle->SetPadBottomMargin(0.15);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetHistLineColor(kRed);
gStyle->SetFuncWidth(2);
gStyle->SetFuncColor(kGreen);
gStyle->SetLineWidth(1);
gStyle->SetLabelSize(0.075,"xyz");
gStyle->SetLabelOffset(0.005,"y");
gStyle->SetLabelOffset(0.005,"x");
gStyle->SetLabelColor(kBlack,"xyz");
gStyle->SetTitleSize(0.09,"xyz");
gStyle->SetTitleOffset(0.7,"y");
gStyle->SetTitleOffset(0.8,"x");
gStyle->SetTitleFillColor(kWhite);
gStyle->SetTextSizePixels(26);
gStyle->SetTextFont(42);
gStyle->SetHistLineWidth(0);
//  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 
gStyle->SetLegendFillColor(kWhite);
//  gStyle->SetFillColor(kWhite);
gStyle->SetLegendFont(42);

// Make a canva and a pad
TCanvas* depo = new TCanvas("depo","Comparison with Huntrup",0,0,1280,960);
TPad* pad = new TPad("pad","The pad",0.02,0.02,0.98,0.98,21);
pad->SetBorderMode(0); pad->SetFillColor(kWhite); pad->SetTickx(); pad->SetTicky(); pad->Draw();
//pad->SetGridx(1);
//pad->SetGridy(1);
pad->SetLogz();
pad->Divide(2,3, 1e-12,1e-12);

double sigma = 5.43;
double norm = 100000;
float b;
double kinEn;
double z_b7 = 0;
double z_m1 = 0;
double z_m2 = 0;
double z_m3 = 0;
double m_imf = 0;
double m_0 = 0;
double m_1 = 0;
double m_2 = 0;
double z_mean = 0;
double x_er = 0;
double y_er = 0;
std::vector<double>* Z;
std::vector<double>* A;

TFile* ReadFile = new TFile("../pbpb158AGeV_20k_Emax8.root");
TTree *tree = (TTree*) ReadFile->Get("Glauber");
tree->SetBranchAddress("impact_parameter", &b);
tree->SetBranchAddress("A_on_A", &A);
tree->SetBranchAddress("Z_on_A", &Z);


TTree *init_par = (TTree*) ReadFile->Get("Conditions");

init_par->SetBranchAddress("Xsect_total", &sigma);
init_par->GetEntry(0);
sigma*=10e3;

double arr_zb[20];
double arr_a12[20];
double arr_a23[20];
double arr_imf[20];
double arr_zmax[20];
double arr_a123[20];
double arr_gamma[20];
double arr_err[20];
double arr_err_12[20];
double arr_err_imf[20];
double arr_err_123[20];
double arr_err_23[20];
double arr_err_gamma[20];
double arr_err_zb[20];
double counter = 0;
double counter_12 = 0;
double counter_23 = 0;
double counter_123 = 0;
double counter_gamma = 0;
double mx2_12 = 0;
double mx2_23 = 0;
double mx2_zmax = 0;
double mx2_gamma = 0;
double mx2_123 = 0;
double mx2_imf = 0;

for (int j = 0; j<20; j++){
  arr_zb[j] = 0;
  arr_a12[j] = 0;
  arr_a23[j] = 0;
  arr_imf[j] = 0;
  arr_a123[j] = 0;
  arr_gamma[j] = 0;
  arr_err[j] = 0;
  arr_err_zb[j] = 0;
  arr_err_gamma[j] = 0;
  arr_err_23[j] = 0;
  arr_err_123[j] = 0;
  arr_err_12[j] = 0;
  arr_zmax[j] = 0;
}

norm = tree->GetEntries();
for (int j = 0; j<20; j++){
  counter = 0;
  counter_123 = 0;
  counter_23 = 0;
  counter_12 = 0;
  counter_gamma = 0;
  mx2_12 = 0;
  mx2_23 = 0;
  mx2_123 = 0;
  mx2_zmax = 0;
  mx2_imf = 0;
  mx2_gamma = 0;
  for(int k = 0; k < tree->GetEntries(); k++){
    tree->GetEntry(k);
    z_b7 = 0;
    z_m1 = 0;
    z_m2 = 0;
    z_m3 = 0;
    m_1 = 0;
    m_2 = 0;
    m_0 = 0;
    m_imf = 0;
    z_mean = 0;
    for(int i = 0; i<Z->size(); i++){
      if (Z->at(i) >= 7){
        z_b7 +=Z->at(i);
      }
      if (Z->at(i) >= z_m1 && Z->at(i) > z_m2 && Z->at(i) > z_m3 && Z->at(i) >= 7){
        z_m3 = z_m2;
        z_m2 = z_m1;
        z_m1 = Z->at(i);
      }
      if (Z->at(i) < z_m1 && Z->at(i) >= z_m2 && Z->at(i) > z_m3 && Z->at(i) >= 7){
        z_m3 = z_m2;
        z_m2 = Z->at(i);
      }
      if (Z->at(i) < z_m1 && Z->at(i) < z_m2 && Z->at(i) >= z_m3 && Z->at(i) >= 7){
        z_m3 = Z->at(i);
      }
      if (Z->at(i) <= 30 && Z->at(i) >= 7){
        m_imf += 1;
      }
    }
    for(int i = 0; i<Z->size(); i++){
      if (Z->at(i) != z_m1 && Z->at(i) >= 7){
        m_0 += 1;
        m_1 += Z->at(i);
        m_2 += pow(Z->at(i),2);
      }
    }
    if (z_b7 > 0 + 4*j && z_b7 < 4 + 4*j){
      if (z_m1 != 0 && z_m2 != 0){
        arr_a12[j] += (z_m1 - z_m2)/(z_m1 + z_m2);
        mx2_12 += pow((z_m1 - z_m2)/(z_m1 + z_m2), 2);
        counter_12 += 1;
      }
      if (z_m2 != 0 && z_m3 != 0){
        arr_a23[j] += (z_m2 - z_m3)/(z_m2 + z_m3);
        mx2_23 += pow((z_m2 - z_m3)/(z_m2 + z_m3), 2);
        counter_23 += 1;
      }
      if (m_0 > 0 && m_2 != 0 && m_1 != 0){
          arr_gamma[j] += m_2*m_0/pow(m_1,2);
          mx2_gamma += pow(m_2*m_0/pow(m_1,2), 2);
          counter_gamma += 1;
      }
      if (z_m1 != 0 && z_m2 != 0 && z_m3 != 0){
        z_mean = (z_m1 + z_m2 + z_m3)/3;
        arr_a123[j] += (pow(pow(z_m1 - z_mean,2) + pow(z_m2 - z_mean,2) + pow(z_m3 - z_mean,2),0.5))/(pow(6, 0.5)*z_mean);
        mx2_123 += pow((pow(pow(z_m1 - z_mean,2) + pow(z_m2 - z_mean,2) + pow(z_m3 - z_mean,2),0.5))/(pow(6, 0.5)*z_mean), 2);
        counter_123 += 1;
        //std::cout << "z_m1: "<< z_m1 << "z_m2: " << z_m2 << "z_m3: " << z_m3 << "m_imf: " << m_imf <<  "z_b7: " << z_b7 << endl;
      }
      arr_zmax[j] += z_b7 - z_m1;
      mx2_zmax += pow((z_b7 - z_m1), 2);
      arr_imf[j] += m_imf;
      mx2_imf += pow(m_imf, 2);
      counter += 1;
      if (z_b7> 50 && z_b7 < 56){
        // std::cout << z_m1 << " " << z_m2 << " " << z_m3 << " " << m_imf <<endl;
      }
    }
  }
  arr_zb[j] = 3 + 4*j;
  if(counter_12 != 0){
    arr_a12[j] = arr_a12[j]/counter_12;
    arr_err_12[j] = pow(mx2_12/counter_12 - pow(arr_a12[j], 2), 0.5)/pow(counter_12 - 1, 0.5);
  }
  if(counter_23 != 0){
    arr_a23[j] = arr_a23[j]/counter_23;
    arr_err_23[j] = pow(mx2_23/counter_23 - pow(arr_a23[j], 2), 0.5)/pow(counter_23 - 1, 0.5);
  }
  if(counter != 0){
    arr_imf[j] = arr_imf[j]/counter;
    arr_err_imf[j] = pow(mx2_imf/counter - pow(arr_imf[j], 2), 0.5)/pow(counter - 1, 0.5);
  }
  if(counter != 0){
    arr_zmax[j] = arr_zmax[j]/counter;
    arr_err_zb[j] = pow(mx2_zmax/counter - pow(arr_zmax[j], 2), 0.5)/pow(counter - 1, 0.5);
  }
  if(counter_123 != 0){
    arr_a123[j] = arr_a123[j]/counter_123;
    arr_err_123[j] = pow(mx2_123/counter_123 - pow(arr_a123[j], 2), 0.5)/pow(counter_123 - 1, 0.5);
  }
  if(counter_gamma != 0){
    arr_gamma[j] = arr_gamma[j]/counter_gamma;
    arr_err_gamma[j] = pow(mx2_gamma/counter_gamma - pow(arr_gamma[j], 2), 0.5)/pow(counter_gamma - 1, 0.5);
  }
}

// Experimental data IMF
double x[11],y[11];
double err_x[11];
for (int j = 0; j<11; j++){
  err_x[j] = 0;
}
int i=0;
ifstream file("imf_exp_data.txt");
while(i<11)
{
  char c;
  file>>x[i]>>y[i];
  if(file.eof())
       break;
  //std::cout << x[i] << y[i] <<endl;
  i++;
}
// Experimental data a_12
double x2[11],y2[11];
double err_x2[11], err_y2[11];
for (int j = 0; j<11; j++){
  err_x2[j] = 0;
  err_y2[j] = 0;
}
i=0;
ifstream file2("a12_exp_data.txt");
while(i<11)
{
  char c;
  file2>>x2[i]>>y2[i];
  if(file2.eof())
       break;
  //std::cout << x2[i] << y2[i] <<endl;
  i++;
}
// errors
i=0;
ifstream file2_b("errors_a12.txt");
while(i<11)
{
  char c;
  file2_b>>x_er>>y_er;
  if(file2_b.eof())
       break;
  err_y2[i] = y_er - y2[i];
  i++;
}
// Experimental data a_23
double x3[9],y3[9];
double err_x3[9], err_y3[9];
for (int j = 0; j<9; j++){
  err_x3[j] = 0;
  err_y3[j] = 0;
}
i=0;
ifstream file3("a23_exp_data.txt");
while(i<9)
{
  char c;
  file3>>x3[i]>>y3[i];
  if(file3.eof())
       break;
  //std::cout << x3[i] << y3[i] <<endl;
  i++;
}
// errors
i=0;
ifstream file3_b("errors_a23.txt");
while(i<9)
{
  char c;
  file3_b>>x_er>>y_er;
  if(file3_b.eof())
       break;
  err_y3[i] = y_er - y3[i];
  i++;
}
// Experimental data a_123
double x4[8],y4[8];
double err_x4[8], err_y4[8];
for (int j = 0; j<8; j++){
  err_x4[j] = 0;
  err_y4[j] = 0;
}
i=0;
ifstream file4("a123_exp_data.txt");
while(i<8)
{
  char c;
  file4>>x4[i]>>y4[i];
  if(file4.eof())
       break;
  //std::cout << x4[i] << y4[i] <<endl;
  i++;
}
// errors
i=0;
ifstream file4_b("errors_a123.txt");
while(i<8)
{
  char c;
  file4_b>>x_er>>y_er;
  if(file4_b.eof())
       break;
  err_y4[i] = y_er - y4[i];
  i++;
}
// Experimental data gamma_2
double x5[10],y5[10];
double err_x5[10], err_y5[10];
for (int j = 0; j<10; j++){
  err_x5[j] = 0;
  err_y5[j] = 0;
}
i=0;
ifstream file5("gamma2_exp_data.txt");
while(i<10)
{
  char c;
  file5>>x5[i]>>y5[i];
  if(file5.eof())
       break;
  //std::cout << x5[i] << y5[i] <<endl;
  i++;
}
// errors
i=0;
ifstream file5_b("errors_gamma.txt");
while(i<10)
{
  char c;
  file5_b>>x_er>>y_er;
  if(file5_b.eof())
       break;
  err_y5[i] = y_er - y5[i];
  i++;
}
// Experimental data zb - zmax
double x6[11],y6[11];
double err_x6[11], err_y6[11];
for (int j = 0; j<11; j++){
  err_x6[j] = 0;
  err_y6[j] = 0;
}
i=0;
ifstream file6("zb_zmax_exp_data.txt");
while(i<11)
{
  char c;
  file6>>x6[i]>>y6[i];
  if(file6.eof())
       break;
  //std::cout << x6[i] << y6[i] <<endl;
  i++;
}
// errors
i=0;
ifstream file6_b("errors_zb.txt");
while(i<11)
{
  char c;
  file6_b>>x_er>>y_er;
  if(file6_b.eof())
       break;
  err_y6[i] = y_er - y6[i];
  i++;
}

pad->Update();
pad->cd(1);

TMultiGraph* mg_mimf = new TMultiGraph();
mg_mimf->GetXaxis()->SetTitle("#it{z}_{b7}");
mg_mimf->GetYaxis()->SetTitle("#it{M}_{imf}");
mg_mimf->GetYaxis()->SetTitleOffset(1.15);
mg_mimf->GetYaxis()->SetTitleSize(10);

TGraphErrors* plot_mimf = new TGraphErrors(20,arr_zb,arr_imf,arr_err,arr_err_imf);
TGraphErrors* plot_exp_imf = new TGraphErrors(11,x,y,err_x,err_x);

plot_exp_imf->SetTitle("");
plot_exp_imf->SetMarkerStyle(kFullSquare);
plot_exp_imf->SetMarkerColor(kCyan+2);
plot_exp_imf->SetMarkerSize(2.2);

plot_mimf->SetTitle("");
plot_mimf->SetMarkerStyle(kFullCircle);
plot_mimf->SetMarkerColor(kMagenta+1);
plot_mimf->SetMarkerSize(2.0);

mg_mimf->Add(plot_exp_imf);
mg_mimf->Add(plot_mimf);
mg_mimf->Draw("AP");
gPad->Modified();
mg_mimf->SetMinimum(0.001);

auto legend = new TLegend(0.15069133,0.155,0.532147,0.407371);
legend->AddEntry(plot_mimf,"AAMCC","p");
legend->AddEntry(plot_exp_imf,"Huntrup et al.","p");
legend->SetTextSize(gStyle->GetTextSize()*1.8);
legend->Draw();

TLatex* sqrtS = new TLatex(49.97466,2.024223, "158A GeV Pb on Pb");
sqrtS->SetTextSize(0.070067);
sqrtS->Draw();

pad->Update();
pad->cd(2);

TMultiGraph* mg_zmax = new TMultiGraph();
mg_zmax->GetXaxis()->SetTitle("#it{z}_{b7}");
mg_zmax->GetYaxis()->SetTitle("#it{z}_{b7} - <#it{Z}_{max}>");
mg_zmax->GetYaxis()->SetTitleOffset(1.15);

TGraphErrors* plot_zmax = new TGraphErrors(20,arr_zb,arr_zmax,arr_err,arr_err_zb);
TGraphErrors* plot_exp_zb_zmax = new TGraphErrors(11,x6,y6,err_x6,err_y6);

plot_exp_zb_zmax->SetTitle("");
plot_exp_zb_zmax->SetMarkerStyle(kFullSquare);
plot_exp_zb_zmax->SetMarkerColor(kCyan+2);
plot_exp_zb_zmax->SetMarkerSize(2.2);

plot_zmax->SetTitle("");
plot_zmax->SetMarkerStyle(kFullCircle);
plot_zmax->SetMarkerColor(kMagenta+1);
plot_zmax->SetMarkerSize(2.0);

mg_zmax->Add(plot_zmax, "Y+AP");
mg_zmax->Add(plot_exp_zb_zmax, "Y+AP");
mg_zmax->Draw("Y+AP");
gPad->Modified();
mg_zmax->SetMinimum(0.001);

pad->Update();
pad->cd(3);

TMultiGraph* mg = new TMultiGraph();
mg->GetXaxis()->SetTitle("#it{z}_{b7}");
mg->GetYaxis()->SetTitle("<#it{a}_{12}>");
mg->GetYaxis()->SetTitleOffset(1.15);

TGraphErrors* plot = new TGraphErrors(20,arr_zb,arr_a12,arr_err,arr_err_12);
TGraphErrors* plot_exp_a12 = new TGraphErrors(11,x2,y2,err_x2,err_y2);

plot_exp_a12->SetTitle("");
plot_exp_a12->SetMarkerStyle(kFullSquare);
plot_exp_a12->SetMarkerColor(kCyan+2);
plot_exp_a12->SetMarkerSize(2.2);

plot->SetTitle("");
plot->SetMarkerStyle(kFullCircle);
plot->SetMarkerColor(kMagenta+1);
plot->SetMarkerSize(2.0);

mg->Add(plot);
mg->Add(plot_exp_a12);
mg->Draw("AP");
gPad->Modified();
mg->SetMinimum(0.001);

pad->Update();
pad->cd(4);

TMultiGraph* mg_a23 = new TMultiGraph();
mg_a23->GetXaxis()->SetTitle("#it{z}_{b7}");
mg_a23->GetYaxis()->SetTitle("<#it{a}_{23}>");
mg_a23->GetYaxis()->SetTitleOffset(1.15);

TGraphErrors* plot_a23 = new TGraphErrors(20,arr_zb,arr_a23,arr_err,arr_err_23);
TGraphErrors* plot_exp_a23 = new TGraphErrors(9,x3,y3,err_x3,err_y3);

plot_exp_a23->SetTitle("");
plot_exp_a23->SetMarkerStyle(kFullSquare);
plot_exp_a23->SetMarkerColor(kCyan+2);
plot_exp_a23->SetMarkerSize(2.2);

plot_a23->SetTitle("");
plot_a23->SetMarkerStyle(kFullCircle);
plot_a23->SetMarkerColor(kMagenta+1);
plot_a23->SetMarkerSize(2.0);

mg_a23->Add(plot_a23, "Y+AP");
mg_a23->Add(plot_exp_a23, "Y+AP");
mg_a23->Draw("Y+AP");
gPad->Modified();
mg_a23->SetMinimum(0.001);

pad->Update();
pad->cd(5);

TMultiGraph* mg_a123 = new TMultiGraph();
mg_a123->GetXaxis()->SetTitle("#it{z}_{b7}");
mg_a123->GetYaxis()->SetTitle("<#it{a}_{123}>");
mg_a123->GetYaxis()->SetTitleOffset(1.15);

TGraphErrors* plot_a123 = new TGraphErrors(20,arr_zb,arr_a123,arr_err,arr_err_123);
TGraphErrors* plot_exp_a123 = new TGraphErrors(8,x4,y4,err_x4,err_y4);

plot_exp_a123->SetTitle("");
plot_exp_a123->SetMarkerStyle(kFullSquare);
plot_exp_a123->SetMarkerColor(kCyan+2);
plot_exp_a123->SetMarkerSize(2.2);

plot_a123->SetTitle("");
plot_a123->SetMarkerStyle(kFullCircle);
plot_a123->SetMarkerColor(kMagenta+1);
plot_a123->SetMarkerSize(2.0);

mg_a123->Add(plot_a123);
mg_a123->Add(plot_exp_a123);
mg_a123->Draw("AP");
gPad->Modified();
mg_a123->SetMinimum(0.001);

pad->Update();
pad->cd(6);

TMultiGraph* mg_gamma = new TMultiGraph();
mg_gamma->GetXaxis()->SetTitle("#it{z}_{b7}");
mg_gamma->GetYaxis()->SetTitle("#it{#gamma}_{2}");
mg_gamma->GetYaxis()->SetTitleOffset(0.85);

TGraphErrors* plot_gamma = new TGraphErrors(20,arr_zb,arr_gamma,arr_err,arr_err_gamma);
TGraphErrors* plot_exp_gamma = new TGraphErrors(10,x5,y5,err_x5,err_y5);

plot_exp_gamma->SetTitle("");
plot_exp_gamma->SetMarkerStyle(kFullSquare);
plot_exp_gamma->SetMarkerColor(kCyan+2);
plot_exp_gamma->SetMarkerSize(2.2);

plot_gamma->SetTitle("");
plot_gamma->SetMarkerStyle(kFullCircle);
plot_gamma->SetMarkerColor(kMagenta+1);
plot_gamma->SetMarkerSize(2.0);

mg_gamma->Add(plot_gamma);
mg_gamma->Add(plot_exp_gamma);
mg_gamma->Draw("APY+");
mg_gamma->GetXaxis()->SetLimits(0.01,80);
mg_gamma->SetMinimum(1.);
mg_gamma->SetMaximum(1.03);

gPad->Modified();
gPad->Update();
pad->Update();

time(&end);
double seconds = difftime(end, start);
std::cout << "Time of calc in minutes: " << seconds/(60) << endl;

}
